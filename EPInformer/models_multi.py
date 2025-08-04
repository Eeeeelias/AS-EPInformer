import copy
import warnings
import numpy as np
import pandas as pd
from scipy import stats

from torch import Tensor
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data_utils
from torch.utils.data import Dataset, DataLoader

warnings.filterwarnings('ignore')

def get_clones(module, N):
    return nn.ModuleList([copy.deepcopy(module) for i in range(N)])

class seq_256bp_encoder(nn.Module):
    def __init__(self, base_size=4, out_dim=128, conv_dim=256):
        super(seq_256bp_encoder, self).__init__()
        self.conv_dim = conv_dim
        self.out_dim = out_dim
        self.base_size = base_size
        # cropped_len = 46
        self.stem_conv = nn.Sequential(
            nn.Conv2d(in_channels = base_size, out_channels = self.conv_dim, kernel_size = (1, 8), stride = 1, padding='same'),
            nn.ELU(),
        )
        self.conv_tower = nn.ModuleList([])
        conv_dim = [self.conv_dim, 128, 64, 64, 128]
        for i in range(4):
            self.conv_tower.append(nn.Sequential(
                nn.Conv2d(in_channels = conv_dim[i], out_channels=conv_dim[i+1], kernel_size=(1, 3), padding=(0, 1)),
                nn.BatchNorm2d(conv_dim[i+1]),
                nn.ELU(),                   
                nn.MaxPool2d(kernel_size=(1, 2), stride=(1, 2)),
            ))
            self.conv_tower.append(nn.Sequential(
                nn.Conv2d(in_channels = conv_dim[i+1], out_channels=conv_dim[i+1], kernel_size=(1, 1)),
                nn.ELU(),
            ))
        
    def forward(self, enhancers_input):
        if enhancers_input.shape[2] == 1:
            x_enhancer = enhancers_input
        else:
            x_enhancer = enhancers_input.permute(0, 3, 1, 2).contiguous()  
        x_enhancer = self.stem_conv(x_enhancer)
#         print(x_enhancer.shape)
        for i in range(0, len(self.conv_tower), 2):
            x_enhancer = self.conv_tower[i](x_enhancer)
            x_enhancer = self.conv_tower[i+1](x_enhancer) + x_enhancer
        return x_enhancer
    

class seq_256bp_encoder_small(nn.Module):
    def __init__(self, base_size=4, out_dim=128, conv_dim=256):
        super(seq_256bp_encoder_small, self).__init__()
        self.conv_dim = conv_dim
        self.out_dim = out_dim
        self.base_size = base_size
        self.stem_conv = nn.Sequential(
            nn.Conv2d(in_channels = base_size, out_channels = self.conv_dim, kernel_size = (1, 8), stride = 1, padding='same'),
            nn.ELU(),
        )
        self.conv_tower = nn.ModuleList([])
        conv_dim = [self.conv_dim, 64, 128]
        for i in range(2):
            conv = nn.Sequential(
                nn.Conv2d(in_channels = conv_dim[i], out_channels=conv_dim[i+1], kernel_size=(1, 3), padding=(0, 1)),
                nn.BatchNorm2d(conv_dim[i+1]),
                nn.ELU(),                   
            )
            conv_by_one = nn.Sequential(
                nn.Conv2d(in_channels = conv_dim[i+1], out_channels=conv_dim[i+1], kernel_size=(1, 1)),
                nn.ELU(),
            )
            pool = nn.MaxPool2d(kernel_size=(1, 4), stride=(1, 4))
            self.conv_tower.append(nn.ModuleList([conv, conv_by_one, pool]))
        
    def forward(self, enhancers_input):
        if enhancers_input.shape[2] == 1:
            x_enhancer = enhancers_input
        else:
            x_enhancer = enhancers_input.permute(0, 3, 1, 2).contiguous()  
        x_enhancer = self.stem_conv(x_enhancer)
        for conv, conv_by_one, pool in self.conv_tower: # type: ignore
            x_enhancer = conv(x_enhancer)
            x_enhancer = conv_by_one(x_enhancer)
            x_enhancer = pool(x_enhancer)
        return x_enhancer
    

class enhancer_predictor_256bp(nn.Module):
    def __init__(self):
        super(enhancer_predictor_256bp, self).__init__()
        self.encoder = seq_256bp_encoder()
        self.embedToAct = nn.Sequential(
            nn.Flatten(start_dim=1),
            nn.Linear(128*16, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(256, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(256, 1),
        )  
    def forward(self, enhancer_seq):
        if len(enhancer_seq.shape) < 4:
            enhancer_seq = enhancer_seq.unsqueeze(2)
        seq_embed = self.encoder(enhancer_seq)
        epi_out = self.embedToAct(seq_embed)
        return epi_out.squeeze(-1)


class MHAttention_encoderLayer(nn.Module):
    def __init__(self, d_model=128, nhead=8, dropout=0.):
        super(MHAttention_encoderLayer, self).__init__()
        # self.activation = activation
        self.self_attn = nn.MultiheadAttention(d_model, nhead, dropout=dropout, batch_first=True)
        
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        # Implementation of Feedforward model
        # self.linear1 = nn.Linear(d_model, 4*d_model) might cause loading problem, this parameter is not neccessary
        # self.linear2 = nn.Linear(4*d_model, d_model) might cause loading problem, this parameter is not neccessary
        # self.dropout = nn.Dropout(dropout)
        
        self.ff = nn.Sequential(
            nn.Linear(d_model, d_model*4),
            nn.ReLU(),
            nn.Linear(d_model*4, d_model)
        )
    # self-attention block
    def _sa_block(self, x, key_padding_mask, attn_mask):
        x, w = self.self_attn(x, x, x, key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        return x, w
         
    def forward(self, x, enhancers_padding_mask=None, attn_mask=None):
        x2 = self.norm1(x)
        x2, attention_w = self._sa_block(x2, key_padding_mask=enhancers_padding_mask, attn_mask=attn_mask)
        x = x2 + x
        x2 = self.norm2(x)
        x = x + self.ff(x2)
        return x, attention_w


class MHAttention_encoderLayer_noLN(nn.Module):
    def __init__(self, d_model=2048, nhead=8, dim_feedforward=256, dropout=0.1, activation=F.relu):
        super(MHAttention_encoderLayer_noLN, self).__init__()
        self.activation = activation
        self.self_attn = nn.MultiheadAttention(d_model, nhead, dropout=dropout, batch_first=True)
        
        # Implementation of Feedforward model
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.dropout = nn.Dropout(dropout)
        self.linear2 = nn.Linear(dim_feedforward, d_model)
    
    # feed forward block
    def _ff_block(self, x):
        x = self.linear2(self.dropout(self.activation(self.linear1(x))))
        return self.dropout(x)
    
    # self-attention block
    def _sa_block(self, x, key_padding_mask, attn_mask):
        x, w = self.self_attn(x, x, x,
                           key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        return x, w
    
    def forward(self, x_pe, enhancers_padding_mask=None, attn_mask=None):
        xt, attention_w = self._sa_block(x_pe, enhancers_padding_mask, attn_mask=attn_mask)
        x_pe = x_pe + xt
        x_pe = x_pe + self._ff_block(x_pe)
        return x_pe, attention_w


class EPInformer_v2(nn.Module):
    def __init__(self, base_size = 4, n_encoder=3, out_dim=128, head = 4, pre_trained_encoder= None, n_enhancer=50, 
                 device='cuda', useBN=True, usePromoterSignal=True, useFeat=True, n_extraFeat=0, useLN=True, exon_data=False,
                 separate_attention=False, use_histones=False):
        super(EPInformer_v2, self).__init__()
        self.n_enhancer = n_enhancer
        self.out_dim = out_dim
        self.useFeat = useFeat
        self.usePromoterSignal = usePromoterSignal
        self.n_extraFeat = n_extraFeat
        self.useBN = useBN
        self.base_size = base_size
        self.useLN = useLN
        self.n_encoder = n_encoder
        self.device = device
        self.use_exon_data = exon_data
        self.separate_attention = separate_attention
        self.histone_dim = 5 if use_histones else 0

        # 1. pre-trained sequence encoder
        if pre_trained_encoder is not None:
            self.seq_encoder = pre_trained_encoder
            self.name = f'EPInformerV2.preTrainedConv.{base_size}base.{out_dim}dim.{n_encoder}Trans.{head}head.{useBN}BN.' \
                        f'{useLN}LN.{useFeat}Feat.{n_extraFeat}extraFeat.{n_enhancer}enh.{exon_data}exon.' \
                        f'{separate_attention}exonAttn.{use_histones}histones'
        else:
            self.seq_encoder = seq_256bp_encoder(base_size=base_size)
            self.name = f'EPInformerV2.{base_size}base.{out_dim}dim.{n_encoder}Trans.{head}head.{useBN}BN.{useLN}LN.' \
                        f'{useFeat}Feat.{n_extraFeat}extraFeat.{n_enhancer}enh.{exon_data}exon.' \
                        f'{separate_attention}exonAttn.{use_histones}histones'
        
        if self.use_exon_data:
            self.event_encoder = seq_256bp_encoder_small(base_size=base_size)

        if useLN: # use layer norm
            self.attn_encoder = get_clones(MHAttention_encoderLayer(d_model=out_dim, nhead=head), self.n_encoder)
        else:
            self.attn_encoder = get_clones(MHAttention_encoderLayer_noLN(d_model=out_dim, nhead=head), self.n_encoder)
       
        # get attention masks for promoter and exon data
        self.attn_mask = self.promoter_attention_mask(additional_dim= 4 if self.use_exon_data else 1)
        if self.separate_attention:
            self.exon_attn_mask = self.exon_attention_mask(additional_dim=4, n_attend_last=3)
        else:
            self.exon_attn_mask = self.attn_mask

        if self.useBN: # use batch norm
            self.conv_out = nn.Sequential(
                nn.Conv2d(in_channels = 128, out_channels=64, kernel_size=(1, 3), dilation=(1, 2)),
                nn.BatchNorm2d(64),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 4)),
                nn.BatchNorm2d(64),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 6)),
                nn.BatchNorm2d(64),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=32, kernel_size=(1, 1)),
                nn.BatchNorm2d(32),
                nn.ELU(),
                nn.Linear(101, int(self.out_dim/32)), 
                 # nn.Linear(38, 8), # 2kb nn.Linear(101, 8)
                nn.ELU(),
            )
            self.event_conv_out = nn.Sequential(
                nn.Conv2d(in_channels = 128, out_channels=64, kernel_size=(1, 3), dilation=(1, 2)),
                nn.BatchNorm2d(64),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 4)),
                nn.BatchNorm2d(64),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 6)),
                nn.BatchNorm2d(64),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=32, kernel_size=(1, 1)),
                nn.BatchNorm2d(32),
                nn.ELU(),
                nn.Linear(40, int(self.out_dim/32)), # added 40 to account for the event length
                nn.ELU(),
            )
        else:
            self.conv_out = nn.Sequential(
                nn.Conv2d(in_channels = 128, out_channels=64, kernel_size=(1, 3), dilation=(1, 2)),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 4)),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 6)),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=32, kernel_size=(1, 1)),
                nn.ELU(),
                nn.Linear(101, int(self.out_dim/32)),
                # nn.Linear(38, 8), # 2kb nn.Linear(101, 8)
                nn.ELU(),
            )
            self.event_conv_out = nn.Sequential(
                nn.Conv2d(in_channels = 128, out_channels=64, kernel_size=(1, 3), dilation=(1, 2)),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 4)),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=64, kernel_size=(1, 3), dilation=(1, 6)),
                nn.ELU(),
                nn.Conv2d(in_channels = 64, out_channels=32, kernel_size=(1, 1)),
                nn.ELU(),
                nn.Linear(40, int(self.out_dim/32)), # added 40 to account for the event length
                nn.ELU(),
            )
        
        if self.useFeat:
            feat_n = 9 if self.usePromoterSignal else 8
        else:
            feat_n = 0

        self.pToExpr = nn.Sequential(
            nn.Linear(self.out_dim+feat_n, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
        )

        self.pToSpliceBinary = nn.Sequential(
            nn.Linear(self.out_dim+feat_n, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            nn.ReLU(),
        )

        self.pToSplice = nn.Sequential(
            nn.Linear(self.out_dim+feat_n, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            # nn.Sigmoid()  
        )
        self.add_pos_conv = nn.Sequential(
                nn.Conv1d(in_channels = self.out_dim + n_extraFeat + self.histone_dim, out_channels=self.out_dim, kernel_size=1),
                nn.ReLU(),
                nn.Conv1d(in_channels = self.out_dim, out_channels=self.out_dim, kernel_size=1),
                nn.ReLU(),
        )

    def shared_parameters(self):
        """
        Returns a list of parameters that are shared across the model.
        This is useful for optimizers that need to know which parameters to update.
        """
        return list(self.attn_encoder.parameters()) + list(self.add_pos_conv.parameters())

    def promoter_attention_mask(self, additional_dim=1):
        # creates attention mask that allows only the first token to attend and be attended to (so promoter sequence)
        attn_mask_additional_dim = additional_dim
        attn_mask = (~np.identity(self.n_enhancer+attn_mask_additional_dim).astype(bool))
        attn_mask[:, 0] = False
        attn_mask[0, :] = False
        attn_mask = torch.from_numpy(attn_mask)
        attn_mask.masked_fill(attn_mask, float('-inf'))
        return attn_mask

    def exon_attention_mask(self, additional_dim=4, n_attend_last=3):
        # attend to the last three tokens (intron, exon, intron sequence) and nothing else
        assert additional_dim >= 4, "additional_dim must be at least 4 for exon attention mask"

        attn_mask = (~np.identity(self.n_enhancer+additional_dim).astype(bool))
        # attend promoter
        attn_mask[:, 0] = False
        attn_mask[0, :] = False
        # attend last three tokens (intron/exon sequences)
        attn_mask[:, -n_attend_last:] = False
        attn_mask[-n_attend_last:, :] = False
        attn_mask = torch.from_numpy(attn_mask)
        attn_mask.masked_fill(attn_mask, float('-inf'))
        return attn_mask


    def forward(self, pe_seq, exon_seq, rna_feat=None, extraFeat=None, cell_line=None):
        # if enhancers_padding_mask is None:
        enhancers_padding_mask = ~(pe_seq.sum(-1).sum(-1) > 0).bool()
        enhancers_padding_mask[:, 0] = False # keep this only for ablation study where pe_seq is zeroed out
        if self.use_exon_data:
            exon_padding_mask = ~(exon_seq.sum(-1).sum(-1) > 0).bool() # added padding mask for event and combined
            enhancers_padding_mask = torch.concat([enhancers_padding_mask, exon_padding_mask], dim=1)
            # do not pad the last three sequences (intron, exon, intron) in the exon_seq
            enhancers_enhanced_padding_mask = enhancers_padding_mask.clone()
            #enhancers_enhanced_padding_mask[:, -3:] = False
    
        # encoding layer
        pe_embed = self.seq_encoder(pe_seq)
        pe_embed = self.conv_out(pe_embed)
        pe_flatten_embed = torch.flatten(pe_embed.permute(0, 2, 1, 3), start_dim=2)

        if self.use_exon_data:
            exon_embed = self.event_encoder(exon_seq)
            exon_embed = self.event_conv_out(exon_embed)        
            exon_flatten_embed = torch.flatten(exon_embed.permute(0, 2, 1, 3), start_dim=2)
            pe_flatten_embed = torch.concat([pe_flatten_embed, exon_flatten_embed], dim=1)

        if extraFeat is not None:
            # fill extraFeat with zeros on dim 1 because we have no promoter signal there
            if self.use_exon_data and extraFeat.shape[1] != pe_flatten_embed.shape[1]:
                extraFeat = F.pad(extraFeat, pad=(0,0,0,3))

            pe_flatten_embed = self.add_pos_conv(torch.concat([pe_flatten_embed, extraFeat], axis=-1)
                                                 .permute(0,2,1)).permute(0,2,1)

        # attention layers
        attn_list = []
        # split self.n_encode in half (e.g. if n_encoder=6 then first 3 layers are for enhancers and last 3 layers for event seqs)
        n_encoder_half = self.n_encoder // 2
        pe_flatten_embed_expr = pe_flatten_embed.clone()
        pe_flatten_embed_splice = pe_flatten_embed.clone()
        for i in range(n_encoder_half):
            pe_flatten_embed_expr, attn = self.attn_encoder[i](pe_flatten_embed_expr, enhancers_padding_mask=enhancers_padding_mask, 
                                                                attn_mask=self.attn_mask.to(self.device))
            attn_list.append(attn.unsqueeze(0))
            neg_i = self.n_encoder - i - 1
            pe_flatten_embed_splice, attn = self.attn_encoder[neg_i](pe_flatten_embed_splice, enhancers_padding_mask=enhancers_padding_mask, 
                                                                    attn_mask=self.attn_mask.to(self.device))
            attn_list.append(attn.unsqueeze(0))

        p_embed_expr = torch.flatten(pe_flatten_embed_expr[:,0,:], start_dim=1)
        p_embed_splice = torch.flatten(pe_flatten_embed_splice[:,0,:], start_dim=1)

        # concat with RNA data (+ promoter)
        if self.useFeat:
            # hijack the first three rna_feat with the cell_line info TODO: change - we can't hijack these anymore
            if cell_line is not None:
                pass
                # rna_feat = torch.cat([cell_line, rna_feat[..., 3:]], dim=-1)
                
            p_embed_expr = torch.cat([p_embed_expr, rna_feat], dim=-1) 
            p_embed_splice = torch.cat([p_embed_splice, rna_feat], dim=-1) 


        # prediction heads
        p_expr = self.pToExpr(p_embed_expr)

        if self.use_exon_data:
            p_splice_binary_logits = self.pToSpliceBinary(p_embed_splice)
            p_splice_regression = self.pToSplice(p_embed_splice)
        else:
            p_splice_binary_logits = None
            p_splice_regression = None

        return p_expr, p_splice_binary_logits, p_splice_regression, torch.cat(attn_list)


class WeightedLoss(nn.Module):
    """
    Class to compute a weighted loss for the EPInformer model that combines expression, splice binary, 
    and splice regression losses. The weights for each loss component are learned parameters.
    """
    def __init__(self):
        super().__init__()
        self.log_weight_expr = nn.Parameter(torch.tensor(1.0), requires_grad=True)
        self.log_weight_splice_binary = nn.Parameter(torch.tensor(1.0), requires_grad=True)
        self.log_weight_splice_reg = nn.Parameter(torch.tensor(1.0), requires_grad=True)

    def forward(self, loss_expr, loss_splice_binary, loss_splice_regression):
        weight_expr = torch.exp(-self.log_weight_expr)
        weight_splice_binary = torch.exp(-self.log_weight_splice_binary)
        weight_splice_regression = torch.exp(-self.log_weight_splice_reg)

        weighted_loss = (
            torch.exp(-self.log_weight_expr) * loss_expr + self.log_weight_expr +
            torch.exp(-self.log_weight_splice_binary) * loss_splice_binary + self.log_weight_splice_binary +
            torch.exp(-self.log_weight_splice_reg) * loss_splice_regression + self.log_weight_splice_reg
        )
        
        return weighted_loss, weight_expr, weight_splice_binary, weight_splice_regression
    