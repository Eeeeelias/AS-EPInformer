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
                 device='cuda', useBN=True, usePromoterSignal=True, useFeat=True, n_extraFeat=0, useLN=True, 
                 exon_data=False, epigen_bp=False, separate_attention=False, use_histones=False, name_add=None, 
                 junctions=False):
        super(EPInformer_v2, self).__init__()
        self.n_enhancer = n_enhancer
        self.out_dim = out_dim
        self.useFeat = useFeat
        self.usePromoterSignal = usePromoterSignal
        self.n_extraFeat = n_extraFeat
        self.useBN = useBN
        self.base_size = base_size
        self.useLN = useLN
        self.junctions = junctions
        self.n_encoder = n_encoder
        self.device = device
        self.use_exon_data = exon_data
        self.use_epigen_bp = epigen_bp
        self.separate_attention = separate_attention
        self.histone_dim = 5 if use_histones else 0

        # 1. pre-trained sequence encoder
        test_name = f"-{name_add}" if name_add is not None else ''
        if pre_trained_encoder is not None:
            self.seq_encoder = pre_trained_encoder
            self.name = f'EPInformerV2{test_name}.preTrainedConv.{base_size}base.{out_dim}dim.{n_encoder}Trans.' \
                        f'{head}head.{useBN}BN.{useLN}LN.{useFeat}Feat.{n_extraFeat}extraFeat.{n_enhancer}enh.' \
                        f'{exon_data}exon.{separate_attention}exonAttn.{use_histones}histones'
        else:
            self.seq_encoder = seq_256bp_encoder(base_size=base_size)
            self.name = f'EPInformerV2{test_name}.{base_size}base.{out_dim}dim.{n_encoder}Trans.{head}head.{useBN}BN.' \
                        f'{useLN}LN.{useFeat}Feat.{n_extraFeat}extraFeat.{n_enhancer}enh.{exon_data}exon.' \
                        f'{separate_attention}exonAttn.{use_histones}histones'
        
        if self.use_exon_data:
            self.event_encoder = seq_256bp_encoder_small(base_size=4)

        if self.use_epigen_bp:
            self.epigen_encoder = seq_256bp_encoder_small(base_size=7)

        attn_encs = self.n_encoder if not self.use_epigen_bp else self.n_encoder+2
        if useLN: # use layer norm
            self.attn_encoder = get_clones(MHAttention_encoderLayer(d_model=out_dim, nhead=head), attn_encs)
        else:
            self.attn_encoder = get_clones(MHAttention_encoderLayer_noLN(d_model=out_dim, nhead=head), attn_encs)

        # get attention masks for promoter and exon data
        add_dim, n_last = (4, 3) if not self.junctions else (3, 2)

        self.attn_mask = self.promoter_attention_mask(additional_dim=add_dim if self.use_exon_data else 1)
        if self.separate_attention:
            self.exon_attn_mask = self.exon_attention_mask(additional_dim=add_dim, n_attend_last=n_last)
        else:
            self.exon_attn_mask = self.attn_mask

        if self.use_epigen_bp:
            self.epigen_mask_1 = self.epigen_attention_mask(is_one=True)
            self.epigen_mask_2 = self.epigen_attention_mask(is_one=False)

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
                nn.Linear(1 if self.junctions else 40, int(self.out_dim/32)), # added 40 to account for the event length
                nn.ELU(),
            )
            self.epigen_conv_out = nn.Sequential(
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
                nn.Linear(1 if self.junctions else 40, int(self.out_dim/32)), # added 40 to account for the event length
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
            nn.Linear(self.out_dim+feat_n + 256*2, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )

        if self.use_epigen_bp:
            self.pToSplice = nn.Sequential(
                nn.Linear(self.out_dim+feat_n+256*2, 128), # usually: 64 for the first token only
                nn.ReLU(),
                nn.Linear(128, 64),
                nn.ReLU(),
                nn.Linear(64, 1),
                # nn.Sigmoid()  
            )
        else:
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
        # assert additional_dim >= 4, "additional_dim must be at least 4 for exon attention mask"

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
    
    def epigen_attention_mask(self, is_one=False):
        attn_mask = (~np.identity(4).astype(bool))
        if is_one:
            attn_mask[:, 0] = False
            attn_mask[0, :] = False
            attn_mask[:, 2] = False
            attn_mask[2, :] = False
        else:
            attn_mask[:, 1] = False
            attn_mask[1, :] = False
            attn_mask[:, 3] = False
            attn_mask[3, :] = False
        attn_mask = torch.from_numpy(attn_mask)
        attn_mask.masked_fill(attn_mask, float('-inf'))
        return attn_mask

    def forward(self, pe_seq, exon_seq, hist_feat=None, epigen_feat=None, extraFeat=None, cell_line=None):
        # if enhancers_padding_mask is None:
        enhancers_padding_mask = ~(pe_seq.sum(-1).sum(-1) > 0).bool()
        enhancers_padding_mask[:, 0] = False # keep this only for ablation study where pe_seq is zeroed out
        if self.use_exon_data:
            exon_padding_mask = ~(exon_seq.sum(-1).sum(-1) > 0).bool() # added padding mask for event and combined
            enhancers_padding_mask = torch.concat([enhancers_padding_mask, exon_padding_mask], dim=1)
            # do not pad the last three sequences (intron, exon, intron) in the exon_seq
            enhancers_enhanced_padding_mask = enhancers_padding_mask.clone()
            #enhancers_enhanced_padding_mask[:, -3:] = False
        
        epigen_padding_mask = None
        epigen_flatten = None
        if epigen_feat is not None and self.use_epigen_bp:
            epigen_padding_mask = ~(epigen_feat.sum(-1).sum(-1) > 0).bool()
            epigen_padding_mask = torch.concat([epigen_padding_mask, exon_padding_mask], dim=1)
            epigen_embed = self.epigen_encoder(epigen_feat)
            epigen_embed = self.epigen_conv_out(epigen_embed)
            epigen_flatten = torch.flatten(epigen_embed.permute(0, 2, 1, 3), start_dim=2)

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
                                                                    attn_mask=self.exon_attn_mask.to(self.device))
            attn_list.append(attn.unsqueeze(0))

        if self.use_epigen_bp and epigen_feat is not None:
            # take last two entries from shape [16, 63, 64] -> [16, 2, 64]
            junction_embed = pe_flatten_embed[:, -2:, :]

            epigen_j1 = torch.cat([epigen_flatten, junction_embed], dim=1) # [16, 4, 64]
            epigen_j2 = torch.cat([epigen_flatten, junction_embed], dim=1)
            epigen_flat_emb_j1, attn = self.attn_encoder[-2](epigen_j1, enhancers_padding_mask=epigen_padding_mask,
                                                          attn_mask=self.epigen_mask_1.to(self.device))
            epigen_flat_emb_j2, attn = self.attn_encoder[-1](epigen_j2, enhancers_padding_mask=epigen_padding_mask,
                                                            attn_mask=self.epigen_mask_2.to(self.device))
            
            #epigen_flat_emb_j1 = torch.flatten(epigen_flat_emb_j1[:, 0, :], start_dim=1)
            #epigen_flat_emb_j2 = torch.flatten(epigen_flat_emb_j2[:, 0, :], start_dim=1)
            epigen_flat_emb_j1 = torch.flatten(epigen_flat_emb_j1, start_dim=1)
            epigen_flat_emb_j2 = torch.flatten(epigen_flat_emb_j2, start_dim=1)


        p_embed_expr = torch.flatten(pe_flatten_embed_expr[:,0,:], start_dim=1)
        p_embed_splice = torch.flatten(pe_flatten_embed_splice[:,0,:], start_dim=1)

        # concat with RNA data (+ promoter)
        if self.useFeat:
            # hijack the first three rna_feat with the cell_line info TODO: change - we can't hijack these anymore
            if cell_line is not None:
                pass
                # rna_feat = torch.cat([cell_line, rna_feat[..., 3:]], dim=-1)
            
            p_embed_expr = torch.cat([p_embed_expr, hist_feat], dim=-1) 
            p_embed_splice = torch.cat([p_embed_splice, hist_feat], dim=-1)

        if self.use_epigen_bp:
            p_embed_splice = torch.cat([p_embed_splice, epigen_flat_emb_j1, epigen_flat_emb_j2], dim=-1)

        # prediction heads
        p_expr = self.pToExpr(p_embed_expr)

        if self.use_exon_data:
            p_splice_binary_logits = self.pToSpliceBinary(p_embed_splice)
            p_splice_regression = self.pToSplice(p_embed_splice)
        else:
            p_splice_binary_logits = None
            p_splice_regression = None

        return p_expr, p_splice_binary_logits, p_splice_regression, torch.cat(attn_list)

### BELOW ONLY USED FOR TESTING BINARY PREDICTION

class ASInformer(nn.Module):
    def __init__(self):
        super().__init__()
        self.name = "ASInformer"

        self.event_encoder = seq_256bp_encoder_small(base_size=11)

        self.channel_red = nn.Sequential(
            nn.Conv2d(128, 64, kernel_size=1),
            nn.ReLU(),
            nn.Conv2d(64, 32, kernel_size=1),
            nn.ReLU(),
            nn.AdaptiveAvgPool2d((3, 8))
        )

        self.pred_head = nn.Sequential(
            nn.Linear(256*3, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, 1)
        )

    def forward(self, x):
        x = self.event_encoder(x)
        x = self.channel_red(x)
        x = torch.flatten(x.permute(0, 2, 1, 3), start_dim=2)
        x = x.reshape(x.shape[0], 256*3) 
        x = self.pred_head(x)
        return x


class ASTrInformer(nn.Module):
    def __init__(self, d_model=256, nhead=2, num_layers=2, dim_feedforward=32, dropout=0.3):
        super().__init__()
        self.name = "ASTrInformer"

        self.event_encoder_1 = seq_256bp_encoder_small(base_size=11)
        self.event_encoder_2 = seq_256bp_encoder_small(base_size=11)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, 
            nhead=nhead, 
            dim_feedforward=dim_feedforward, 
            dropout=dropout,
            batch_first=True
        )
        self.combination = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        self.final = nn.Sequential(
            nn.Conv1d(25, 16, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),  # global pooling
            nn.Flatten(),
            nn.Linear(16, 1),              # binary output
        )

    def forward(self, x):
        x1, x2 = torch.split(x, 1, dim=1)
        x1 = self.event_encoder_1(x1)
        x2 = self.event_encoder_2(x2)

        x1 = x1.permute(0, 2, 1, 3) # torch.Size([16, 1, 128, 25])
        x2 = x2.permute(0, 2, 1, 3) # torch.Size([16, 1, 128, 25])


        # combine
        x1 = x1.squeeze(1)
        x2 = x2.squeeze(1)
        combined = torch.cat([x1, x2], dim=1) # torch.Size([16, 256, 25])
        combined = combined.permute(0, 2, 1)

        # Pass through the combination transformer
        out = self.combination(combined) # torch.Size([16, 25, 256])

        # out = out.permute(0, 2, 1)  # torch.Size([16, 256, 25])

        # Final prediction
        return self.final(out)


class ASTransformer(nn.Module):
    def __init__(self, feature_dim=10, d_model=16, nhead=2, num_layers=2, dim_feedforward=32, dropout=0.3):
        super().__init__()
        self.name = "ASTransformer"
        
        # Project input features to model dimension
        self.input_proj = nn.Linear(feature_dim, d_model)
        
        # Shared Transformer Encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, 
            nhead=nhead, 
            dim_feedforward=dim_feedforward, 
            dropout=dropout,
            batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # Final classifier
        self.fc = nn.Linear(d_model, 1)
    
    def forward(self, seq):
        seq1 = seq[:, 0, :]
        seq2 = seq[:, 1, :]

        # Project features
        seq1 = self.input_proj(seq1)
        seq2 = self.input_proj(seq2)
        
        # Pass through shared transformer
        emb1 = self.transformer(seq1)[:, 0, :]  # take CLS-like first token
        emb2 = self.transformer(seq2)[:, 0, :]

        # Combine embeddings
        combined = emb1 + emb2  # or torch.cat([...], dim=-1) + linear

        # Output probability
        return self.fc(combined)


class ASLSTM(nn.Module):
    def __init__(self, hidden_size=2, dropout=0.3):
        super().__init__()
        self.name = "ASLSTM"
        
        # LSTMs for each sequence
        self.lstm_intron_exon = nn.LSTM(
            input_size=10,
            hidden_size=hidden_size,
            batch_first=True,
            bidirectional=False
        )
        
        self.lstm_exon_intron = nn.LSTM(
            input_size=10,
            hidden_size=hidden_size,
            batch_first=True,
            bidirectional=False
        )
        
        # After concatenation → feature dim = hidden_size * 2
        self.lstm_merged = nn.LSTM(
            input_size=hidden_size * 2,
            hidden_size=hidden_size * 2,
            batch_first=True,
            bidirectional=False
        )
        
        self.dropout = nn.Dropout(dropout)
        
        self.fc = nn.Linear(hidden_size * 2, 1)
        
    def forward(self, seq):
        # intron_exon, exon_intron: [batch, 240, 10]
        intron_exon = seq[:, 0, :]
        exon_intron = seq[:, 1, :]
        
        intron_exon_out, _ = self.lstm_intron_exon(intron_exon)  # [batch, 240, hidden]
        exon_intron_out, _ = self.lstm_exon_intron(exon_intron)  # [batch, 240, hidden]
        
        merged = torch.cat([intron_exon_out, exon_intron_out], dim=2)  # concat along features → [batch, 240, hidden*2]
        
        merged_out, (h_n, _) = self.lstm_merged(merged)  # merged_out: [batch, 240, hidden]
        
        # Keras `return_sequences=False` means take the last time step's hidden state
        last_hidden = merged_out[:, -1, :]  # [batch, hidden]
        
        dropped = self.dropout(last_hidden)
        
        out = self.fc(dropped)  # [batch, 1]
        
        return out


class TemporalEncoder(nn.Module):
    """
    1D temporal encoder over time axis.
    Input per stream: (B, 10, 240)  -> Output embedding: (B, emb_dim)
    """
    def __init__(self, in_ch=10, emb_dim=64, p_drop=0.2):
        super().__init__()
        self.conv1 = nn.Conv1d(in_ch, 64, kernel_size=5, padding=2)
        self.bn1   = nn.BatchNorm1d(64)
        self.conv2 = nn.Conv1d(64, 64, kernel_size=5, dilation=2, padding=4)
        self.bn2   = nn.BatchNorm1d(64)
        self.conv3 = nn.Conv1d(64, emb_dim, kernel_size=5, dilation=4, padding=8)
        self.bn3   = nn.BatchNorm1d(emb_dim)
        self.drop  = nn.Dropout(p_drop)
        # global average pooling over time at the end

    def forward(self, x):            # x: (B, 10, 240)
        x = self.drop(F.relu(self.bn1(self.conv1(x))))
        x = self.drop(F.relu(self.bn2(self.conv2(x))))
        x = F.relu(self.bn3(self.conv3(x)))
        x = F.adaptive_avg_pool1d(x, 1).squeeze(-1)  # (B, emb_dim)
        return x


class ASdCNNsmall(nn.Module):
    def __init__(self, emb_dim=64, shared_encoders=True, fusion="mlp", p_drop=0.2):
        super().__init__()
        self.name = "ASdCNNsmall"
        self.shared = shared_encoders
        self.encA = TemporalEncoder(in_ch=10, emb_dim=emb_dim, p_drop=p_drop)
        self.encB = self.encA if shared_encoders else TemporalEncoder(in_ch=10, emb_dim=emb_dim, p_drop=p_drop)

        # multiple fusion choices; default: concat + MLP
        self.fusion = fusion
        if fusion == "mlp":
            in_dim = emb_dim * 2
        elif fusion == "tricks":  # concat + (diff, prod)
            in_dim = emb_dim * 4
        elif fusion == "gated":   # learn a gate to weight A vs B
            in_dim = emb_dim * 3
            self.gate = nn.Sequential(nn.Linear(emb_dim*2, emb_dim), nn.ReLU(), nn.Linear(emb_dim, 1))
        else:
            raise ValueError("fusion must be one of: 'mlp', 'tricks', 'gated'")

        self.cls = nn.Sequential(
            nn.Dropout(p_drop),
            nn.Linear(in_dim, 64),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, x):            # x: (B, 2, 240, 10)
        b, two, t, f = x.shape
        assert two == 2 and f == 10, f"expected (B,2,240,10), got {x.shape}"
        a = x[:, 0].permute(0, 2, 1)  # (B, 240, 10) -> (B, 10, 240)
        b_ = x[:, 1].permute(0, 2, 1)

        za = self.encA(a)             # (B, emb_dim)
        zb = self.encB(b_)

        if self.fusion == "mlp":
            z = torch.cat([za, zb], dim=1)
        elif self.fusion == "tricks":
            z = torch.cat([za, zb, torch.abs(za - zb), za * zb], dim=1)
        elif self.fusion == "gated":
            gate = torch.sigmoid(self.gate(torch.cat([za, zb], dim=1)))  # (B,1)
            z = torch.cat([gate * za + (1 - gate) * zb, torch.cat([za, zb], dim=1)], dim=1)  # emb + context
        else:
            raise ValueError("fusion must be one of: 'mlp', 'tricks', 'gated'")
        out = self.cls(z)
        return out


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
    