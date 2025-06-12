import torch
import os
import numpy as np
import pandas as pd
# torch
import torch.nn as nn
# import torch.nn.functional as F
import torch.optim
import torch.utils.data as data_utils
from torch.utils.data import Subset, Dataset

# from model.EPInformer import EPInformer_v2, enhancer_predictor_256bp
import h5py


class promoter_enhancer_dataset(Dataset):
    def __init__(self, data_folder = 'data/', expr_type='CAGE', usePromoterSignal=True, first_signal='distance', signal_type='H3K27ac', 
                 cell_type='K562', distance_threshold=None, hic_threshold=None, n_enhancers=50, n_extraFeat=1,
                 rna_seq_source='xpresso', tpm='gene', single_event_train=False):
        self.expr_type = expr_type
        self.cell_type = cell_type
        self.data_folder = data_folder
        self.first_signal = first_signal
        self.n_enhancers = n_enhancers
        self.signal_type = signal_type
        self.n_extraFeat = n_extraFeat
        self.usePromoterSignal = usePromoterSignal
        self.distance_threshold = distance_threshold
        self.hic_threshold = hic_threshold
        self.rna_seq_source = rna_seq_source
        self.use_normalized_psi = False
        self.tpm_level = "_summed_tpm" if tpm == 'transcript' else "_gene_level_tpm"
        self.gene_sequences = h5py.File(self.data_folder + '/event_sequences.h5', 'r')
        self.event_keys = list(self.gene_sequences.keys())
        self.psi_response = pd.read_csv(self.data_folder + '/psi_response.csv', index_col=0)
        if cell_type == 'K562':
            # promoter_df = pd.read_csv('/content/drive/MyDrive/EPInformer/EPInformer_activity/data/K562/DNase_ENCFF257HEE_Neighborhoods/GeneList.txt', sep='\t', index_col='symbol')
            promoter_df = pd.read_csv(self.data_folder + '/K562_DNase_ENCFF257HEE_hic_4DNFITUOMFUQ_1MB_ABC_nominated/DNase_ENCFF257HEE_Neighborhoods/GeneList.txt', sep='\t', index_col='symbol')
            promoter_df['PromoterActivity'] = np.sqrt(promoter_df['H3K27ac.RPM.TSS1Kb']*promoter_df['DHS.RPM.TSS1Kb'])
            self.promoter_df = promoter_df
            self.data_h5 = h5py.File(self.data_folder + '/K562_DNase_ENCFF257HEE_2kb_4DNFITUOMFUQ_enhancer_promoter_encoding.h5', 'r')
            # self.data_h5 = h5py.File('/content/drive/MyDrive/EPInformer/EPInformer_activity/data/K562/K562_DNase_ENCFF257HEE_2kb_noCutOff_hic_noFlankSeq_150kb60e_AllPutative_signals_False_v2.h5')
        elif cell_type == 'GM12878':
            promoter_df = pd.read_csv(self.data_folder + '/GM12878_DNase_ENCFF020WZB_hic_4DNFI1UEG1HD_1MB_ABC_nominated/DNase_ENCFF020WZB_Neighborhoods/GeneList.txt', sep='\t', index_col='symbol')
            promoter_df['PromoterActivity'] = np.sqrt(promoter_df['H3K27ac.RPM.TSS1Kb']*promoter_df['DHS.RPM.TSS1Kb'])
            self.promoter_df = promoter_df 
            self.data_h5 = h5py.File(self.data_folder + '/GM12878_DNase_ENCFF020WZB_2kb_4DNFI1UEG1HD_promoter_enhancer_encoding.h5', 'r')

        self.expr_df = pd.read_csv(self.data_folder + '/GM12878_K562_18377_gene_expr_fromXpresso.csv', index_col='ENSID')
        self.present_genes = self.expr_df.index
        self.valid_events = self.get_valid_genes()
        self.idx_map = None
        if single_event_train and (self.expr_type == 'multi' or self.expr_type == 'transcript'):
            self.idx_map, self.event_keys = self.map_idx_single_genes()

        if self.rna_seq_source == 'epiatlas':
            self.epiatlas_expr_df = pd.read_csv(self.data_folder + '/GM12878_K562_18360_gene_expr_epiatlas.csv', index_col='ENSID')
            self.present_genes = self.epiatlas_expr_df.index

        
    def __len__(self): # changed to filter for events 
        return len(self.event_keys)

    def __getitem__(self, idx):
        
        if self.expr_type == 'multi' or self.expr_type == 'transcript':
            if self.idx_map is not None:
                gene_id = self.event_keys[idx].split(";")[0]
                idx = self.idx_map[gene_id]
            event = self.valid_events[idx]
            gene_id = event.split(";")[0]
            # find idx where gene_id is in the data_h5
            idx = np.where(self.data_h5['ensid'][:] == gene_id.encode())[0][0]

        sample_ensid = self.data_h5['ensid'][idx].decode()
        seq_code = self.data_h5['pe_code'][idx]
        enhancer_distance = self.data_h5['distance'][idx,1:]
        enhancer_intensity = self.data_h5['activity'][idx,1:]
        enhancer_contact = self.data_h5['hic'][idx,1:]

        # added exon & intron sequences
        segment_tensor = torch.Tensor([])
        if self.expr_type == 'multi' or self.expr_type == 'transcript':
            upstream, downstream, exon = None, None, None
            sequences = self.gene_sequences[event]
            for key in sequences.keys():
                if 'upstream' in key:
                    upstream = sequences[key][()].decode()
                    upstream = upstream[-1024:]
                elif 'downstream' in key:
                    downstream = sequences[key][()].decode()
                    downstream = downstream[:1024]
                elif 'exon' in key:
                    exon = sequences[key][()].decode()
                    if len(exon) > 1024:
                        exon_start = exon[:512]
                        exon_end = exon[-512:]
                        exon = exon_start + exon_end
        
            # one hot encode the sequences
            vocab = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
            segment_tensor = torch.stack([self.one_hot_encode(upstream, vocab), self.one_hot_encode(exon, vocab), 
                                        self.one_hot_encode(downstream, vocab)])
        
        if self.signal_type == 'H3K27ac':
            promoter_activity = self.promoter_df.loc[sample_ensid]['PromoterActivity']
        elif self.signal_type == 'DNase':
            promoter_activity = self.promoter_df.loc[sample_ensid]['normalized_dhs']
            # enhancer_intensity = dhs_intensity
        promoter_code = seq_code[:1]
        enhancers_code = seq_code[1:]

        rnaFeat = list(self.expr_df.loc[sample_ensid][['UTR5LEN_log10zscore','CDSLEN_log10zscore',
                                                       'INTRONLEN_log10zscore','UTR3LEN_log10zscore','UTR5GC','CDSGC',
                                                       'UTR3GC', 'ORFEXONDENSITY']].values.astype(float))
        pe_activity = np.concatenate([[0], enhancer_intensity]).flatten()

        if self.usePromoterSignal and self.n_extraFeat > 1:
            rnaFeat = np.array(rnaFeat + [promoter_activity])
        else:
            rnaFeat = np.array(rnaFeat + [0])

        if self.distance_threshold is not None:
            enhancer_distance = enhancer_distance.flatten()
            enhancers_zero = np.zeros_like(enhancers_code)
            enhancers_zero[abs(enhancer_distance) < self.distance_threshold] = enhancers_code[abs(enhancer_distance) < self.distance_threshold]
            enhancers_code = enhancers_zero

            enhancer_distance_zero = np.zeros_like(enhancer_distance)
            enhancer_distance_zero[abs(enhancer_distance) < self.distance_threshold] = enhancer_distance[abs(enhancer_distance) < self.distance_threshold]
            enhancer_distance = enhancer_distance_zero

        if self.hic_threshold is not None:
            enhancer_contact = enhancer_contact.flatten()
            enhancers_zero = np.zeros_like(enhancers_code)
            enhancers_zero[enhancer_contact > self.hic_threshold] = enhancers_code[enhancer_contact > self.hic_threshold]
            enhancers_code = enhancers_zero

            enhancer_contact_zero = np.zeros_like(enhancer_contact)
            enhancer_contact_zero[enhancers_code[enhancer_contact > self.hic_threshold ]] = enhancer_contact[enhancers_code[enhancer_contact > self.hic_threshold]]
            enhancer_contact = enhancer_contact_zero

        pe_hic = np.concatenate([[0], enhancer_contact]).flatten()
        pe_hic = np.log10(1+pe_hic)
        pe_distance = np.concatenate([[0], enhancer_distance/1000]).flatten()
        # print(pe_distance)
        if self.n_extraFeat == 1:
            pe_feat = np.concatenate([pe_distance[:,np.newaxis]],axis=-1)
        elif self.n_extraFeat == 2:
            pe_feat = np.concatenate([pe_distance[:,np.newaxis], pe_activity[:,np.newaxis]],axis=-1)
        elif self.n_extraFeat == 3:
            pe_feat = np.concatenate([pe_distance[:,np.newaxis], pe_hic[:,np.newaxis], pe_activity[:,np.newaxis], ],axis=-1)
        else:
            pe_feat = np.concatenate([pe_distance[:,np.newaxis]],axis=-1)

        promoter_code_tensor = torch.from_numpy(promoter_code).float()
        pe_feat_tensor = torch.from_numpy(pe_feat[:self.n_enhancers+1])
        if self.n_extraFeat == 0: # Use promoter only
            enhancers_code = np.zeros_like(enhancers_code[:self.n_enhancers, :])
        enhancers_code_tensor = torch.from_numpy(enhancers_code[:self.n_enhancers, :]).float()
        pe_code_tensor = torch.concat([promoter_code_tensor, enhancers_code_tensor])
        rnaFeat_tensor = torch.from_numpy(rnaFeat).float()
        # print(pe_distance_tensor)
        psi_tensor = 0
        normal = "_normal" if self.use_normalized_psi else ""

        if self.expr_type == 'CAGE':
            cage_expr = np.log10(self.expr_df.loc[sample_ensid][self.cell_type + '_CAGE_128*3_sum']+1)
            expr_tensor = torch.from_numpy(np.array([cage_expr])).float()
        elif self.expr_type == 'RNA' and self.rna_seq_source == 'xpresso':
            rna_expr = self.expr_df.loc[sample_ensid]['Actual_' + self.cell_type]
            expr_tensor = torch.from_numpy(np.array([rna_expr])).float()
        elif self.expr_type == 'RNA' and self.rna_seq_source == 'epiatlas':
            rna_expr = self.epiatlas_expr_df.loc[sample_ensid][self.cell_type]
            expr_tensor = torch.from_numpy(np.array([rna_expr])).float()
            
        elif self.expr_type == 'multi':
            event_expr = self.psi_response.loc[event]
            expr_tensor = torch.from_numpy(np.array(event_expr[f'{self.cell_type}{self.tpm_level}{normal}'])).float()
            psi_tensor = torch.from_numpy(np.array(event_expr[f'{self.cell_type}_SE_psi{normal}'])).float()
        
        elif self.expr_type == 'transcript':
            event_expr = self.psi_response.loc[event]
            expr_tensor = torch.from_numpy(np.array(event_expr[f'{self.cell_type}{self.tpm_level}{normal}'])).float()
            psi_tensor = torch.from_numpy(np.array(event_expr[f'{self.cell_type}_SE_psi{normal}'])).float()
        
        else:
            assert False, 'Label does not exist!'
        
        return pe_code_tensor, segment_tensor, rnaFeat_tensor, pe_feat_tensor, expr_tensor, psi_tensor, sample_ensid

    def get_valid_genes(self):
        if not self.expr_type == 'multi' and not self.expr_type == 'transcript':
            return [x.decode() for x in self.data_h5['ensid'][:]]
        gene_ids = [x.split(";")[0] for x in self.event_keys]
        print(f"Found {len(gene_ids)} unique genes in the dataset, {len(set(gene_ids))} after removing duplicates.")
        ensid_list = set([x.decode() for x in self.data_h5['ensid'][:]])
        found_gene_ids = set([x for x in gene_ids if x in ensid_list])
        return [x for x in self.event_keys if x.split(";")[0] in found_gene_ids]

    def one_hot_encode(self, seq, vocab, length=1024):
        indices = [vocab[item] for item in seq]
        tensor = torch.tensor(indices)
        one_hot = torch.nn.functional.one_hot(tensor, num_classes=len(vocab)).float()
        # add padding
        if len(seq) < length:
            one_hot = torch.cat([one_hot, torch.zeros(length - len(seq), 4)], dim=0)
        elif len(seq) > length:
            one_hot = one_hot[:length]
        return one_hot
    
    def map_idx_single_genes(self):
        """
        Maps the indices of single genes in the dataset.
        Returns a dictionary with gene IDs as keys and their corresponding indices as values.
        """
        idx_map = {}
        seen_genes = set()
        for idx, event in enumerate(self.valid_events):
            gene_id = event.split(";")[0]
            if gene_id not in seen_genes:
                seen_genes.add(gene_id)
                idx_map[gene_id] = idx
        # only the idx we now have in the map should be valid events
        valid_events = np.array(self.valid_events)[list(idx_map.values())]
        return idx_map, list(valid_events)

    def z_score_normalize(self, train_idx):
        """
        Z-score normalize the psi responses for the training, validation, and test sets.
        """
        # get all psi responses for train_idx
        train_data = self.psi_response.iloc[train_idx]
        psi_responses = np.array(train_data[f'{self.cell_type}_SE_psi'])
        expr_responses = np.array(train_data[f'{self.cell_type}{self.tpm_level}'])
        mean_psi = np.mean(psi_responses)
        std_psi = np.std(psi_responses)
        mean_expr = np.mean(expr_responses)
        std_expr = np.std(expr_responses)
        self.use_normalized_psi = True

        self.psi_response[f'{self.cell_type}_SE_psi_normal'] = (
            self.psi_response[f'{self.cell_type}_SE_psi'] - mean_psi) / std_psi
        
        self.psi_response[f'{self.cell_type}{self.tpm_level}_normal'] = (
            self.psi_response[f'{self.cell_type}{self.tpm_level}'] - mean_expr) / std_expr

        return {'mean_psi': mean_psi, 'std_psi': std_psi, 
                'mean_expr': mean_expr, 'std_expr': std_expr}