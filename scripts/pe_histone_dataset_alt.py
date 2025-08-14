## ALTERNATIVE VERSION USED ONLY FOR EXTENDED SE EVENTS NOT CONNECTED WITH EPINFORMER 

import torch
import os
import numpy as np
import pandas as pd
# torch
from torch.nn.functional import one_hot
from torch.utils.data import Dataset
import torch.nn.functional as F

# from model.EPInformer import EPInformer_v2, enhancer_predictor_256bp
import h5py


class PEHistoneDataset(Dataset):
    def __init__(self, data_folder = 'data/', expr_type='CAGE', usePromoterSignal=True, first_signal='distance', signal_type='H3K27ac', 
                 cell_type='K562', distance_threshold=None, hic_threshold=None, n_enhancers=50, n_extraFeat=1,
                 rna_seq_source='xpresso', tpm='gene', single_event_train=False, event_genes=False, include_exons=False,
                 include_enhancers=True,
                 set_exon_zero=False, set_pe_zero=False, set_histones_zero=False, set_extra_feat_zero=False, 
                 set_promoter_zero=False, remove_ar=False, one_tpm_ar=False, **kwargs):
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
        self.filter_for_event_genes = event_genes
        self.include_exons = include_exons
        self.include_enhancers = include_enhancers
        self.include_exon_enhancers = include_enhancers
        self.one_tpm_ar = one_tpm_ar
        self.promoter_dict = {}
        self.data_dict = {}
        # ablation test params
        self.zero_out_exons = set_exon_zero
        self.zero_out_pe_data = set_pe_zero
        self.zero_out_histone_data = set_histones_zero
        self.zero_out_feat_data = set_extra_feat_zero
        self.zero_out_promoter = set_promoter_zero
        self.remove_ar = remove_ar
        
        if not self.include_exons and (self.expr_type == 'multi' or self.expr_type == 'splice'):
            print("INFO: Exons will be included automatically for multi or transcript expression types.")
            self.include_exons = True

        self.use_normalized_psi = False
        self.tpm_level = "_summed_tpm" if tpm == 'transcript' else "_gene_level_tpm"
        self.gene_sequences = h5py.File(self.data_folder + f'/extended_events_{self.cell_type}.h5', 'r')
        self.event_keys = [x.decode() for x in self.gene_sequences['event_id'][:]] # type: ignore
        if self.remove_ar:
            self.event_keys = [x for x in self.event_keys if ';AR:' not in x]
        self.psi_response = pd.read_csv(self.data_folder + f'/extended_psi_response_{self.cell_type}.csv', index_col=0)
        # copy index to a new column for easier access
        self.psi_response['event_id'] = self.psi_response.index
        self.psi_response['event_type'] = "SE"

        if self.one_tpm_ar: # remove ar events with less than 1 tpm in the cell line
            self.event_keys = [x for x in self.event_keys if self.psi_response.loc[x, f'{self.cell_type}{self.tpm_level}'] >= np.log10(1+0.1)]

        self.all_event_genes = set(self.psi_response['gene_id'].unique())
        # K562 promoter data
        promoter_df = pd.read_csv(self.data_folder + 'IHEC-ChIP-Seq-Histone-Signals/Promoter_Combined_K562_Histone_Signals.csv', index_col='gene_id')
        self.promoter_dict['K562'] = promoter_df
        # GM12878 promoter data
        promoter_df = pd.read_csv(self.data_folder + 'IHEC-ChIP-Seq-Histone-Signals/Promoter_Combined_GM12878_Histone_Signals.csv', index_col='gene_id')
        self.promoter_dict['GM12878'] = promoter_df

        self.data_dict['K562'] = h5py.File(self.data_folder + '/K562_histone_appended_pe_encoding.h5', 'r')
        self.data_dict['GM12878'] = h5py.File(self.data_folder + '/GM12878_histone_appended_pe_encoding.h5', 'r')

        self.expr_df = pd.read_csv(self.data_folder + '/GM12878_K562_18377_gene_expr_fromXpresso.csv', index_col='ENSID')
        self.present_genes = self.expr_df.index
        self.valid_events = self.get_valid_genes()
        self.idx_map = None # map various indices to gene ids
        if single_event_train and (self.expr_type == 'multi' or self.expr_type == 'splice'):
            self.idx_map, self.event_keys = self.map_idx_single_genes()

        if self.rna_seq_source == 'epiatlas':
            self.epiatlas_expr_df = pd.read_csv(self.data_folder + '/GM12878_K562_18360_gene_expr_epiatlas.csv', index_col='ENSID')
            self.present_genes = self.epiatlas_expr_df.index

        if self.filter_for_event_genes:
            self.idx_map = {}
            # use idx map to filter for event genes
            c = 0
            for i, gene in enumerate(self.data_dict[cell_type]['gene_id'][:]):
                gene_id = gene.decode()
                if gene_id in self.all_event_genes:
                    self.idx_map[c] = i
                    c += 1

    def __len__(self): # changed to filter for events 
        if self.include_exons and self.cell_type == 'both':
            return len(self.event_keys) * 2
        if self.include_exons:
            return len(self.event_keys)
        if self.filter_for_event_genes:
            return len(self.all_event_genes)
        
        return len(self.data_dict[self.cell_type]['gene_id'])

    def __getitem__(self, idx):       
        histone_marks = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']

        # added exon & intron sequences
        segment_tensor = torch.Tensor([])

        sequences = self.gene_sequences['event_data'][idx] # (3, 1024, 10)
        event_name = self.gene_sequences['event_id'][idx].decode() # get event name

        ### BINARY TESTING
        segment_tensor = torch.from_numpy(np.array(sequences)).float() # convert to tensor
        # for the second sequence, find at what idx the padding starts
        ex_mask = (segment_tensor[1] == 0).all(dim=1)
        ex_pad_start = ex_mask.nonzero(as_tuple=True)[0][0] if ex_mask.any() else segment_tensor[1].shape[0]
        in_mask = (segment_tensor[0] == 0).all(dim=1)
        in_pad_start = in_mask.nonzero(as_tuple=True)[0][0] if in_mask.any() else segment_tensor[0].shape[0]
        junction1 = torch.cat([segment_tensor[0, in_pad_start-120:in_pad_start, :], segment_tensor[1, :0+120, :]], dim=0)
        junction2 = torch.cat([segment_tensor[1, ex_pad_start-120:ex_pad_start, :], segment_tensor[2, :0+120, :]], dim=0)
        # if junctions less than [240, 10], pad first dim with 0
        if junction1.shape[0] < 240:
            pad = 240 - junction1.shape[0]
            junction1 = F.pad(junction1, (0, 0, pad, 0))
        if junction2.shape[0] < 240:
            pad = 240 - junction2.shape[0]
            junction2 = F.pad(junction2, (0, 0, pad, 0))
        segment_tensor = torch.stack([junction1, junction2], dim=0) # 2, 240, 10
        ### END

        if self.zero_out_exons:
            segment_tensor = torch.zeros_like(segment_tensor)

        # print(pe_distance_tensor)
        ## Retrieve the true values for the expression and psi tensors
        psi_tensor = 0

        event_expr = self.psi_response.loc[event_name]
        psi_tensor = torch.from_numpy(np.array(event_expr[f'{self.cell_type}_SE_psi'])).float()
        sample_ensid = event_expr['gene_id']
        
        return segment_tensor, psi_tensor, sample_ensid


    def get_valid_genes(self):
        if not self.include_exons:
            return [x.decode() for x in self.data_dict[self.cell_type]['gene_id'][:]]
        gene_ids = [x.split(";")[0] for x in self.event_keys]
        print(f"Found {len(gene_ids)} unique genes in the dataset, {len(set(gene_ids))} after removing duplicates.")
        if self.cell_type == 'both':
            ensid_list = set([x.decode() for x in self.data_dict['K562']['gene_id'][:]]) # just use K562 ensids since both cell lines have the same genes
        else:
            ensid_list = set([x.decode() for x in self.data_dict[self.cell_type]['gene_id'][:]])
        found_gene_ids = set([x for x in gene_ids if x in ensid_list])
        return [x for x in self.event_keys if x.split(";")[0] in found_gene_ids]

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
        # get boolean mask for fake events
        fake_events_train = train_data['event_type'] == 'AR'
        fake_events_all = self.psi_response['event_type'] == 'AR'
        psi_responses_real = psi_responses[~fake_events_train]
        expr_responses = np.array(train_data[f'{self.cell_type}{self.tpm_level}'])
        mean_psi = np.mean(psi_responses_real)
        std_psi = np.std(psi_responses_real)
        mean_expr = np.mean(expr_responses)
        std_expr = np.std(expr_responses)
        self.use_normalized_psi = True

        #self.psi_response[f'{self.cell_type}_SE_psi_normal'] = (
        #    self.psi_response[f'{self.cell_type}_SE_psi'] - mean_psi) / std_psi

        # first set placeholder value for fake events, then replace non-fake events with normalized values
        normalised_psi = (self.psi_response.loc[~fake_events_all, f'{self.cell_type}_SE_psi'] - mean_psi) / std_psi
        self.psi_response[f'{self.cell_type}_SE_psi_normal'] = np.ones_like(self.psi_response[f'{self.cell_type}_SE_psi'])
        self.psi_response.loc[~fake_events_all, f'{self.cell_type}_SE_psi_normal'] = normalised_psi
        
        
        self.psi_response[f'{self.cell_type}{self.tpm_level}_normal'] = (
            self.psi_response[f'{self.cell_type}{self.tpm_level}'] - mean_expr) / std_expr

        return {'mean_psi': mean_psi, 'std_psi': std_psi, 
                'mean_expr': mean_expr, 'std_expr': std_expr}
