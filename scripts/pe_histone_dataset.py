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
                 epigen_bp=False, enhancer_bp=False, include_enhancers=True, use_junctions=False, junction_length=250,
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
        self.use_epigen_bp = epigen_bp
        self.use_enhancer_bp = enhancer_bp
        self.include_enhancers = include_enhancers
        self.include_exon_enhancers = include_enhancers
        self.one_tpm_ar = one_tpm_ar
        self.use_junctions = use_junctions
        self.junction_length = junction_length
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
        self.gene_sequences = h5py.File(self.data_folder + '/event_encoding.h5', 'r')
        self.event_histone_seqs = h5py.File(self.data_folder + '/event_bp_histone_marks.h5', 'r')
        self.event_keys = [x.decode() for x in self.gene_sequences['event_id'][:]] # type:ignore
        if self.remove_ar:
            self.event_keys = [x for x in self.event_keys if ';AR:' not in x]
        self.psi_response = pd.read_csv(self.data_folder + '/psi_response.csv', index_col=0)
        # copy index to a new column for easier access
        self.psi_response['event_id'] = self.psi_response.index
        self.psi_response['event_type'] = self.psi_response['event_id'].apply(lambda x: 'AR' if ';AR:' in x else 'SE')

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
        self.data_dict['K562_bp'] = h5py.File(self.data_folder + '/K562_histone_enhancer_marks.h5', 'r')
        self.data_dict['GM12878_bp'] = h5py.File(self.data_folder + '/GM12878_histone_enhancer_marks.h5', 'r')

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
        if self.cell_type == 'K562':
            data_h5 = self.data_dict['K562']
            promoter_df = self.promoter_dict['K562']
            curr_cell_type = 'K562'
        elif self.cell_type == 'GM12878':
            data_h5 = self.data_dict['GM12878']
            promoter_df = self.promoter_dict['GM12878']
            curr_cell_type = 'GM12878'
        elif self.cell_type == 'both':
            if idx < len(self.event_keys):
                data_h5 = self.data_dict['K562']
                promoter_df = self.promoter_dict['K562']
                curr_cell_type = 'K562'
            else:
                idx -= len(self.event_keys)
                data_h5 = self.data_dict['GM12878']
                promoter_df = self.promoter_dict['GM12878']
                curr_cell_type = 'GM12878'
        else:
            raise ValueError(f"Cell type {self.cell_type} not supported. Choose 'K562' or 'GM12878'.")
        
        cell_tensor = torch.tensor([0, 1, 0]) if curr_cell_type == 'K562' else torch.tensor([1, 0, 0])
        histone_marks = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3', 'methylation']

        if self.include_exons:
            if self.idx_map is not None:
                gene_id = self.event_keys[idx].split(";")[0]
                idx = self.idx_map[gene_id]
            event = self.valid_events[idx]
            gene_id = event.split(";")[0]
            # find idx where gene_id is in the data_h5
            idx = np.where(data_h5['gene_id'][:] == gene_id.encode())[0][0]

        if self.filter_for_event_genes and self.expr_type != 'multi' and self.expr_type != 'splice':
            idx = self.idx_map[idx]

        sample_ensid = data_h5['gene_id'][idx].decode()
        seq_code = data_h5['pe_seqs'][idx]
        enhancer_distance = data_h5['distance'][idx,1:] # from the old h5, thus 1:
        enhancer_data = {key: data_h5[key][idx,:] for key in histone_marks if key in data_h5}
        enhancer_intensity = enhancer_data.pop('H3K27ac', None)
        enhancer_histones = np.stack([x for x in enhancer_data.values()], axis=-1)
        enhancer_epigen_tensor = torch.Tensor([])
        if self.use_enhancer_bp:
            enhancer_epigen = [self.data_dict[f"{curr_cell_type}_bp"][x][idx] for x in histone_marks]
            enhancer_epigen = np.concatenate(enhancer_epigen, axis=-1)
            enhancer_epigen_tensor = torch.from_numpy(enhancer_epigen)

        # added exon & intron sequences
        segment_tensor = torch.Tensor([])
        event_hist_tensor = torch.Tensor([])
        if self.include_exons:
            # get the idx at which event == self.gene_sequences['event_id'][idx].decode()
            event_idx = np.where(self.gene_sequences['event_id'][:] == event.encode())[0][0]
            assert event == self.gene_sequences['event_id'][event_idx].decode(), "Event ID mismatch!"
            segment_tensor = torch.from_numpy(self.gene_sequences['event_seq'][event_idx])
            event_histones = [self.gene_sequences[f"{x}_{self.cell_type}"][event_idx,:] for x in histone_marks 
                              if f"{x}_{self.cell_type}" in self.gene_sequences]
            event_hist_tensor = torch.from_numpy(np.stack(event_histones, axis=-1))
            # add a zero tensor to emulate distance
            event_hist_tensor = torch.cat([torch.zeros((event_hist_tensor.shape[0], 1)), event_hist_tensor], dim=-1)

        bp_epigen_marks = None
        if self.use_epigen_bp:
            if not self.include_exons:
                raise ValueError("Epigenetic data on base pair level can only be used when including exons.")
            bp_epigen_marks = [self.event_histone_seqs[self.cell_type][x][event_idx] for x in histone_marks]
            bp_epigen_marks = np.concatenate(bp_epigen_marks, axis=-1)
            bp_epigen_marks = torch.from_numpy(bp_epigen_marks)

        if self.use_junctions and self.include_exons:
            segment_tensor, bp_epigen_marks = self.get_junctions(segment_tensor, bp_epigen_marks)
            # combine 1+2 and 2+3 and get averages to get from shape [3, 7] to [2, 7]
            junction1 = (event_hist_tensor[0] + event_hist_tensor [1]) / 2
            junction2 = (event_hist_tensor[1] + event_hist_tensor [2]) / 2
            event_hist_tensor = torch.stack([junction1, junction2], dim=0)

        promoter_code = seq_code[:1]
        enhancers_code = seq_code[1:]

        pe_activity = np.concatenate([[0], enhancer_intensity]).flatten()
        enhancer_histones = np.vstack([np.zeros((1, enhancer_histones.shape[1])), enhancer_histones])

        promoter_histones = self.promoter_dict[curr_cell_type].loc[sample_ensid][['H3K27me3.mean','H3K36me3.mean',
                                                                                 'H3K4me1.mean','H3K4me3.mean',
                                                                                 'H3K9me3.mean', 'H3K27ac.mean']]

        promoter_histones = [0, 0, 0] + list(promoter_histones.values.astype(float)) # add three zeros to keep dims

        if not self.usePromoterSignal or self.n_extraFeat <= 1:
            promoter_histones = promoter_histones[:-1] # remove last feature which is promoter activity
        promoter_histones = np.array(promoter_histones, dtype=np.float32)

        if self.distance_threshold is not None:
            enhancer_distance = enhancer_distance.flatten()
            enhancers_zero = np.zeros_like(enhancers_code)
            enhancers_zero[abs(enhancer_distance) < self.distance_threshold] = enhancers_code[abs(enhancer_distance) < self.distance_threshold]
            enhancers_code = enhancers_zero

            enhancer_distance_zero = np.zeros_like(enhancer_distance)
            enhancer_distance_zero[abs(enhancer_distance) < self.distance_threshold] = enhancer_distance[abs(enhancer_distance) < self.distance_threshold]
            enhancer_distance = enhancer_distance_zero

        pe_distance = np.concatenate([[0], enhancer_distance/1000]).flatten()
        # print(pe_distance)
        if self.n_extraFeat == 1:
            pe_feat = np.concatenate([pe_distance[:,np.newaxis]],axis=-1)
        elif self.n_extraFeat == 2:
            pe_feat = np.concatenate([pe_distance[:,np.newaxis], pe_activity[:,np.newaxis]],axis=-1)
        else:
            pe_feat = np.concatenate([pe_distance[:,np.newaxis]],axis=-1)

        if self.include_enhancers:
            pe_feat = np.concatenate([pe_feat, enhancer_histones], axis=-1)
        
        if self.include_exon_enhancers:
            pe_feat = np.concatenate([pe_feat, event_hist_tensor], axis=0)

        promoter_code_tensor = torch.from_numpy(promoter_code).float()
        # zero out promoter code if set_promoter_zero is True
        if self.zero_out_promoter:
            promoter_code_tensor = torch.zeros_like(promoter_code_tensor)

        extra_dims = 4 if self.include_exon_enhancers else 1
        pe_feat_tensor = torch.from_numpy(pe_feat[:self.n_enhancers + extra_dims])
        if self.n_extraFeat == 0: # Use promoter only
            enhancers_code = np.zeros_like(enhancers_code[:self.n_enhancers, :])
        enhancers_code_tensor = torch.from_numpy(enhancers_code[:self.n_enhancers, :]).float()
        pe_code_tensor = torch.concat([promoter_code_tensor, enhancers_code_tensor])
        histone_features_tensor = torch.from_numpy(promoter_histones).float()

        # zero out promoter/enhancer data for ablation study
        if self.zero_out_pe_data:
            pe_code_tensor = torch.zeros_like(pe_code_tensor)
        if self.zero_out_histone_data:
            # replace everything but the promoter activity with zeros
            if self.usePromoterSignal and self.n_extraFeat > 1:
                histone_features_tensor[:-1] = 0
            else:
                histone_features_tensor = torch.zeros_like(histone_features_tensor)
        if self.zero_out_feat_data:
            pe_feat_tensor = torch.zeros_like(pe_feat_tensor)
        if self.zero_out_exons:
            segment_tensor = torch.zeros_like(segment_tensor)

        # print(pe_distance_tensor)
        ## Retrieve the true values for the expression and psi tensors
        psi_tensor = 0
        normal = "_normal" if self.use_normalized_psi else ""

        if self.expr_type == 'CAGE':
            cage_expr = np.log10(self.expr_df.loc[sample_ensid][curr_cell_type + '_CAGE_128*3_sum']+1)
            expr_tensor = torch.from_numpy(np.array([cage_expr])).float()
        elif self.expr_type == 'RNA' and self.rna_seq_source == 'xpresso':
            rna_expr = self.expr_df.loc[sample_ensid]['Actual_' + curr_cell_type]
            expr_tensor = torch.from_numpy(np.array([rna_expr])).float()
        elif self.expr_type == 'RNA' and self.rna_seq_source == 'epiatlas':
            rna_expr = self.epiatlas_expr_df.loc[sample_ensid][curr_cell_type]
            expr_tensor = torch.from_numpy(np.array([rna_expr])).float()
            
        elif self.expr_type == 'multi' or self.expr_type == 'splice':
            event_expr = self.psi_response.loc[event]
            expr_tensor = torch.from_numpy(np.array(event_expr[f'{curr_cell_type}{self.tpm_level}{normal}'])).float()
            psi_tensor = torch.from_numpy(np.array(event_expr[f'{curr_cell_type}_SE_psi{normal}'])).float()
        
        else:
            assert False, 'Label does not exist!'

        return {'pe_seq': pe_code_tensor, 'segment_seq': segment_tensor, 'histone_features': histone_features_tensor,
                'ev_epigen_seq': bp_epigen_marks, 'en_epigen_seq': enhancer_epigen_tensor, 'pe_feat': pe_feat_tensor, 
                'cell': cell_tensor, 'expr': expr_tensor, 'psi': psi_tensor, 'sample': sample_ensid}
    

    def get_junctions(self, segment_tensor, bp_epigen_marks=None):
        ex_mask = (segment_tensor[1] == 0).all(dim=1)
        ex_pad_start = ex_mask.nonzero(as_tuple=True)[0][0] if ex_mask.any() else segment_tensor[1].shape[0]
        in_mask = (segment_tensor[0] == 0).all(dim=1)
        in_pad_start = in_mask.nonzero(as_tuple=True)[0][0] if in_mask.any() else segment_tensor[0].shape[0]
        half_len = self.junction_length // 2
        junction1 = torch.cat([segment_tensor[0, in_pad_start-half_len:in_pad_start, :], segment_tensor[1, :half_len, :]], dim=0)
        junction2 = torch.cat([segment_tensor[1, ex_pad_start-half_len:ex_pad_start, :], segment_tensor[2, :half_len, :]], dim=0)
        epigen_junction1, epigen_junction2 = torch.zeros_like(junction1), torch.zeros_like(junction2)
        if bp_epigen_marks is not None:
            epigen_junction1 = torch.cat([bp_epigen_marks[0, in_pad_start-half_len:in_pad_start, :], 
                                        bp_epigen_marks[1, :half_len, :]], dim=0)
            epigen_junction2 = torch.cat([bp_epigen_marks[1, ex_pad_start-half_len:ex_pad_start, :],
                                         bp_epigen_marks[2, :half_len, :]], dim=0)
        # if junctions less than [240, 10], pad first dim with 0
        if junction1.shape[0] < self.junction_length:
            pad = self.junction_length - junction1.shape[0]
            junction1 = F.pad(junction1, (0, 0, pad, 0))
            epigen_junction1 = F.pad(epigen_junction1, (0, 0, pad, 0))
        if junction2.shape[0] < self.junction_length:
            pad = self.junction_length - junction2.shape[0]
            junction2 = F.pad(junction2, (0, 0, pad, 0))
            epigen_junction2 = F.pad(epigen_junction2, (0, 0, pad, 0))
        # # 2, self.junction_length, 10
        return torch.stack([junction1, junction2], dim=0), torch.stack([epigen_junction1, epigen_junction2], dim=0)


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
