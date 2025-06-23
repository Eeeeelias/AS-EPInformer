from EPInformer.models_multi import EPInformer_v2, enhancer_predictor_256bp
from scripts.utils import prepare_input
import scripts.utils_forTraining as train
import scripts.promoter_enhancer_dataset as ped
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from tqdm import tqdm
import torch
from torch.utils.data import Subset, Dataset
from sklearn.metrics import r2_score
import torch.utils.data as data_utils
import scripts.promoter_enhancer_dataset as pe_dataset
from collections import defaultdict
import random

split_df = pd.read_csv('./data/leave_chrom_out_crossvalidation_split_18377genes.csv', index_col=0)

def split_multitask_ids(ids: list[str], train_frac: float = 0.7, val_frac: float = 0.15, test_frac: float = 0.15, 
                        seed: int = 42) -> tuple[list[str], list[str], list[str]]:
    """
    Splits event IDs into train/val/test sets, grouping by gene of each ID.
    
    Parameters:
    - ids: List of ID strings, each with parts separated by ';'.
    - train_frac, val_frac, test_frac: Fractions for each split (should sum to 1).
    - seed: Random seed for reproducibility.
    
    Returns:
    - A tuple of three lists: (train_ids, val_ids, test_ids)
    """
    print(f"Splitting gene-aware")
    assert abs(train_frac + val_frac + test_frac - 1.0) < 1e-6, "Fractions must sum to 1."

    # Group full IDs by their shared key (0-th element of the split)
    key_to_indices = defaultdict(list)
    for idx, id_ in enumerate(ids):
        key = id_.split(';')[0]
        key_to_indices[key].append(idx)

    # Shuffle keys deterministically
    random.seed(seed)
    all_keys = list(key_to_indices.keys())
    random.shuffle(all_keys)

    # Compute split cutoffs
    n = len(all_keys)
    train_cutoff = int(train_frac * n)
    val_cutoff = int((train_frac + val_frac) * n)

    train_keys = all_keys[:train_cutoff]
    val_keys = all_keys[train_cutoff:val_cutoff]
    test_keys = all_keys[val_cutoff:]

    # Gather indices
    train_indices = [idx for k in train_keys for idx in key_to_indices[k]]
    val_indices = [idx for k in val_keys for idx in key_to_indices[k]]
    test_indices = [idx for k in test_keys for idx in key_to_indices[k]]

    return train_indices, val_indices, test_indices

def filter_id_lists(existing_ids, train_ids, valid_ids, test_ids):
    """
    Filter the existing IDs to only include those exist in the data
    """
    filtered_train_ids = [i for i in train_ids if i in existing_ids]
    filtered_valid_ids = [i for i in valid_ids if i in existing_ids]
    filtered_test_ids = [i for i in test_ids if i in existing_ids]
    
    return filtered_train_ids, filtered_valid_ids, filtered_test_ids

cell = 'K562'
distance_threshold = 100_000
n_enhancers = 60
device = 'xpu'
# num_feature == 1: distance; num_feature == 2: distance + enhancer activity; num_feature == 3: distance + enhancer activity + hic contacts
n_extraFeat = 3
batch_size = 16
expr_type = 'RNA'
prediction_res = []
fi = 1
fold_i = 'fold_1'

activations = {}

def save_activation(name):
    def hook(module, input, output):
        if isinstance(output, tuple):
            activations[name + "_output"] = output[0].detach().cpu()
            activations[name + "_attn"] = output[1].detach().cpu()
        else:
            activations[name] = output.detach().cpu()
    return hook

train_ensid = split_df[split_df[fold_i] == 'train'].index
valid_ensid = split_df[split_df[fold_i] == 'valid'].index
test_ensid = split_df[split_df[fold_i] == 'test'].index

all_ds = pe_dataset.promoter_enhancer_dataset(data_folder= './data/', expr_type=expr_type, cell_type=cell,
                                                  n_extraFeat=n_extraFeat, usePromoterSignal=True,
                                                  n_enhancers=n_enhancers, hic_threshold=None,
                                                  distance_threshold=distance_threshold, include_exons=True,
                                                  rna_seq_source='xpresso', tpm='gene',
                                                  single_event_train=False, event_genes=False)


train_idx, valid_idx, test_idx = split_multitask_ids(all_ds.event_keys, train_frac=0.8, val_frac=0.1, 
                                                             test_frac=0.1, seed=42+int(fi))

train_ds = Subset(all_ds, train_idx)
valid_ds = Subset(all_ds, valid_idx)
test_ds = Subset(all_ds, test_idx)
pretrained_convNet = enhancer_predictor_256bp()
pt_model_name = '{}_seq2activityLog2_leaveChrOut_combinedRS_2bins_bs64_H3K27ac_adamW_erisxdl_r0'.format(cell)
checkpoint = torch.load(f"./trained_models/pretrained_enhancer_encoder/{fold_i}_best_{pt_model_name}_checkpoint.pt", 
                        weights_only=False, map_location=device)
model = EPInformer_v2(n_encoder=6, pre_trained_encoder=pretrained_convNet.encoder,
                              n_enhancer=n_enhancers, out_dim=64, n_extraFeat=n_extraFeat, device=device, exon_data=True)
model = model.to(device)
# model_name= 'tunedEPInformerV2.preTrainedConv.4base.64dim.3Trans.4head.TrueBN.TrueLN.TrueFeat.3extraFeat.60enh.K562.rmEnhNone.bs16.seq_feat_dist.DNaseH.distanceDist100k.hic0.len2k.distance.{}'.format(expr_type)
# checkpoint = torch.load("./trained_models/EPInformer_PE_Activity_HiC/K562/fold_{}_EPInformer_PE_Activity_HiC_{}_K562_checkpoint.pt".format(fi, expr_type), weights_only=False)
checkpoint = torch.load("trained_models/2025-06-21-11/fold_1_best_EPInformer-PE-Activity-HiC.preTrainedConv.4base.64dim.6Trans.8head.TrueBN.TrueLN.TrueFeat.3extraFeat.60enh.Trueexon.K562.multi_checkpoint.pt", weights_only=False, map_location=device)

model.load_state_dict(checkpoint['model_state_dict'])

# register all attn_encoders
for i in range(model.n_encoder):
    model.attn_encoder[i].register_forward_hook(save_activation(f"encoder_layer_{i}"))

testloader = data_utils.DataLoader(test_ds, batch_size=batch_size, pin_memory=True, num_workers=0)
model.eval()
with torch.no_grad():
    preds = []
    actual = []
    preds_psi = []
    actual_psi = []
    uncorr_pred_expr = []
    uncorr_pred_splice = []
    uncorr_actual_expr = []
    uncorr_actual_psi = []

    ensid_list = []
    for data in testloader:
        input_PE, input_seg, input_feat, input_dist, y_expr, y_psi, eid = data
        input_PE = input_PE.float().to(device)
        input_seg = input_seg.float().to(device)
        input_feat = input_feat.float().to(device)
        # input_dist = input_dist.long().to(device)
        input_dist = input_dist.float().to(device)
        y_expr = y_expr.float().to(device)
        y_psi = y_psi.float().to(device)
        print(f"input_PE shape: {input_PE.shape}, input_seg shape: {input_seg.shape}, input_feat shape: {input_feat.shape}, input_dist shape: {input_dist.shape}")
        _ = model(input_PE, input_seg, input_feat, input_dist)
        break


for name, act in activations.items():
    print(f"{name}: mean={act.mean().item():.4f}, std={act.std().item():.4f}, min={act.min().item():.4f}, max={act.max().item():.4f}")
