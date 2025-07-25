import argparse
from datetime import datetime
import random
from collections import defaultdict
import os
import shutil
import sys
import yaml
from types import SimpleNamespace

from scipy import stats
from tqdm import tqdm
import pandas as pd
import numpy as np
import torch
from torch.utils.data import Subset, Dataset
import scripts.utils_forTraining as utils
import scripts.promoter_enhancer_dataset as pe_dataset
import scripts.pe_histone_dataset as pe_histone_dataset

def parse_yaml_config(config_path):
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def str_to_list(value):
    if isinstance(value, list):
        return value
    return value.split(',')

def dict_to_namespace(d):
    """Recursively converts a dict into a SimpleNamespace."""
    if isinstance(d, dict):
        return SimpleNamespace(**{k: dict_to_namespace(v) for k, v in d.items()})
    elif isinstance(d, list):
        return [dict_to_namespace(i) for i in d]
    else:
        return d

def filter_id_lists(existing_ids, train_ids, valid_ids, test_ids):
    """
    Filter the existing IDs to only include those exist in the data
    """
    filtered_train_ids = [i for i in train_ids if i in existing_ids]
    filtered_valid_ids = [i for i in valid_ids if i in existing_ids]
    filtered_test_ids = [i for i in test_ids if i in existing_ids]
    
    return filtered_train_ids, filtered_valid_ids, filtered_test_ids

def create_set_indices(all_ids, train_ratio=0.8, valid_ratio=0.1, events=False, seed=0, splits=None, fold_i='fold_1'):
    """
    Create train, valid, test indices based on the given ratios
    """
    np.random.seed(seed)
    if splits is None:
        all_ids = np.arange(len(all_ids))
        np.random.shuffle(all_ids)
    
    n_total = len(all_ids)
    n_train = int(n_total * train_ratio)
    n_valid = int(n_total * valid_ratio)
    
    if splits is not None:
        print("Using predefined splits")
        ensid_list = [eid.decode() for eid in all_ids.data_h5['ensid'][:]]
        ensid_df = pd.DataFrame(ensid_list, columns=['ensid'])
        ensid_df['idx'] = np.arange(len(ensid_list))
        ensid_df = ensid_df.set_index('ensid')

        train_ensid = splits[splits[fold_i] == 'train'].index
        valid_ensid = splits[splits[fold_i] == 'valid'].index
        test_ensid = splits[splits[fold_i] == 'test'].index

        train_ids = ensid_df.loc[train_ensid]['idx']
        valid_ids = ensid_df.loc[valid_ensid]['idx']
        test_ids = ensid_df.loc[test_ensid]['idx']
    else:
        train_ids = all_ids[:n_train]
        valid_ids = all_ids[n_train:n_train + n_valid]
        test_ids = all_ids[n_train + n_valid:]
    
    return train_ids, valid_ids, test_ids

def split_multitask_ids(ids: list[str], train_frac: float = 0.7, val_frac: float = 0.15, test_frac: float = 0.15, 
                        seed: int = 42, short_run=False, pre_splits=None, fold_i=None, both_cell_lines=False) -> tuple[list[str], list[str], list[str]]:
    """
    Splits event IDs into train/val/test sets, grouping by gene of each ID.
    
    Parameters:
    - ids: List of ID strings, each with parts separated by ';'.
    - train_frac, val_frac, test_frac: Fractions for each split (should sum to 1).
    - seed: Random seed for reproducibility.
    
    Returns:
    - A tuple of three lists: (train_ids, val_ids, test_ids)
    """
    print("Splitting gene-aware")
    assert abs(train_frac + val_frac + test_frac - 1.0) < 1e-6, "Fractions must sum to 1."

    # Group full IDs by their shared key (0-th element of the split)
    key_to_indices = defaultdict(list)
    for idx, id_ in enumerate(ids):
        key = id_.split(';')[0]
        key_to_indices[key].append(idx)

    # Shuffle keys deterministically
    if pre_splits is not None:
        print("Using predefined chromosome splits")
        # Use predefined splits from the DataFrame
        train_keys = pre_splits[pre_splits[fold_i] == 'train'].index.tolist()
        val_keys = pre_splits[pre_splits[fold_i] == 'valid'].index.tolist()
        test_keys = pre_splits[pre_splits[fold_i] == 'test'].index.tolist()

        train_indices = [idx for k in train_keys for idx in key_to_indices[k]]
        val_indices = [idx for k in val_keys for idx in key_to_indices[k]]
        test_indices = [idx for k in test_keys for idx in key_to_indices[k]]
        print(f"Using predefined splits: {len(train_indices)} train")

    else:
        # Shuffle keys randomly
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

    if both_cell_lines:
        # add the len(idx) to all indices since the len of the dataset is len(event_keys) * 2
        train_indices = train_indices + [idx + len(ids) for idx in train_indices]
        val_indices = val_indices + [idx + len(ids) for idx in val_indices]
        test_indices = test_indices + [idx + len(ids) for idx in test_indices]
        print(f"Both cell lines included, total indices: {len(train_indices)}")

    if short_run:
        # Limit the number of samples for a short run
        train_indices = train_indices[:2000]
        val_indices = val_indices[:100]
        test_indices = test_indices[:100]

    return train_indices, val_indices, test_indices

##### init ######
parser = argparse.ArgumentParser()
parser.add_argument('--config_path', type=str, required=True)
parser.add_argument('--print_config_help', action='store_true', help='Print help info for config keys')

args = parser.parse_args()
yml_config = parse_yaml_config(args.config_path)

if args.print_config_help:
    help_path = os.path.join(os.path.dirname(args.config_path), 'config_help.yaml')
    if os.path.exists(help_path):
        help_info = parse_yaml_config(help_path)
        for key, val in help_info.items():
            print(f"{key}: {val}")
    else:
        print("No config_help.yaml found.")
    sys.exit(0)

config = dict_to_namespace(yml_config)

#### import the right EPinformer model
if config.base.expr_assay == 'multi' or config.base.expr_assay == 'splice' or config.base.include_exons:
    from EPInformer.models_multi import EPInformer_v2, enhancer_predictor_256bp, WeightedLoss
else:
    from EPInformer.models import EPInformer_v2, enhancer_predictor_256bp, WeightedLoss


cell = config.base.cell
if config.hardware.cuda:
    device = 'cuda'
elif config.hardware.xpu:
    device = 'xpu'
else:
    device = 'cpu'

distance_threshold = config.base.distance_threshold
n_epoch = config.train.max_epochs
hic_threshold = config.base.hic_threshold
if hic_threshold == -1:
    hic_threshold = None

if config.base.model_type == 'EPInformer-PE': 
    n_extraFeat = 1
elif config.base.model_type == 'EPInformer-PE-Activity':
    n_extraFeat = 2
elif config.base.model_type == 'EPInformer-PE-Activity-HiC':
    n_extraFeat = 3
else:
    raise ValueError(f"Unsupported model type: {config.base.model_type}")

use_pretrained = config.hardware.use_pretrained_encoder
fold_list = str_to_list(config.base.fold)
n_encoder = config.train.n_interact_enc
batch_size = config.train.batch_size
expr_type = config.base.expr_assay
n_enhancers = 60
#################

today = datetime.now()   # Get date

datetime_str = today.strftime("%Y-%m-%d-%H-%M")
day_str = today.strftime("%Y-%m-%d")

if config.debug.short_run:
    saved_model_path = None
else:
    saved_model_path = f'./trained_models/{day_str}/{datetime_str}/'
    if not os.path.exists(saved_model_path):
        os.makedirs(saved_model_path, exist_ok=True)
    shutil.copy(args.config_path, os.path.join(saved_model_path, 'config.yaml'))

if 'all' in fold_list:
    fold_list = list(range(1, 13))

if config.debug.short_run:
    fold_list = fold_list[:1]  # For testing, only use the first fold

split_df = pd.read_csv('./data/leave_chrom_out_crossvalidation_split_18377genes.csv', index_col=0)


for fi in fold_list:
    print("-"*10, 'fold', fi, '-'*10)
    fold_i = 'fold_' + str(fi)
    
    # removed pre-set indices for train, valid, test since our data is now event-based, not gene-based
    ablation_tests = {"set_exon_zero": config.debug.set_exon_zero, "set_pe_zero": config.debug.set_pe_zero,
                      "set_histones_zero": config.debug.set_histones_zero, "set_extra_feat_zero": config.debug.set_extra_feat_zero,
                      "set_promoter_zero": config.debug.set_promoter_zero, "remove_ar": config.debug.remove_ar_events,
                      "one_tpm_ar": False}

    all_ds = pe_histone_dataset.PEHistoneDataset(data_folder= './data/', expr_type=expr_type, cell_type=cell,
                                                  n_extraFeat=n_extraFeat, usePromoterSignal=True,
                                                  n_enhancers=n_enhancers, hic_threshold=hic_threshold,
                                                  distance_threshold=distance_threshold, include_exons=config.optim.include_exons,
                                                  rna_seq_source=config.optim.rna_seq_source, tpm=config.optim.tpm_level,
                                                  single_event_train=config.optim.single_events, event_genes=config.optim.event_genes,
                                                  **ablation_tests)
    # create train, valid, test indices
    #train_idx, valid_idx, test_idx = create_set_indices(np.arange(len(all_ds)), train_ratio=0.8, valid_ratio=0.1, 
    #                                                    events=True, seed=42+int(fi))
    if config.optim.include_exons:
        train_idx, valid_idx, test_idx = split_multitask_ids(all_ds.event_keys, train_frac=0.8, val_frac=0.1, 
                                                             test_frac=0.1, seed=42+int(fi), short_run=config.debug.short_run, 
                                                             pre_splits=split_df, fold_i=fold_i, both_cell_lines=(cell == 'both'))
    else:
        train_idx, valid_idx, test_idx = create_set_indices(all_ds, train_ratio=0.8, valid_ratio=0.1,
                                                            events=False, seed=42+int(fi), splits=split_df, fold_i=fold_i)
    
    normals = None
    if config.optim.z_score_normalise:
        normals = all_ds.z_score_normalize(train_idx)
        print("Z-score normalization applied successfully.")

    train_ds = Subset(all_ds, train_idx)
    valid_ds = Subset(all_ds, valid_idx)
    test_ds = Subset(all_ds, test_idx)

    if use_pretrained:
        pretrained_convNet = enhancer_predictor_256bp()
        pt_model_name = f'{cell}_seq2activityLog2_leaveChrOut_combinedRS_2bins_bs64_H3K27ac_adamW_erisxdl_r0'
        checkpoint = torch.load(f"./trained_models/pretrained_enhancer_encoder/{fold_i}_best_{pt_model_name}_checkpoint.pt", 
                                weights_only=False, map_location=device)
        print('Loading pretrained model ...', pt_model_name)
        model = EPInformer_v2(n_encoder=n_encoder, pre_trained_encoder=pretrained_convNet.encoder,
                              n_enhancer=n_enhancers, out_dim=64, n_extraFeat=n_extraFeat, device=device, 
                              exon_data=config.optim.include_exons, separate_attention=True).to(device)
    else:
        model = EPInformer_v2(n_encoder=n_encoder, pre_trained_encoder=None, n_enhancer=n_enhancers, 
                              out_dim=64, n_extraFeat=n_extraFeat, device=device, exon_data=config.optim.include_exons, 
                              separate_attention=True).to(device)

    if config.optim.learn_loss_weights:
        print("Learning loss weights for the splicing and expression tasks.")
        weighted_loss = WeightedLoss().to(device)
    else:
        weighted_loss = None

    model = model.to(device)
    model.name = model.name.replace('EPInformerV2', config.base.model_type) + '.' +  cell + '.' + expr_type

    utils.train(model, train_ds, valid_dataset=valid_ds, epochs=n_epoch, model_name = model.name, fold_i=fi, 
                batch_size=batch_size, device=device, saved_model_path=saved_model_path, predict=expr_type, 
                loss_class=weighted_loss, weigh_samples=config.optim.weigh_samples)

    test_df = utils.test(model, test_ds, model_name = model.name, saved_model_path=saved_model_path, fold_i=fi, 
                         batch_size=batch_size, normals=normals, device=device, predict=config.base.expr_assay)
