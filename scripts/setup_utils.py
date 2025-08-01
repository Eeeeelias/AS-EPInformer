import os
import random
import yaml
from types import SimpleNamespace
from collections import defaultdict

import numpy as np
import pandas as pd

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

def init_dataset(dataset_class, input_config, ablation_tests=None):
    dataset_args = {
        'data_folder': input_config['data_folder'],
        'expr_type': input_config['expr_type'],
        'cell_type': input_config['cell_type'],
        'n_extraFeat': input_config['n_extraFeat'],
        'usePromoterSignal': input_config['usePromoterSignal'],
        'n_enhancers': input_config['n_enhancers'],
        'hic_threshold': input_config['hic_threshold'],
        'distance_threshold': input_config['distance_threshold'],
        'include_exons': input_config['include_exons'],
        'include_enhancers': input_config.get('include_enhancers', True),
        'rna_seq_source': input_config['rna_seq_source'],
        'tpm': input_config['tpm_level'],
        'single_event_train': input_config['single_events'],
        'event_genes': input_config['event_genes'],
    }

    # Merge in ablation tests if any
    if ablation_tests:
        dataset_args.update(ablation_tests)

    if dataset_class == 'pe_histone':
        cls = pe_histone_dataset.PEHistoneDataset
    else:
        cls = pe_dataset.PromoterEnhancerDataset

    return cls(**dataset_args)
