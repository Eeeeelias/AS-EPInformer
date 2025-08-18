import argparse
from datetime import datetime
import random
from collections import defaultdict
import os
import shutil
import sys

from scipy import stats
from tqdm import tqdm
import pandas as pd
import numpy as np
import torch
import torch.utils.data as data_utils
from torch.utils.data import Subset, Dataset
import scripts.train_utils as reg_utils
import scripts.train_utils_binary as bin_utils
import scripts.setup_utils as sp
from EPInformer.models_multi import ASInformer, ASTransformer, ASLSTM, ASdCNNsmall
from scripts.pe_utils import plot_loss_curve
from sklearn.metrics import matthews_corrcoef, mean_squared_error, r2_score, roc_auc_score, balanced_accuracy_score
from sklearn.ensemble import RandomForestClassifier



def filter_psi(X, y):
    # filter out samples y is >0.2 and < 0.8
    mask = (y <= 0.2) | (y >= 0.8)
    y = y[mask]
    y = np.where(y > 0.5, 0, 1) 
    return X[mask], y


def loader_to_numpy(loader):
    X_list, y_list = [], []
    for data in loader:
        X, y, _ = data   # like your tuple structure
        X_list.append(X.view(X.size(0), -1).numpy())  # flatten images/features
        y_list.append(y.numpy())
    X = np.concatenate(X_list, axis=0)
    y = np.concatenate(y_list, axis=0)
    X, y = filter_psi(X,y)
    return X, y


##### init ######
parser = argparse.ArgumentParser()
parser.add_argument('--config_path', type=str, required=True)
parser.add_argument('--print_config_help', action='store_true', help='Print help info for config keys')

args = parser.parse_args()
yml_config = sp.parse_yaml_config(args.config_path)

if args.print_config_help:
    help_path = os.path.join(os.path.dirname(args.config_path), 'config_help.yaml')
    if os.path.exists(help_path):
        help_info = sp.parse_yaml_config(help_path)
        for key, val in help_info.items():
            print(f"{key}: {val}")
    else:
        print("No config_help.yaml found.")
    sys.exit(0)

config = sp.dict_to_namespace(yml_config)


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

if config.base.goal == 'regression':
    print("Training for regression task.")
    utils = reg_utils
else:
    print("Training for binary classification task.")
    utils = bin_utils

use_pretrained = config.hardware.use_pretrained_encoder
fold_list = sp.str_to_list(config.base.fold)
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

# split_df = pd.read_csv('./data/leave_chrom_out_crossvalidation_split_18377genes.csv', index_col=0)
split_df = None

for fi in fold_list:
    print("-"*10, 'fold', fi, '-'*10)
    fold_i = 'fold_' + str(fi)
    
    # removed pre-set indices for train, valid, test since our data is now event-based, not gene-based
    ablation_tests = {"set_exon_zero": config.debug.set_exon_zero, 
                      "set_pe_zero": config.debug.set_pe_zero,
                      "set_histones_zero": config.debug.set_histones_zero, 
                      "set_extra_feat_zero": config.debug.set_extra_feat_zero,
                      "set_promoter_zero": config.debug.set_promoter_zero, 
                      "remove_ar": config.debug.remove_ar_events,
                      "one_tpm_ar": False}

    dataset_options = {'data_folder': './data/',
                        'expr_type': expr_type,
                        'cell_type': cell,
                        'n_extraFeat': n_extraFeat,
                        'usePromoterSignal': True,
                        'n_enhancers': n_enhancers,
                        'hic_threshold': hic_threshold,
                        'distance_threshold': distance_threshold,
                        'include_exons': config.optim.include_exons,
                        'include_enhancers': config.optim.enhancer_histones,
                        'rna_seq_source': config.optim.rna_seq_source,
                        'tpm_level': config.optim.tpm_level,
                        'single_events': config.optim.single_events,
                        'event_genes': config.optim.event_genes}

    all_ds = sp.init_dataset(config.base.dataset, dataset_options, ablation_tests)

    # create train, valid, test indices
    #train_idx, valid_idx, test_idx = create_set_indices(np.arange(len(all_ds)), train_ratio=0.8, valid_ratio=0.1, 
    #                                                    events=True, seed=42+int(fi))
    if config.optim.include_exons and config.base.goal != "binary":
        train_idx, valid_idx, test_idx = sp.split_multitask_ids(all_ds.event_keys, train_frac=0.8, val_frac=0.1, 
                                                             test_frac=0.1, seed=42+int(fi), short_run=config.debug.short_run, 
                                                             pre_splits=split_df, fold_i=fold_i, both_cell_lines=(cell == 'both'))
    elif config.optim.include_exons:
        train_idx, valid_idx, test_idx = sp.split_binary(all_ds.event_keys, train_frac=0.8, val_frac=0.1,test_frac=0.1, 
                                                         seed=42+int(fi), short_run=config.debug.short_run, 
                                                         split_simple=True)
    else:
        train_idx, valid_idx, test_idx = sp.create_set_indices(all_ds, train_ratio=0.8, valid_ratio=0.1,
                                                            events=False, seed=42+int(fi), splits=split_df, fold_i=fold_i)
    
    normals = None
    if config.optim.z_score_normalise:
        normals = all_ds.z_score_normalize(train_idx)
        print("Z-score normalization applied successfully.")

    train_ds = Subset(all_ds, train_idx)
    trainloader = data_utils.DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=5, pin_memory=True)
    valid_ds = Subset(all_ds, valid_idx)
    valloader = data_utils.DataLoader(valid_ds, batch_size=batch_size, shuffle=False, num_workers=5, pin_memory=True)
    test_ds = Subset(all_ds, test_idx)
    testloader = data_utils.DataLoader(test_ds, batch_size=batch_size, shuffle=False, num_workers=5, pin_memory=True)

    # use a different model
    model = RandomForestClassifier(n_estimators=1000, random_state=42)

    # print number of parameters of the model
    # RandomForestClassifier does not have parameters in the same way as PyTorch models
    print(f"Number of trees in the forest: {model.n_estimators}")

    X_train, y_train = loader_to_numpy(trainloader)
    X_val, y_val     = loader_to_numpy(valloader)
    X_test, y_test   = loader_to_numpy(testloader)

    model.fit(X_train, y_train)

    val_preds = model.predict(X_val)
    val_acc = balanced_accuracy_score(y_val, val_preds)
    print(f"Validation accuracy: {val_acc:.4f}")

    # test
    test_preds = model.predict(X_test)
    test_acc = balanced_accuracy_score(y_test, test_preds)
    test_mcc = matthews_corrcoef(y_test, test_preds)
    test_auroc = roc_auc_score(y_test, test_preds)
    print(f"Test accuracy: {test_acc:.4f}")
    print(f"Test MCC: {test_mcc:.4f}")
    print(f"Test AUROC: {test_auroc:.4f}")