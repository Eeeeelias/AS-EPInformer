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
from torch.utils.data import Subset, Dataset
import scripts.train_utils as utils
import scripts.promoter_enhancer_dataset as pe_dataset
import scripts.pe_histone_dataset as pe_histone_dataset
import scripts.setup_utils as sp

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
        train_idx, valid_idx, test_idx = sp.split_multitask_ids(all_ds.event_keys, train_frac=0.8, val_frac=0.1, 
                                                             test_frac=0.1, seed=42+int(fi), short_run=config.debug.short_run, 
                                                             pre_splits=split_df, fold_i=fold_i, both_cell_lines=(cell == 'both'))
    else:
        train_idx, valid_idx, test_idx = sp.create_set_indices(all_ds, train_ratio=0.8, valid_ratio=0.1,
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
