import argparse
from datetime import datetime
import random
from collections import defaultdict

from scipy import stats
from tqdm import tqdm
import pandas as pd
import numpy as np
import torch
from torch.utils.data import Subset, Dataset
import scripts.utils_forTraining as utils
import scripts.promoter_enhancer_dataset as pe_dataset

parser = argparse.ArgumentParser()
def list_of_strings(arg):
    return arg.split(',')
parser.add_argument('--cell', type=str, help='cell line (support K562 and GM12878)', choices=['K562', 'GM12878'])  
parser.add_argument("--fold", type=list_of_strings, help="test fold", default='all')
parser.add_argument("--model_type", type=str, help='EPInformer type', default='EPInformer-PE-Activity', 
                    choices=['EPInformer-PE', 'EPInformer-PE-Activity', 'EPInformer-PE-Activity-HiC'])  
parser.add_argument('--distance_threshold', type=int, help='max distance to TSS', default=100_000) 
parser.add_argument('--hic_threshold', type=int, help='hic loop thresold', default=-1) 
parser.add_argument('--expr_assay', type=str, help='expression_assay', choices=['CAGE', 'RNA', 'multi', 'splice'])
parser.add_argument('--batch_size', type=int, help='batch size', default=16)
parser.add_argument('--n_interact_enc',type=int, help='layers of interaction encoder', default=3)
parser.add_argument('--epochs',type=int, help='training epochs', default=100)
parser.add_argument('--cuda', help='use cuda', action='store_true')
parser.add_argument('--xpu', help='use xpu', action='store_true')
parser.add_argument('--use_pretrained_encoder', help='use pretrained sequence encoder', action='store_true')
parser.add_argument('--rna_seq_source', type=str, help='Which RNA-seq source to use', choices=['xpresso', 'epiatlas'], default='xpresso')
parser.add_argument('--tpm_level', type=str, help='TPM level for RNA-seq', choices=['gene', 'transcript'], default='gene')
parser.add_argument('--include_exons', help='Include exons in the input data', action='store_true')
parser.add_argument('--single_events', help='Use single events per gene for training', action='store_true')
parser.add_argument('--z_score_normalise', help='Apply z-score normalization to the training data', action='store_true')
parser.add_argument('--event_genes', action='store_true', help='Use only genes that also have events in the training set')
parser.add_argument('--learn_loss_weights', action='store_true', help='Learn loss weights for the splicing and expression tasks')
parser.add_argument('--weigh_samples', action='store_true', help='Weigh samples based on their frequency in the training set')
parser.add_argument('--expr_loss', type=str, default='mse', choices=['mse', 'smoothl1'], help='Loss function for expression task')
parser.add_argument('--splice_loss', type=str, default='bce', choices=['bce', 'smoothl1'], help='Loss function for splicing task')
# ablation/debug arguments
parser.add_argument('--short_run', action='store_true', help='Run a short version for testing purposes')
parser.add_argument('--set_exon_zero', action='store_true', help='Set exon data to zero for testing purposes')
parser.add_argument('--set_pe_zero', action='store_true', help='Set promoter/enhancer data to zero for testing purposes')
parser.add_argument('--set_rna_zero', action='store_true', help='Set RNA half life data to zero for testing purposes')
parser.add_argument('--set_extra_feat_zero', action='store_true', help='Set extra features (HiC, Activity, Distance) to zero for testing purposes')
parser.add_argument('--set_promoter_zero', action='store_true', help='Set promoter data to zero for testing purposes')

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
                        seed: int = 42, short_run=False) -> tuple[list[str], list[str], list[str]]:
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

    if short_run:
        # Limit the number of samples for a short run
        train_indices = train_indices[:2000]
        val_indices = val_indices[:100]
        test_indices = test_indices[:100]

    return train_indices, val_indices, test_indices

# example
# python train_EPInformer.py --cell K562  --model_type EPInformer-PE-Activity --expr_assay CAGE 
# --use_pretrained_encoder --batch_size 16

##### parameter ######
args = parser.parse_args()

#### import the right EPinformer model
if args.expr_assay == 'multi' or args.expr_assay == 'splice' or args.include_exons:
    from EPInformer.models_multi import EPInformer_v2, enhancer_predictor_256bp, WeightedLoss
else:
    from EPInformer.models import EPInformer_v2, enhancer_predictor_256bp, WeightedLoss


cell = args.cell

if args.cuda:
    device = 'cuda'
elif args.xpu:
    device = 'xpu'
else:
    device = 'cpu'
distance_threshold = args.distance_threshold
n_epoch = args.epochs
hic_threshold = args.hic_threshold
if hic_threshold == -1:
    hic_threshold = None

if args.model_type == 'EPInformer-PE': 
    n_extraFeat = 1
elif args.model_type == 'EPInformer-PE-Activity':
    n_extraFeat = 2
elif args.model_type == 'EPInformer-PE-Activity-HiC':
    n_extraFeat = 3
else:
    raise ValueError(f"Unsupported model type: {args.model_type}")

use_pretrained = args.use_pretrained_encoder
fold_list = args.fold 
n_encoder = args.n_interact_enc
batch_size = args.batch_size 
expr_type = args.expr_assay
n_enhancers = 60
#################

today = datetime.now()   # Get date

datetime_str = today.strftime("%Y-%m-%d-%H-%M")
split_df = pd.read_csv('./data/leave_chrom_out_crossvalidation_split_18377genes.csv', index_col=0)
if args.short_run:
    saved_model_path = None
else:
    saved_model_path = f'./trained_models/{datetime_str}/'

if 'all' in fold_list:
    fold_list = list(range(1, 13))

if args.short_run:
    fold_list = fold_list[:1]  # For testing, only use the first fold

for fi in fold_list:
    print("-"*10, 'fold', fi, '-'*10)
    fold_i = 'fold_' + str(fi)
    
    # removed pre-set indices for train, valid, test since our data is now event-based, not gene-based

    all_ds = pe_dataset.promoter_enhancer_dataset(data_folder= './data/', expr_type=expr_type, cell_type=cell,
                                                  n_extraFeat=n_extraFeat, usePromoterSignal=True,
                                                  n_enhancers=n_enhancers, hic_threshold=hic_threshold,
                                                  distance_threshold=distance_threshold, include_exons=args.include_exons,
                                                  rna_seq_source=args.rna_seq_source, tpm=args.tpm_level,
                                                  single_event_train=args.single_events, event_genes=args.event_genes,
                                                  set_exon_zero=args.set_exon_zero, set_pe_zero=args.set_pe_zero,
                                                  set_rna_zero=args.set_rna_zero, set_extra_feat_zero=args.set_extra_feat_zero,
                                                  set_promoter_zero=args.set_promoter_zero)
    # create train, valid, test indices
    #train_idx, valid_idx, test_idx = create_set_indices(np.arange(len(all_ds)), train_ratio=0.8, valid_ratio=0.1, 
    #                                                    events=True, seed=42+int(fi))
    if args.include_exons:
        train_idx, valid_idx, test_idx = split_multitask_ids(all_ds.event_keys, train_frac=0.8, val_frac=0.1, 
                                                             test_frac=0.1, seed=42+int(fi), short_run=args.short_run)
    else:
        train_idx, valid_idx, test_idx = create_set_indices(all_ds, train_ratio=0.8, valid_ratio=0.1,
                                                            events=False, seed=42+int(fi), splits=split_df, fold_i=fold_i)
    
    normals = None
    if args.z_score_normalise:
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
                              exon_data=args.include_exons, separate_attention=True).to(device)
    else:
        model = EPInformer_v2(n_encoder=n_encoder, pre_trained_encoder=None, n_enhancer=n_enhancers, 
                              out_dim=64, n_extraFeat=n_extraFeat, device=device, exon_data=args.include_exons, 
                              separate_attention=True).to(device)

    if args.learn_loss_weights:
        print("Learning loss weights for the splicing and expression tasks.")
        weighted_loss = WeightedLoss().to(device)
    else:
        weighted_loss = None

    model = model.to(device)
    model.name = model.name.replace('EPInformerV2', args.model_type) + '.' +  cell + '.' + expr_type

    utils.train(model, train_ds, valid_dataset=valid_ds, epochs=n_epoch, model_name = model.name, fold_i=fi, 
                batch_size=batch_size, device=device, saved_model_path=saved_model_path, predict=expr_type, 
                loss_class=weighted_loss, weigh_samples=args.weigh_samples)
    
    test_df = utils.test(model, test_ds, model_name = model.name, saved_model_path=saved_model_path, fold_i=fi, 
                         batch_size=batch_size, normals=normals, device=device, predict=args.expr_assay)
