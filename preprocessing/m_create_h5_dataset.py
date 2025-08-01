import os
import sys

import pandas as pd
import numpy as np
import h5py as h5
import tqdm 
# Suppress stderr temporarily
stderr_fileno = sys.stderr.fileno()
with open(os.devnull, 'w') as fnull:
    old_stderr = os.dup(stderr_fileno)
    os.dup2(fnull.fileno(), stderr_fileno)
    import torch  # import after suppressing
    from torch.nn.functional import one_hot, pad
    os.dup2(old_stderr, stderr_fileno)  # Restore stderr
    os.close(old_stderr)

VOCAB = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def one_hot_encode_event(sequences):
    upstream, downstream, exon = None, None, None
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
    segment_tensor = torch.stack([one_hot_encode(upstream, VOCAB), one_hot_encode(exon, VOCAB), 
                                one_hot_encode(downstream, VOCAB)])

    return segment_tensor


def cut_to_len(seq, length):
    if len(seq) == length:
        return seq
    elif len(seq) > length:
        # get the center slice
        start = (len(seq) - length) // 2
        return seq[start:start + length]
    else:
        # pad with zeros at both ends
        pad_length = length - len(seq)
        if len(seq) % 2 == 0:        
            pad_start = pad_length // 2
            pad_end = pad_length - pad_start
        else: # just so it works with uneven lengths
            pad_start = (pad_length // 2) + 1
            pad_end = pad_length - pad_start 
        return np.pad(seq, ((pad_start, pad_end), (0, 0)), mode='constant', constant_values=0)


def one_hot_encode(seq, vocab, length=1024):
    indices = [vocab.get(item, -1) for item in seq]
    tensor = torch.tensor(indices)
    valid_mask = tensor >= 0
    one_hot_tensor = torch.zeros((len(tensor), len(vocab)), dtype=torch.float)
    one_hot_tensor[valid_mask] = one_hot(tensor[valid_mask], num_classes=len(vocab)).float() # pylint: disable=not-callable
    if not length:
        return one_hot_tensor
    # add padding if length is specified
    if len(seq) < length:
        one_hot_tensor = torch.cat([one_hot_tensor, torch.zeros(length - len(seq), 4)], dim=0)
    elif len(seq) > length:
        one_hot_tensor = one_hot_tensor[:length]
    return one_hot_tensor


def adjust_histone_info(histone_info, split_segments=True):
    segment_order = ['intron1', 'exon', 'intron2']
    if split_segments:
        histone_info['event_id'] = histone_info['gene_id'].apply(lambda x: x.split('|')[1])
        histone_info['segment'] = histone_info['gene_id'].apply(lambda x: x.split('|')[0])
        histone_info['segment'] = pd.Categorical(histone_info['segment'], categories=segment_order, ordered=True)
        return histone_info.set_index('event_id')
    else:
        return histone_info.set_index('gene_id')


def create_gene_h5(gene_h5, current_h5, enhancer_links, enhancer_histone_info):
    # existing keys: ['activity', 'distance', 'ensid', 'hic', 'pe_code']
    # enhancer_links columns: chr, start, end, name, geneid, tss, sequence
    # enhancer_histone_info: dataframe cols: h3k27ac, h3k24me3, rows 
    dt = h5.string_dtype(encoding='utf-8')
    gene_ids = []
    pe_seqs = []
    histone_tensors = {key.split(".")[0]: [] for key in enhancer_histone_info.columns if key not in ['covered', 'gene_id']}
    for idx, gene in tqdm.tqdm(enumerate(current_h5['ensid']), total=len(current_h5['ensid'])):
        gene_id = gene.decode('utf-8')
        gene_ids.append(gene_id)
        promoter_seq = current_h5['pe_code'][idx][0].astype(np.int16).reshape(1, 2000, 4)
        # get all enhancer links for this gene
        gene_links = enhancer_links[enhancer_links['gene_id'] == gene_id]

        # quickly fill with zero if no links
        if gene_links.empty:
            pe_one_hot = np.concatenate([promoter_seq, np.zeros((60, 2000, 4), dtype=np.int16)], axis=0)
            pe_seqs.append(pe_one_hot)
            for col in enhancer_histone_info.columns:
                if col in ['covered', 'gene_id']:
                    continue
                histone_tensors[col.split(".")[0]].append(torch.zeros(60))
            continue

        # first get the sequences
        # we need np array of shape (61, 2000, 4)
        enhancer_seqs = []
        for _, row in gene_links.iterrows():
            enhancer_seq = row['sequence']
            enhancer_one_hot = one_hot_encode(enhancer_seq, VOCAB, length=None) # we want the full sequence here
            enhancer_one_hot = cut_to_len(enhancer_one_hot, 2000)
            enhancer_seqs.append(enhancer_one_hot)
        
        enhancer_seqs = np.stack(enhancer_seqs)
        pe_one_hot = np.concatenate([promoter_seq, enhancer_seqs], axis=0)
        if pe_one_hot.shape[0] < 61:
            padding = np.zeros((61 - pe_one_hot.shape[0], 2000, 4), dtype=np.int16)
            pe_one_hot = np.concatenate([pe_one_hot, padding], axis=0)
        assert pe_one_hot.shape == (61, 2000, 4), f"Expected shape (61, 2000, 4), got {pe_one_hot.shape}"
        pe_seqs.append(pe_one_hot)

        # then do histone information
        histone_rows = enhancer_histone_info.loc[enhancer_histone_info.index.isin(list(gene_links['name']))]
        for col in enhancer_histone_info.columns:
            if col in ['covered', 'gene_id']:
                continue
            histone_tensor = torch.tensor(histone_rows[col].values, dtype=torch.float32)
            if len(histone_tensor) < 60: # pad if needed
                histone_tensor = pad(histone_tensor, (0, 60 - histone_tensor.size(0)))
            histone_tensors[col.split(".")[0]].append(histone_tensor)

    # convert list of numpy arrays to a single tensor
    pe_seqs = np.stack(pe_seqs)
    pe_seqs = torch.tensor(pe_seqs, dtype=torch.bool) # save as bool since its more memory efficient
    histone_tensors = {key: torch.stack(value) for key, value in histone_tensors.items()}

    gene_h5.create_dataset('gene_id', data=gene_ids, dtype=dt)
    gene_h5.create_dataset('pe_seqs',  data=pe_seqs.numpy(), dtype='bool')
    gene_h5.create_dataset('distance', data=current_h5['distance'][:])
    for key, value in histone_tensors.items():
        gene_h5.create_dataset(key, data=value, dtype='float32')
    return gene_h5


def create_event_h5(event_h5, gene_sequences, histone_info):
    gene_ids = []
    event_ids = list(gene_sequences.keys())
    one_hot_seqs = []
    histone_tensors_k = {key.split(".")[0]: [] for key in histone_info['K562'].columns if key not in ['covered', 'segment', 'gene_id']}
    histone_tensors_g = {key.split(".")[0]: [] for key in histone_info['GM12878'].columns if key not in ['covered', 'segment', 'gene_id']}
    for event_id in tqdm.tqdm(gene_sequences.keys(), total=len(gene_sequences)):
        gene_id = event_id.split(";")[0]
        gene_ids.append(gene_id)

        seqs = gene_sequences[event_id]
        one_hot_seqs.append(one_hot_encode_event(seqs))

        for cell in histone_info.keys():
            if event_id not in histone_info[cell].index:
                print("Skipping event_id", event_id, "as it does not exist in histone info")
                continue
            histone_rows = histone_info[cell].loc[event_id].sort_values('segment')
            for col in histone_rows.columns:
                if col in ['covered', 'segment', 'gene_id']:
                    continue
                histone_tensor = torch.tensor(histone_rows[col].values, dtype=torch.float32)
                if cell == 'K562':
                    histone_tensors_k[col.split(".")[0]].append(histone_tensor)
                else:
                    histone_tensors_g[col.split(".")[0]].append(histone_tensor)

    # combine back to tensors
    one_hot_seqs = torch.stack(one_hot_seqs)
    histone_tensors_k = {key: torch.stack(value) for key, value in histone_tensors_k.items()}
    histone_tensors_g = {key: torch.stack(value) for key, value in histone_tensors_g.items()}
    
    dt = h5.string_dtype(encoding='utf-8')
    event_h5.create_dataset('gene_id', data=gene_ids, dtype=dt)
    event_h5.create_dataset('event_id', data=event_ids, dtype=dt)
    event_h5.create_dataset('event_seq', data=one_hot_seqs, dtype='float32')
    for key, value in histone_tensors_k.items():
        event_h5.create_dataset(f'{key}_K562', data=value, dtype='float32')
    for key, value in histone_tensors_g.items():
        event_h5.create_dataset(f'{key}_GM12878', data=value, dtype='float32')
    
    return event_h5

    # h5[encoded seq] h5[histone col][cell line] h5[event id]


if __name__ == "__main__":
    cell_line = 'K562'

    # load all the h5 files
    if cell_line == 'K562':
        current_h5 = h5.File(f'../data/{cell_line}_DNase_ENCFF257HEE_2kb_4DNFITUOMFUQ_enhancer_promoter_encoding.h5', 'r')
    else:
        current_h5 = h5.File(f'../data/{cell_line}_DNase_ENCFF020WZB_2kb_4DNFI1UEG1HD_promoter_enhancer_encoding.h5', 'r')
    new_gene_h5 = h5.File(f'../data/{cell_line}_histone_appended_pe_encoding.h5', 'w')
    # event_h5 = h5.File(f'../data/event_encoding.h5', 'w')
    #gene_sequences = h5.File('../data/event_sequences.h5', 'r')

    # load histone info
    try:
        event_histone_info = {
            'K562': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Events_Combined_K562_Histone_Signals.csv'),
            'GM12878': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Events_Combined_GM12878_Histone_Signals.csv')
        }
        enhancer_histone_info = {
            'K562': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Enhancer_Combined_K562_Histone_Signals.csv'),
            'GM12878': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Enhancer_Combined_GM12878_Histone_Signals.csv')
        }
    except FileNotFoundError as exc:
        raise FileNotFoundError("Histone info files not found. Please run m_read_histone_signals_events.py " \
                                "and m_read_histone_signals_promoter.py first.") from exc

    # get the gene enhancer links
    try:
        enhancer_links = pd.read_csv(f'../data/{cell_line}_gene_enhancer_links_with_sequences.csv')
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Gene enhancer links file for {cell_line} not found." \
                                 "Please run m_create_gene_enhancer_links.py first.") from exc

    event_histone_info = {key: adjust_histone_info(df) for key, df in event_histone_info.items()}
    enhancer_histone_info = {key: adjust_histone_info(df, split_segments=False) for key, df in enhancer_histone_info.items()}

    create_gene_h5(new_gene_h5, current_h5, enhancer_links, enhancer_histone_info[cell_line])
    # create_event_h5(event_h5, gene_sequences, event_histone_info)

    current_h5.close()
    new_gene_h5.close()
    # event_h5.close()
    #gene_sequences.close()