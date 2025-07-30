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
    from torch.nn.functional import one_hot
    os.dup2(old_stderr, stderr_fileno)  # Restore stderr
    os.close(old_stderr)



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
    vocab = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    segment_tensor = torch.stack([one_hot_encode(upstream, vocab), one_hot_encode(exon, vocab), 
                                one_hot_encode(downstream, vocab)])

    return segment_tensor


def one_hot_encode(seq, vocab, length=1024):
    indices = [vocab[item] for item in seq]
    tensor = torch.tensor(indices)
    one_hot_tensor = one_hot(tensor, num_classes=len(vocab)).float() # pylint: disable=not-callable
    # add padding
    if len(seq) < length:
        one_hot_tensor = torch.cat([one_hot_tensor, torch.zeros(length - len(seq), 4)], dim=0)
    elif len(seq) > length:
        one_hot_tensor = one_hot_tensor[:length]
    return one_hot_tensor


def adjust_histone_info(histone_info):
    segment_order = ['intron1', 'exon', 'intron2']
    histone_info['event_id'] = histone_info['gene_id'].apply(lambda x: x.split('|')[1])
    histone_info['segment'] = histone_info['gene_id'].apply(lambda x: x.split('|')[0])
    histone_info['segment'] = pd.Categorical(histone_info['segment'], categories=segment_order, ordered=True)
    return histone_info.set_index('event_id')


def create_gene_h5(gene_h5, current_h5, enhancer_links, enhancer_histone_info):
    # existing keys: ['activity', 'distance', 'ensid', 'hic', 'pe_code']
    dt = h5.string_dtype(encoding='utf-8')
    # TODO: first need to get the sequences for the enhancer regions
    pass


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
    current_h5 = h5.File(f'../data/{cell_line}_DNase_ENCFF257HEE_2kb_4DNFITUOMFUQ_enhancer_promoter_encoding.h5', 'r')
    new_gene_h5 = h5.File(f'../data/{cell_line}_histone_appended_pe_encoding.h5', 'w')
    event_h5 = h5.File(f'../data/{cell_line}_event_encoding.h5', 'w')
    gene_sequences = h5.File('../data/event_sequences.h5', 'r')

    # load histone info
    try:
        event_histone_info = {
            'K562': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Events_Combined_K562_Histone_Signals.csv'),
            'GM12878': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Events_Combined_GM12878_Histone_Signals.csv')
        }
        enhancer_histone_info = {
            'K562': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Enhancers_Combined_K562_Histone_Signals.csv'),
            'GM12878': pd.read_csv('../data/IHEC-ChIP-Seq-Histone-Signals/Enhancers_Combined_GM12878_Histone_Signals.csv')
        }
    except FileNotFoundError as exc:
        raise FileNotFoundError("Histone info files not found. Please run m_read_histone_signals_events.py" \
                                "and m_read_histone_signals_promoter.py first.") from exc

    # get the gene enhancer links
    try:
        enhancer_links = pd.read_csv(f'../data/{cell_line}_gene_enhancer_links.csv')
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Gene enhancer links file for {cell_line} not found." \
                                 "Please run m_create_gene_enhancer_links.py first.") from exc

    event_histone_info = {key: adjust_histone_info(df) for key, df in event_histone_info.items()}
    enhancer_histone_info = {key: adjust_histone_info(df) for key, df in enhancer_histone_info.items()}

    create_gene_h5(new_gene_h5, current_h5, enhancer_links, enhancer_histone_info[cell_line])
    create_event_h5(event_h5, gene_sequences, event_histone_info)

    current_h5.close()
    new_gene_h5.close()
    event_h5.close()
    gene_sequences.close()