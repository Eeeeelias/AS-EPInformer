import pandas as pd
import pysam
import numpy as np
from tqdm import tqdm
import h5py as h5
from Bio.Seq import Seq
import pyBigWig

from m_read_histone_signals_events import file_uuids, get_file_name
from m_extract_bp_bigwig import resize_seq


psi_file = pd.read_csv('../data/transcript_SE_f1.psi', sep='\t')
fasta_file = pysam.FastaFile('/nfs/data/references/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')

# I need to get sequence info & histone mark info

def get_event_data(event_id):
    # split event_id ENSG00000000419.12;SE:chr20:50940933-50941129:50941209-50942031:-
    gene_id = event_id.split(";")[0].split(".")[0]
    chrom = event_id.split(":")[1]
    strand = event_id.split(":")[-1]
    intron_seq1 = (event_id.split(":")[2].split("-")[0], event_id.split(":")[2].split("-")[1])
    exon_seq = (event_id.split(":")[2].split("-")[1], event_id.split(":")[3].split("-")[0])
    intron_seq2 = (event_id.split(":")[3].split("-")[0], event_id.split(":")[3].split("-")[1])
    return {'event': event_id, 'gene_id': gene_id, 'chrom': chrom, 'strand': strand, 
            'intron_seq1': intron_seq1, 'exon_seq': exon_seq, 'intron_seq2': intron_seq2}


def get_sequence(event_info):
    chromosome = event_info['chrom']
    strand = event_info['strand']
    regions = [event_info['intron_seq1'], event_info['exon_seq'], event_info['intron_seq2']]
    region_names = ['upstream', 'exon', 'downstream']
    sequences = {}
    contains_n = False
    is_poisoned = False
    for region, name in zip(regions, region_names):
        start = int(region[0])
        end = int(region[1])
        if name != 'exon':
            end -= 1
        else:
            start -= 1
        seq = fasta_file.fetch(chromosome, start, end)

        if strand == '-':
            seq = str(Seq(seq).reverse_complement())
        if "N" in seq:
            contains_n = True
        sequences[f"{region}:{name}"] = seq
    if contains_n:
        is_poisoned = True

    return {'sequences': sequences, 'is_poisoned': is_poisoned}


def one_hot_encode_event(sequences):
    upstream, downstream, exon = None, None, None
    for key in sequences.keys():
        if 'upstream' in key:
            upstream = sequences[key]
            upstream = upstream[-1024:]
        elif 'downstream' in key:
            downstream = sequences[key]
            downstream = downstream[:1024]
        elif 'exon' in key:
            exon = sequences[key]
            if len(exon) > 1024:
                exon_start = exon[:512]
                exon_end = exon[-512:]
                exon = exon_start + exon_end

    return {'upstream': one_hot_encode(upstream), 'exon': one_hot_encode(exon), 'downstream': one_hot_encode(downstream)}


def one_hot_encode(sequence, seq_len=1024):
    vocab = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    one_hot = np.zeros((seq_len, len(vocab)), dtype=int)
    for i, nucleotide in enumerate(sequence):
        if nucleotide in vocab:
            one_hot[i, vocab[nucleotide]] = 1
    return one_hot


def extract_bp_histone_signals(bigwig_file, full_event_info, seq_len=1024):
    bw = pyBigWig.open(bigwig_file)
    outputs = {}

    for event in full_event_info:
        # {'event': event_id, 'gene_id': gene_id, 'chrom': chrom, 'strand': strand, 
        #    'intron_seq1': intron_seq1, 'exon_seq': exon_seq, 'intron_seq2': intron_seq2}
        chrom = event['chrom']
        event_id = event['event']
        for part_name in ['intron_seq1', 'exon_seq', 'intron_seq2']:
            part = event[part_name]
            start = int(part[0])
            end = int(part[1])
            values = bw.values(chrom, start, end)
            values = resize_seq(values, seq_len)
            outputs[f"{part_name.replace('_seq', '')}|{event_id}"] = values

    bw.close()
    return outputs


def extract_5mc_signal(bigwig_file, full_event_info, seq_len=1024):
    bw = pyBigWig.open(bigwig_file)
    outputs = {}

    for event in full_event_info:
        chrom = event['chrom']
        event_id = event['event']
        for part_name in ['intron_seq1', 'exon_seq', 'intron_seq2']:
            part = event[part_name]
            start = int(part[0])
            end = int(part[1])
            values = bw.values(chrom, start, end)
            values = resize_seq(values, seq_len)
            outputs[f"{part_name.replace('_seq', '')}|{event_id}"] = values

    bw.close()
    return outputs


def get_all_histone_marks(uuids, cell_line, full_event_info):
    all_histone_marks = {}
    for histone_type, tissues in tqdm(uuids.items()):
        for tissue, uuid in tissues.items():
            if tissue != cell_line:
                continue
            bigwig_file = get_file_name(uuid)
            histone_marks = extract_bp_histone_signals(bigwig_file, full_event_info)
            all_histone_marks[histone_type] = histone_marks
    return all_histone_marks


def transform_histones(all_histone_marks):
    # for each event we should have histone marks for intron1, exon, intron2
    # go from mark1 1024, mark2 1024, ... to intron1|event (1024, 6)
    all_events = {}
    for histone_type, marks in all_histone_marks.items():
        for event_id, mark in marks.items():
            if event_id not in all_events:
                all_events[event_id] = {}
            all_events[event_id][histone_type] = mark
    # combine the marks in numpy so we have (1024, 6)
    all_events_adj = {e_id: np.stack([mark for _, mark in all_events[e_id].items()], axis=-1) 
                      for e_id in all_events.keys()}
    return all_events_adj



def assemble_event(event_seqs, event_hist):
    # stack sequences from event [(1024, 4), (1024, 4), (1024, 4)] -> (3, 1024, 4)
    stacked_seqs = np.stack([event_seqs['upstream'], event_seqs['exon'], event_seqs['downstream']], axis=0)
    stacked_hist = np.stack(event_hist, axis=0)
    stacked_all = np.concatenate((stacked_seqs, stacked_hist), axis=-1)  # (3, 1024, 4 + 6)
    return stacked_all


def create_h5_ds(h5_file, dataset, names):
    dt = h5.string_dtype(encoding='utf-8')
    with h5.File(h5_file, 'w') as f:
        f.create_dataset('event_data', data=dataset, compression='gzip', compression_opts=5, chunks=(1, 3, 1024, 10))
        f.create_dataset('event_id', data=names, dtype=dt)

        print(f"Dataset created in {h5_file} with {len(names)} events.")


def create_response(response_file, psi_events, cell_line):
    with open(response_file, 'w') as f:
        f.write(f"event_id,gene_id,{cell_line}_SE_psi\n")
        for _, event in psi_events.iterrows():
            psi_value = event[f"{cell_line}_TPM"]
            gene_id = event['event_id'].split(";")[0].split(".")[0]
            f.write(f"{event['event_id']},{gene_id},{psi_value}\n")


if __name__ == "__main__":
    cell_line = 'GM12878'
    h5_path = f'../data/extended_events_{cell_line}.h5'
    response_path = f"../data/extended_psi_response_{cell_line}.csv"
    # filter out psi where K562_TPM is nan
    print(f"{len(psi_file)} events in total")
    psi_file = psi_file.reset_index()
    psi_file.rename(columns={'index': 'event_id'}, inplace=True)
    psi_file = psi_file[~psi_file[f"{cell_line}_TPM"].isna()]
    print(f"{len(psi_file)} events after filtering for {cell_line} TPM")
    # get event data
    psi_events = psi_file['event_id'].apply(get_event_data)

    event_seqs = {}
    for event in tqdm(psi_events):
        event_info = get_sequence(event)
        if event_info['is_poisoned']:
            continue
        event_seqs[event['event']] = one_hot_encode_event(event_info['sequences'])
    print(len(event_seqs), "events after filtering for poisoned sequences (i.e. containing N)")

    histone_marks = transform_histones(get_all_histone_marks(file_uuids, cell_line, psi_events))
    print(int(len(histone_marks) / 3), "histone marks extracted")

    event_ds = np.zeros((len(event_seqs), 3, 1024, 4 + 6))
    event_names = []
    for i, event_id in enumerate(event_seqs.keys()):
        hist_marks = [histone_marks[f"intron1|{event_id}"], histone_marks[f"exon|{event_id}"], histone_marks[f"intron2|{event_id}"]]
        event_seq = event_seqs[event_id]
        assembled_event = assemble_event(event_seq, hist_marks)
        event_ds[i] = assembled_event
        event_names.append(event_id)

    create_response(response_path, psi_file, cell_line)
    create_h5_ds(h5_path, event_ds, event_names)