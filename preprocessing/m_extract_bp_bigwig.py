import pyBigWig
import pandas as pd
import h5py as h5
import numpy as np

from m_read_histone_signals_events import file_uuids, get_file_name


def process_histone_signals(bed_file, existing_h5=None, proc_events=True, only_cell=None):
    # need a dict that per tissue (cell line) holds a np array per histone signal
    dataset_dict = {cell: {} for cell in file_uuids[list(file_uuids.keys())[0]].keys()}

    for histone_type, tissues in file_uuids.items():
        for tissue, uuid in tissues.items():
            if only_cell and tissue != only_cell:
                print(f"Skipping {tissue}...")
                continue
            bigwig_file = get_file_name(uuid)
            print(f"Processing {histone_type} for {tissue} ({'events' if proc_events else 'enhancer'})...")
            if proc_events:
                hist_outs = extract_bp_histone_signals(bigwig_file, bed_file, proc_events, seq_len=1024)
                hist_arr = combine_events(hist_outs)
            else:
                hist_outs = extract_bp_histone_signals(bigwig_file, bed_file, proc_events, seq_len=2000)
                hist_arr = combine_enhancers(hist_outs, existing_h5)
            dataset_dict[tissue][histone_type] = hist_arr

    return dataset_dict


def resize_seq(seq, length):
    if len(seq) < length:
        seq = np.pad(seq, (0, length - len(seq)), mode='constant')
    elif len(seq) > length:
        seq = seq[:length]
    return seq


def resize_seq_center(seq, length):
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
        seq = np.array(seq)
        return np.pad(seq, (pad_start, pad_end), mode='constant', constant_values=0)


def extract_bp_histone_signals(bigwig_file, bed_file, proc_events, seq_len=1024):
    bw = pyBigWig.open(bigwig_file)
    if proc_events:
        bed = pd.read_csv(bed_file, sep='\t', header=None)
    else:
        bed = pd.read_csv(bed_file, sep=',')
    
    outputs = {} if proc_events else {x: {} for x in bed['gene_id'].unique()} # only pre fill gene id if enhancers

    for chrom, start, end, name, gene, _ in bed.itertuples(index=False, name=None): # gene colum only makes sense for enhancers
        values = bw.values(chrom, start, end)
        if proc_events:
            values = resize_seq(values, seq_len) 
            outputs[name] = values
        else:
            values = resize_seq_center(values, seq_len)
            outputs[gene][name] = values
    bw.close()
    return outputs


def combine_events(event_parts):
    """
    Combine the different event parts into a np.array of shape (n_events, 3, n_bp, 1)

    :event_parts: A dictionary containing the sequence info per event
    """
    events = sorted([x.split("|")[1] for x in event_parts.keys()]) # events in event_encoding.h5 are sorted
    n_events = len(events)

    histone_arr = np.zeros((n_events, 3, 1024, 1))

    for i, ev in enumerate(events):
        intron1_data = event_parts[f"intron1|{ev}"]
        exon_data = event_parts[f"exon|{ev}"]
        intron2_data = event_parts[f"intron2|{ev}"]

        histone_arr[i, 0, :, 0] = intron1_data
        histone_arr[i, 1, :, 0] = exon_data
        histone_arr[i, 2, :, 0] = intron2_data

    return histone_arr


def combine_enhancers(enhancer_parts, existing_h5):
    # I get a dictionary of enhancer parts i.e., {'gene1': {enhancer: [2000], 'gene2': {enhancer: [2000]}}}
    # and I need a np array per gene with the final shape of (n_genes, 61, 2000, 1) per histone mark
    gene_order = [x.decode() for x in existing_h5['gene_id']]
    full_arr = np.zeros((len(gene_order), 61, 2000, 1))
    for idx, gene in enumerate(gene_order):
        enhancers = enhancer_parts.get(gene, {})
        gene_array = np.zeros((61, 2000, 1))
        for _, values in enhancers.items():
            gene_array[1:len(values), :, 0] = values # start at idx to be consistent with orignal EPInformer
        full_arr[idx] = gene_array
    return full_arr


def create_h5_dataset(h5_file, histone_arr, proc_events=True):
    # per histone mark i need an np.array of shape  (n_events, 3, n_bp, 1)
    print(f"Creating HDF5 dataset at {h5_file}...")
    with h5.File(h5_file, "w") as f:
        # per cell line create a group, in that group each histone mark gets its np array
        for cell, histone_data in histone_arr.items():
            if not proc_events and len(histone_data.keys()) > 1:
                for histone_mark, data in histone_data.items():
                    f.create_dataset(histone_mark, data=data, compression="gzip", compression_opts=5, chunks=(1, 61, 2000, 1))
                continue
            group = f.create_group(cell)
            for histone_mark, data in histone_data.items():
                group.create_dataset(histone_mark, data=data, compression="gzip", compression_opts=5, chunks=(1, 3, 1024, 1))


if __name__ == "__main__":
    process_events = False # if False, process enhancers
    if process_events:
        bed_file = "../data/extracted_events.bed"
        h5_file = "../data/event_bp_histone_marks.h5"
        dataset_raw = process_histone_signals(bed_file, proc_events=process_events)
        create_h5_dataset(h5_file, dataset_raw)
    else:
        # K562
        k562_link_file = "../data/K562_gene_enhancer_links.csv"
        new_k562_h5 = "../data/K562_histone_enhancer_marks.h5"
        existing_k562_h5 = "../data/K562_histone_appended_pe_encoding.h5"
        existing_k562_h5 = h5.File(existing_k562_h5, "r")
        dataset_raw = process_histone_signals(k562_link_file, existing_k562_h5, proc_events=process_events, only_cell='K562')
        create_h5_dataset(new_k562_h5, dataset_raw, proc_events=process_events)
        existing_k562_h5.close()

        # GM12878 repeat as above
        gm12878_link_file = "../data/GM12878_gene_enhancer_links.csv"
        new_gm12878_h5 = "../data/GM12878_histone_enhancer_marks.h5"
        existing_gm12878_h5 = "../data/GM12878_histone_appended_pe_encoding.h5"
        existing_gm12878_h5 = h5.File(existing_gm12878_h5, "r")
        dataset_raw = process_histone_signals(gm12878_link_file, existing_gm12878_h5, proc_events=process_events, only_cell='GM12878')
        create_h5_dataset(new_gm12878_h5, dataset_raw, proc_events=process_events)
        existing_gm12878_h5.close()