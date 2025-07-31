import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from Bio.Seq import Seq
import h5py
import random
import tqdm

# RUN THIS SCRIPT AFTER m_create_gene_enhancer_links.py

# first get the coordinates of the enhancers
enhancer_links = pd.read_csv('../data/GM12878_gene_enhancer_links.csv')

# get the fasta
fasta_file = pysam.FastaFile('/nfs/data/references/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')



# create a tuple list of enhancer chr, start, end, name 
enhancer_tuples = set(zip(enhancer_links['chrom'], enhancer_links['start'], enhancer_links['end'], enhancer_links['name']))

# create a dictionary to store the sequences
enhancer_sequences = {}
# iterate over the enhancer tuples and get the sequences
for chrom, start, end, name in tqdm.tqdm(enhancer_tuples):
    seq = fasta_file.fetch(chrom, start, end)
    enhancer_sequences[(chrom, start, end, name)] = seq

# Add the sequences to the enhancer_links DataFrame
enhancer_links['sequence'] = enhancer_links.apply(lambda row: enhancer_sequences.get((row['chrom'], row['start'], row['end'], row['name'])), axis=1)

# Save the updated DataFrame to a new CSV file
enhancer_links.to_csv('../data/GM12878_gene_enhancer_links_with_sequences.csv', index=False)