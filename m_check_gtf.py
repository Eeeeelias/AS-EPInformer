import pandas as pd
import numpy as np
import os


isoform_file_k562 = '/nfs/data/IHEC/RNAseq/RNA-Seq/' \
    'ihec.rna-seq.ihec-grapenf-containerv1.1.0.IHECRE00001887.4.b2b4ded3-20fb-4626-b1c9-2afcd18bd6ec.isoforms.results'

isoform_file_gm12878 = '/nfs/data/IHEC/RNAseq/RNA-Seq/' \
    'ihec.rna-seq.ihec-grapenf-containerv1.1.0.IHECRE00001892.5.c5893651-dcc4-47b1-a647-a853e04e006d.isoforms.results'

gtf_file = '/nfs/home/students/e.albrecht/SUPPA/gencode.v29.TSL12.gtf'


# isoform information
isoform_data_k562 = pd.read_csv(isoform_file_k562, sep='\t')
isoform_data_gm12878 = pd.read_csv(isoform_file_gm12878, sep='\t')


### THIS IS JUST TO CHECK THE VERSION OF THE GTF FILE.
# For this project: Ensembl release 95 fits very well.

# split transcript_id version into separate column and sort by that column

isoform_data_k562['transcript_version'] = isoform_data_k562['transcript_id'].str.split('.').str[1]
isoform_data_gm12878['transcript_version'] = isoform_data_gm12878['transcript_id'].str.split('.').str[1]

# sort by transcript_id and transcript_version
# isoform_data = isoform_data.sort_values(by=['transcript_version'], ascending=False)
# print(isoform_data[['transcript_id', 'transcript_version']].head(n=10))

# retrieve transcript id and version from gtf file

transcripts = set()
with open(gtf_file, 'r') as f:
    for line in f.readlines(): 
        if line.startswith("#"):
            continue
        if not "transcript_id" in line:
            continue
        line = line.split(";")
        transcript_name = ""
        for item in line:
            if item.startswith(" transcript_id"):
                transcript_name = item.split("\"")[1]
        if transcript_name != "":
            transcripts.add(transcript_name)

# print the first 10 transcript ids

print("Number of transcripts in gtf file: ", len(transcripts))
print("Number of transcripts in isoform file: ", len(isoform_data_k562['transcript_id'].unique()))
print(f"Matching transcripts: {len(transcripts.intersection(set(isoform_data_k562['transcript_id'].unique())))}")
print(f"Missing transcripts: {len(transcripts.difference(set(isoform_data_k562['transcript_id'].unique())))}")

# check the tpm percentile of missing transcripts
if len(transcripts.difference(set(isoform_data_k562['transcript_id'].unique()))) > 0:
    quant_transcripts = set(isoform_data_k562['transcript_id'].unique())
    gtf_transcripts = set(transcripts) 
    mismatching_transcripts = quant_transcripts.difference(gtf_transcripts)

    missing_tpm = isoform_data_k562[isoform_data_k562['transcript_id'].isin(mismatching_transcripts)]['TPM']
    tpm_threshold = isoform_data_k562['TPM'].quantile(0.95)

    high_expr_missing_info = isoform_data_k562[
        isoform_data_k562['transcript_id'].isin(mismatching_transcripts) & 
        (isoform_data_k562['TPM'] > tpm_threshold)
    ]

    print(high_expr_missing_info)

### write transcripts + tpm to file

reduced_isoform = isoform_data_k562[['transcript_id', 'TPM']]
reduced_isoform = reduced_isoform.rename(columns={'TPM': 'K562_TPM'})
reduced_isoform = pd.merge(
    reduced_isoform,
    isoform_data_gm12878[['transcript_id', 'TPM']],
    on='transcript_id',
    how='left'
)
reduced_isoform = reduced_isoform.rename(columns={'TPM': 'GM12878_TPM'})
print(reduced_isoform.columns)

# write to file but include header line with 'sample1'
with open(f'data/K562_GM12878_transcript_tpm.txt', 'w') as f:
    f.write('K562_TPM\tGM12878_TPM\n')
    reduced_isoform.to_csv(f, sep='\t', index=False, header=False, mode='a')