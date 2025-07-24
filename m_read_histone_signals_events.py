import glob
import subprocess
import pandas as pd

file_uuids = {
    'H3K27ac': {
        'K562': '92e33eb4-60fb-4aba-99b3-db60d30e791f',
        'GM12878': '1a48d1e3-a2d4-434f-8cee-d284ec0da60a'
    },
    'H3K27me3': {
        'K562': '59a9b7d0-5484-46d1-b8f0-117dbd401378',
        'GM12878': '902ae139-ccdf-4a1c-a16a-ec321c33c554'
    },
    'H3K36me3': {
        'K562': '8d92c560-11be-4c39-83b5-4989c5ac9c68',
        'GM12878': 'f92b9e80-4bf7-4519-837e-1a5c78d2d2b2'
    },
    'H3K4me1': {
        'K562': 'e03cd402-356e-42a8-a5d3-47467ef89fb6',
        'GM12878': 'd37ba2fd-c103-41c2-a854-82a474c3bdc9'
    },
    'H3K4me3': {
        'K562': 'aa825751-1900-4a19-81df-fe3780ad3aed',
        'GM12878': '68f6c400-c106-443e-8636-56b5afbcfe6d'
    },
    'H3K9me3': {
        'K562': 'a6734a5b-3b03-45c9-8bf9-b4b30b035291',
        'GM12878': 'a3994b5d-df9e-4d2b-8d1f-2a68f9eea142'
    }
}

bigwig_base_dir = '/nfs/data3/IHEC/ChIP-Seq/'
output_dir = '/nfs/proj/EPInformer-AS/AS-Epinformer-Dup/data/IHEC-ChIP-Seq-Histone-Signals'

def get_events(file_path):
    event_parts = {'event_id': [], 'intron1': [], 'intron2': [], 'exon': [], 'chr': [], 'strand': []}
    events = pd.read_csv(file_path)
    longest_exon = 0
    for _, row in events.iterrows():
        event = row['event_id'] # ENSG00000224051;SE:chr1:1325102-1326836:1327032-1327241:+
        event_parts['event_id'].append(event)
        parts = event.split(':')
        event_parts['chr'].append(parts[1]) # chromosome
        event_parts['exon'].append(f"{parts[2].split('-')[1]}-{parts[3].split('-')[0]}") # exon start
        event_parts['strand'].append(parts[4]) # strand

        # We take the last 1024 bases from intron1 and the first 1024 bases from intron2 if they're longer than that
        intron1_len = int(parts[2].split('-')[1]) - int(parts[2].split('-')[0])
        intron2_len = int(parts[3].split('-')[1]) - int(parts[3].split('-')[0])
        if not intron1_len > 1024:
            event_parts['intron1'].append(parts[2]) # first intron
        else:
            event_parts['intron1'].append(f"{int(parts[2].split('-')[1]) - 1024}-{int(parts[2].split('-')[1])}")
        if not intron2_len > 1024:
            event_parts['intron2'].append(parts[3])
        else:
            event_parts['intron2'].append(f"{int(parts[3].split('-')[0])}-{int(parts[3].split('-')[0]) + 1024}")

        # exon_len = int(parts[3].split('-')[0]) - int(parts[2].split('-')[1])
    return pd.DataFrame(event_parts)


def create_bed_file(events, output_file):
    with open(output_file, 'w') as f:
        for idx, row in events.iterrows():  # chr1	11623	12123	intergenic|chr1:11623-12123
            intron1_start, intron1_end = row['intron1'].split("-") 
            f.write(f"{row['chr']}\t{intron1_start}\t{intron1_end}\tintron1|{row['event_id']}\t0\t{row['strand']}\n")

            exon_start, exon_end = row['exon'].split("-")
            f.write(f"{row['chr']}\t{exon_start}\t{exon_end}\texon|{row['event_id']}\t0\t{row['strand']}\n")

            intron2_start, intron2_end = row['intron2'].split("-") 
            f.write(f"{row['chr']}\t{intron2_start}\t{intron2_end}\tintron2|{row['event_id']}\t0\t{row['strand']}\n")


    print(f"BED file created at {output_file}")


def extract_histone_signals(bigwig_file, bed_file, output_file):
    cmd = [
        'bigWigAverageOverBed',
        bigwig_file,
        bed_file,
        output_file
    ]
    subprocess.run(cmd, check=True)


def get_file_name(uuid):
    # glob bigwig base dir for uuid and return the first match
    pattern = f"{bigwig_base_dir}*{uuid}*fc.signal.bw"
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    else:
        raise FileNotFoundError(f"No file found for UUID: {uuid}")
    

def process_histone_signals():
    output_files = []
    for histone_type, tissues in file_uuids.items():
        for tissue, uuid in tissues.items():
            bigwig_file = get_file_name(uuid)
            bed_file = f"data/extracted_events.bed"
            name_suffix = "Events"
            output_file = f"{output_dir}/{name_suffix}{tissue}.{histone_type}.tab"

            print(f"Processing {histone_type} for {tissue}...")
            extract_histone_signals(bigwig_file, bed_file, output_file)
            print(f"Output saved to {output_file}")
            output_files.append(output_file)
    return output_files


def combine_signal_files(output_files, column='mean'):
    output_dfs = {
        'K562': [],
        'GM12878': []
    }
    for file in output_files:
        tissue = 'K562' if 'K562' in file else 'GM12878'
        df = pd.read_csv(file, sep='\t', names=['gene_id', 'size', 'covered', 'sum', 'mean0', 'mean'])
        df = df[['gene_id', 'covered', column]]
        df.rename(columns={column: f"{file.split('.')[-2]}.{column}"}, inplace=True)
        output_dfs[tissue].append(df)

    name_suffix = "Events"
    for tissue, dfs in output_dfs.items():
        combined_df = pd.concat(dfs, axis=1)
        combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]
        combined_df.to_csv(f"{output_dir}/{name_suffix}_Combined_{tissue}_Histone_Signals.csv", index=False)
        print(f"Combined data for {tissue} saved to {output_dir}/{name_suffix}_Combined_{tissue}_Histone_Signals.csv")


if __name__ == "__main__":
    events = get_events("data/psi_response.csv")
    create_bed_file(events, 'data/extracted_events.bed')
    output_files = process_histone_signals()
    combine_signal_files(output_files)
    print("All histone signals for events extracted successfully.")