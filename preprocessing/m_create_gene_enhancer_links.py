import pandas as pd
import tqdm
import os


def filter_enhancers_by_distance(enhancer_distance, max_distance):
    return enhancer_distance < max_distance


def filter_enhancers_by_promoter_region(enhancer_distance, region_len=2000):
    # enhancers around 2kb of the tss are considered promoters and thus excluded
    return enhancer_distance > region_len

# NOT USED
def filter_enhancers_by_promoter(enhancer_name):
    # if the enhancer name starts with "promoter", it is considered a promoter
    return not enhancer_name.startswith('promoter')


def get_gene_info(gene_starts):
    gene_starts = gene_starts[['chr', 'TargetGene', 'TargetGeneTSS']].copy()
    gene_starts = gene_starts.drop_duplicates().sort_values(by=['TargetGene'])
    return gene_starts


def create_links(enhancer_df, gene_tss, max_dist, promoter_region):
    links = [] 
    for _, gene in tqdm.tqdm(gene_tss.iterrows(), total=len(gene_tss)):
        tss = gene['TargetGeneTSS']
        chrom = gene['chr']
        gene_id = gene['TargetGene']  
        chrom_df = enhancer_df[enhancer_df['chrom'] == chrom].copy()
        if len(chrom_df) == 0:
            raise ValueError(f"chromosome {chrom} not found in enhancer_df")
        chrom_df['tss_distance'] = (((chrom_df['start'] + chrom_df['end']) / 2) - tss).abs()
        chrom_df = chrom_df.sort_values('tss_distance', ascending=True)
        n_enhancers = 0
        for _, row in chrom_df.iterrows():
            enhancer_start = row['tss_distance']
            if (filter_enhancers_by_distance(enhancer_start, max_dist) and 
                filter_enhancers_by_promoter_region(enhancer_start, promoter_region)):
                filtered_tuple = (row['chrom'], row['start'], row['end'], row['name'], gene_id, tss)
                links.append(filtered_tuple)
                n_enhancers += 1
            if n_enhancers >= 60 or row['tss_distance'] > 100_000:
                break
    return pd.DataFrame(links, columns=['chrom', 'start', 'end', 'name', 'gene_id', 'tss'])

if __name__ == "__main__":
    cell_lines = ['K562', 'GM12878']
    cell_line = cell_lines[1]  # Change this to 'GM12878' for the other cell line
    max_distance = 100_000
    # TODO: combine the gene starts since this might result in bugs for other bed files
    if cell_line == 'K562':
        enhancer_bed = '../data/K562_DNase_ENCFF257HEE_hic_4DNFITUOMFUQ_1MB_ABC_nominated/DNase_ENCFF257HEE_Neighborhoods/EnhancerList.bed'
        gene_starts = '../data/K562_DNase_ENCFF257HEE_hic_4DNFITUOMFUQ_1MB_ABC_nominated/Gene-enhancer links/EnhancerPredictions.txt'
    else:
        enhancer_bed = '../data/GM12878_DNase_ENCFF020WZB_hic_4DNFI1UEG1HD_1MB_ABC_nominated/DNase_ENCFF020WZB_Neighborhoods/EnhancerList.bed'
        gene_starts = '../data/GM12878_DNase_ENCFF020WZB_hic_4DNFI1UEG1HD_1MB_ABC_nominated/Gene-enhancer links/EnhancerPredictions.txt'
    
    gene_starts_df = pd.read_csv(gene_starts, sep='\t')
    enhancer_df = pd.read_csv(enhancer_bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'name']) #chr1	11623	12123	intergenic|chr1:11623-12123

    gene_tss = get_gene_info(gene_starts_df) # ['chr', 'TargetGene', 'TargetGeneTSS']
    link_df = create_links(enhancer_df, gene_tss, 100_000, 1000)
    
    link_df.to_csv(f'../data/{cell_line}_gene_enhancer_links.csv', index=False)
