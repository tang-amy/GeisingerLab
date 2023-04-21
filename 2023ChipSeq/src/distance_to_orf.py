## Yunfei Dai
## 16Feb2022

import os, re, sys
import pandas as pd, numpy as np
from matplotlib import pyplot as plt
#import matplotlib.mlab as mlab
from Bio import SeqIO
from scipy.stats import norm
from bisect import bisect_left

def parse_gb(gb):
    # get all sequence records for the specified genbank file
    recs = [rec for rec in SeqIO.parse(gb, "genbank")]
    # print the number of sequence records that were extracted
    # print(len(recs))
    # print annotations for each sequence record
    for rec in recs:
        annotation_type = rec.annotations
    # print the gene sequence feature summary information for each feature in each
    # sequence record
    dict = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "gene"]
        for feat in feats:
            locus_tag = feat.qualifiers['locus_tag']
            locus_tag = ";".join(locus_tag)
            location = feat.location
            strand = int(location.strand)
            # flip start and end if gene is on negative strand
            if strand == 1:
                start_pos = int(location.start)
                end_pos = int(location.end)
            elif strand == -1:
                start_pos = int(location.end)
                end_pos = int(location.start)
            dict.update({start_pos: (locus_tag, start_pos, end_pos, strand)})
    return dict

def is_intra(summit, dict, chrom):
    # return whether a peak summit is within any CDS region
    match = False
    for key, value in dict[chrom].items():
        start = value[1]
        end = value[2]
        if start <= summit <= end: # genes can overlap, thus there might be more than 1 hit
            match = True
            break
    return match

def nearest_orf(summit, gene_dict, start_codon_dict, chrom):
    # for a given peak summit, return its nearest start codon and gene info
    dict_annotation = gene_dict[chrom]
    start_codon_list = start_codon_dict[chrom]
    match_index = bisect_left(start_codon_list, summit)
    if match_index == 0:
        nearest_match = start_codon_list[0]
    elif match_index == len(start_codon_list):
        nearest_match = start_codon_list[-1]
    else:
        left = start_codon_list[match_index - 1]
        right = start_codon_list[match_index]
        if right - summit < summit - left:
            nearest_match = right
        else:
            nearest_match = left
    accession = dict_annotation[nearest_match][0]
    start = dict_annotation[nearest_match][1]
    end = dict_annotation[nearest_match][2]
    strand = dict_annotation[nearest_match][3]
    # 21Apr2023 changed the calculation of distance so that:
    # distance is positive: for a peak left to a start codon on the + strand
    # distance is negative: for a peak right to a start codon on the + strand
    # distance is positive: for a peak right to a start codon on the - starnd
    # distance is negative: for a peak left to a start codon on the - strand
    # distance means "distance from start codon"
    if strand ==  1:
        distance_to_match = nearest_match - summit
    elif strand == -1:
        distance_to_match = summit - nearest_match
    match_info = [accession, start, end, strand, summit, distance_to_match]
    return match_info 

def make_histogram(infile, outfile, distance_plot, gene_dict, start_codon_dict):
    # make histogram, and output stats to a tsv file
    df_peak = pd.read_csv(infile, sep='\t')
    distance_list_intergenic = []
    distance_list_coding = []
    match_stats = []
    for index, peak in df_peak.iterrows():
        chrom = peak['chrom']
        summit = peak['average_summit']
        match_info = nearest_orf(summit, gene_dict, start_codon_dict, chrom)
        distance = match_info[-1]
        fold_enrich = peak['average_enrichment']
        if is_intra(summit, gene_dict, chrom):
            match_info.append('coding')
            distance_list_coding.append(distance)
        else:
            match_info.append('intergenic')
            distance_list_intergenic.append(distance)
        match_info.append(fold_enrich)
        match_stats.append(match_info)
    
    # write output to csv
    col_names = ['locus_tag', 'start', 'end', 'strand', 'summit_pos', 'distance_to_match', 'match_type', 'average_fold_enrichment']
    df_out = pd.DataFrame(match_stats, columns=col_names)
    df_out.to_csv(outfile, sep='\t', index=True)
    
    # 23Mar2023 changed default bin size to 20
    bins = np.arange(np.min(distance_list_coding), np.max(distance_list_coding), 20)
    fig, ax = plt.subplots()
    ax.hist(distance_list_coding, bins=bins, alpha=0.5, label='coding')
    ax.hist(distance_list_intergenic, bins=bins, alpha=0.5, label='intergenic')
    ax.set_xlabel('Distance from nearest start codon (nt)')
    ax.set_ylabel('Count')
    ax.set_xlim(-500,500) # 23Mar2023 changed x axis limit to (-500, 500)
    ax.legend()
    ax.set_title(os.path.basename(infile))
    ax.spines[['right','top']].set_visible(False)
    fig.tight_layout()
    plt.savefig(distance_plot)
'''
    # add a 'best fit' line
    mu1, sigma1 = norm.fit(distance_list_coding)
    y1 = mlab.normpdf(bins, mu1, sigma1)
    ax.plot(bins, y1, 'b--', linewidth=1)
    
    mu2, sigma2 = norm.fit(distance_list_intergenic)
    y2 = mlab.normpdf(bins, mu2, sigma2)
    ax.plot(bins, y2, 'r--', linewidth=1)
'''

def main():
    # input is .intersect.bed (from multiinter)
    infile = sys.argv[1]
    outfile = sys.argv[2]
    distance_plot = sys.argv[3]

    annotations = "/work/geisingerlab/Yunfei/2022_Nov_ChipSeq_Analysis/References/annotations"
    genome = os.path.join(annotations, "NZ_CP012004.1.gb")
    pAb1 = os.path.join(annotations, "pAb1.gb")
    pAb2 = os.path.join(annotations, "pAb2.gb")
    pAb3 = os.path.join(annotations, "pAb3.gb")

    dict_genome = parse_gb(genome)
    dict_pAb1 = parse_gb(pAb1)
    dict_pAb2 = parse_gb(pAb2)
    dict_pAb3 = parse_gb(pAb3)
    gene_dict = {'NZ_CP012004.1': dict_genome, 'NC_009083.1': dict_pAb1, 'NC_009084.1': dict_pAb2, 'NZ_CP012005.1': dict_pAb3}

    start_codon_dict = {}
    for key in gene_dict:
        start_codon_list = sorted(list(gene_dict[key].keys()))
        start_codon_dict.update({key: start_codon_list})

    make_histogram(infile, outfile, distance_plot, gene_dict, start_codon_dict)


if __name__ == "__main__":
    main()
