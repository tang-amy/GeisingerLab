#! /usr/bin/python3

## This script searches for ORFs within 1000 bp of peak position
## Input is bed file with averaged peak sumit positions (average_summit.bed)

import os
import pandas as pd
from sys import argv
from Bio import SeqIO
from scipy.stats import norm
from bisect import bisect_left

def parse_gb(gb):
    # get all sequence records for the specified genbank file
    recs = [rec for rec in SeqIO.parse(gb, "genbank")]
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

def make_orf_list(annotations, chrom):
    gb = os.path.join(annotations, chrom + '.gb')
    dict_chrom = parse_gb(gb)
    start_codon_list = sorted(list(dict_chrom.keys()))
    return dict_chrom, start_codon_list

def is_intra(summit, dict, chrom):
    # return whether a peak summit is within any CDS region
    match = False # intergenic
    for key, value in dict[chrom].items():
        start = value[1]
        end = value[2]
        chrom = value[3]
        if chrom == 1:
            if start <= summit <= end: # genes can overlap, thus there might be more than 1 hit
                match = True # coding
                break
        elif chrom == -1:
            if end <= summit <= start: # on the negative strand gene positions are [end, start]
                match = True # coding
                break
    return match

class NearestORF:
    def __init__(self, summit, gene_dict, start_codon_dict, chrom):
        self.summit = summit
        self.chrom = chrom
        self.dict_annotation = gene_dict[chrom]
        self.start_codon_list = start_codon_dict[chrom]

    def find_nearest_start_codon(self):
        # find the start codon closest to the peak
        match_index = bisect_left(self.start_codon_list, self.summit)
        if match_index == 0:
            nearest_match = self.start_codon_list[0]
        elif match_index == len(self.start_codon_list):
            nearest_match = self.start_codon_list[-1]
        else:
            left = self.start_codon_list[match_index - 1]
            right = self.start_codon_list[match_index]
            if right - self.summit < self.summit - left:
                nearest_match = right
            else:
                nearest_match = left
        return match_index, nearest_match

    def get_match_info(self, match_ORF, strand):
        # for any matched ORF, get its corresponding gene annotation info (accession, gene coordinates, strand, distance to match)
        accession, start, end, ORF_strand = self.dict_annotation[match_ORF]
        if strand == 1:
            distance_to_match = self.summit - match_ORF
        elif strand == -1:
            distance_to_match = match_ORF - self.summit
        match_info = [self.summit, 0, accession, self.chrom, start, end, ORF_strand, distance_to_match, 0] # last column: intergenic distance to last match
        return match_info

    def find_next_ORFs(self, ORF_list_to_search, strand, nearest_match_info):
        # for a given match, search for the next nearest ORFs and output in a list
        # the search will be in one direction based on the strand info
        counter = 0
        all_ORF_match = [nearest_match_info]
        for next_match_index in ORF_list_to_search:
            counter += 1
            next_match_strand = self.dict_annotation[next_match_index][3]
            if next_match_strand == strand:
                match_info = self.get_match_info(next_match_index, next_match_strand)
                match_info[1] = counter # Update Nth nearest ORF
                next_match_distance = match_info[7]
                if abs(next_match_distance) <= 1000:
                    last_ORF_pos = all_ORF_match[-1][4] # start position of last ORF match
                    intergenic_distance = abs(match_info[4] - last_ORF_pos) # calculate intergenic distance
                    match_info[8] = intergenic_distance # update intergenic distance
                    all_ORF_match.append(match_info)
                else:
                    break
            else:
                break
        return all_ORF_match

    def nearest_orf(self):
        nearest_match_index, nearest_match = self.find_nearest_start_codon()
        nearest_match_info = self.get_match_info(nearest_match, self.dict_annotation[nearest_match][3])
        nearest_match_strand = self.dict_annotation[nearest_match][3]
        if  nearest_match_strand == 1:  # If nearest ORF is on + strand, search right
            ORF_list = self.start_codon_list[(nearest_match_index + 1):]
        elif nearest_match_strand == -1:  # If nearest ORF is on - strand, search left
            if self.summit != nearest_match:
                ORF_list = self.start_codon_list[:nearest_match_index-1][::-1]  # List is reversed
            else: # if summit position is the same as first ORF, include this position in ORF_list
                ORF_list = self.start_codon_list[:nearest_match_index][::-1]
        all_ORF_match = self.find_next_ORFs(ORF_list, self.dict_annotation[nearest_match][3], nearest_match_info)    
        return all_ORF_match
    
def make_match_table(infile, outfile, gene_dict, start_codon_dict):
    # write ORF matches to a tsv file
    df_peak = pd.read_csv(infile, sep='\t')
    distance_list_intergenic = []
    distance_list_coding = []
    match_stats = []
    for index, peak in df_peak.iterrows():
        chrom = peak['chrom']
        summit = peak['average_summit']
        fold_enrich = peak['average_enrichment']
        orf_finder = NearestORF(summit, gene_dict, start_codon_dict, chrom)
        all_ORFs = orf_finder.nearest_orf()
        for match_info in all_ORFs:
            distance = match_info[-1]
            if is_intra(summit, gene_dict, chrom):
                match_info.append('coding')
                distance_list_coding.append(distance)
            else:
                match_info.append('intergenic')
                distance_list_intergenic.append(distance)
            match_info.append(fold_enrich)
            match_stats.append(match_info)
    
    # write output to tsv
    # add column 'Nth nearest ORF'
    col_names = ['summit_pos', 'Nth nearest ORF', 'locus_tag', 'chrom', 'start', 'end', 'strand', 'distance_to_match', 'intergenic_distance', 'match_type', 'average_fold_enrichment']
    df_out = pd.DataFrame(match_stats, columns=col_names)
    df_out = df_out.sort_values(by = ['summit_pos', 'Nth nearest ORF'], ascending = [True, True])
    df_out.to_csv(outfile, sep='\t', index=False)
    
def main():
    infile = argv[1]
    outfile = argv[2]
    annotations = argv[3] # directory containing gb files 

    gene_dict = {}
    start_codon_dict = {}
    for fname, chrom in {'NZ_CP012004.1':'NZ_CP012004.1', 
                        'pAb1' : 'NC_009083.1', 
                        'pAb2' : 'NC_009084.1', 
                        'pAb3' : 'NZ_CP012005.1'}.items():
        dict_chrom, start_codon_list = make_orf_list(annotations, fname)
        gene_dict.update({chrom: dict_chrom})
        start_codon_dict.update({chrom : start_codon_list})

    make_match_table(infile, outfile, gene_dict, start_codon_dict)

if __name__ == "__main__":
    main()