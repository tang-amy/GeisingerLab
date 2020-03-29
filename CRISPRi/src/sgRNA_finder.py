## Yunfei Dai
## 03/27/2020

"""

This script identifies potential sgRNAs that target regions in proximity to transcription start sites (TSSs).

The default range queried is -50 to +100 of the TSSs (available from Kroger et al. 2017).

The script also checks for specificity of sgRNAs, only unique target sequences are kept.

Input file is .bed generated with PAM_finder.py

The output is a tab delimited txt file with 5 columns.
    Column 1: chromosome ID
    Column 2: "+" or "-" strand
    Column 3: sgRNA start position
    Column 4: sgRNA end position
    Column 5: sgRNA sequence
    Column 6: ACX60_ locus_tag (could be multiple for one sgRNA)

[mandatory options]
    -i input file (.bed)
    -t TSS list (.csv)
    -g reference genome file (.fasta)
    -o output file

[optional]
    -u acceptable upstream range to tss (default is 50)
    -d acceptable downstream range to tss (default is 100)

"""
from optparse import OptionParser
import os
from typing import Dict, Any

import pandas as pd
from Bio import SeqIO
import timeit

options = OptionParser()
options.add_option("-i", "--infile", dest="infile",
                   help="input list of PAM in .bed format)")
options.add_option("-t", "--tss", dest="tss",
                   help="input tss list (tab delimited)")
options.add_option("-g", "--reference", dest="reference",
                   help="reference genome in .fasta format")
options.add_option("-u", "--up_range", dest="up_range", default="50")
options.add_option("-d", "--down_range", dest="down_range", default="100"),
options.add_option("-o", "--output_name", dest="outfile")


# Reads in pam_list.bed and TSS list.
# Identify pam sequences in proximity to at least one TSS.
# Output as "sg_candidate_[infile].bed"

def tss_pam(pam, tss, up_range, down_range, shortlist):
    df_PAM = pd.read_csv(pam, sep='\t', names=["ID", "start_pos", "end_pos", "strand"])
    df_TSS = pd.read_csv(tss, sep='\t')
    pam_index = df_PAM.index.tolist()  # index for all pam sequences
    all_pam = df_PAM["start_pos"].tolist()  # start_pos of pam sequences
    pam_chrom = df_PAM["ID"].tolist()
    pam_end = df_PAM["end_pos"].tolist()
    pam_strand = df_PAM["strand"].tolist()
    tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()[:-1]]

    close_pam_index = []
    for m in pam_index:
        pam_to_search = all_pam[m]
        v = binary_search(tss_coordinate, 0, len(tss_coordinate) - 1, pam_to_search, up_range, down_range)
        if v != -1:
            close_pam_index.append(m)

    short_list = []
    for i in close_pam_index:
        chrom_id = pam_chrom[i]
        pam_start = all_pam[i]
        pam_end_pos = pam_end[i]
        strand = pam_strand[i]
        short_list.append([chrom_id, pam_start, pam_end_pos, strand])
    df = pd.DataFrame(short_list, columns=['ID', 'pam_start_pos', 'pam_end_pos', 'strand'])
    print("Found " + str(len(short_list)) + " PAM sequences in proximity to at least one tss" + "\n")

    df.to_csv(shortlist, sep='\t', index=False, header=False)
    print("PAM shortlist saved as: " + shortlist + "\n")
    return shortlist  # returns path of the shortlist of candidates


def sgRNA_finder(sg_candidate, tss, reference, up_range, down_range, sg_outfile):
    df_candidate = pd.read_csv(sg_candidate, sep='\t', names=["ID", "start_pos", "end_pos", "strand"])
    candidate_start = df_candidate["start_pos"].tolist()  # start_pos of pam sequences
    candidate_chrom = df_candidate["ID"].tolist()
    candidate_end = df_candidate["end_pos"].tolist()
    candidate_strand = df_candidate["strand"].tolist()
    pam_index = df_candidate.index.tolist()
    t2 = timeit.default_timer()

    df_TSS = pd.read_csv(tss, sep='\t')
    tss_index = df_TSS.index.tolist()[:-1]
    tss_locus = df_TSS["Locus_tag"][:-1].tolist()
    tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()[:-1]]

    # screen seed regions to avoid off-target effect (12 bp seed regions of the sgRNAs are screened)
    with open(reference, "r") as f:
        for record in SeqIO.parse(f, "fasta"):  # read sequence(s) from fasta
            gDNA = str(record.seq)
            candidate_seq = {}
            sg_list = []
            for n in pam_index:
                sequence = str(record.seq[candidate_start[n]:candidate_end[n]])  # sgRNA sequence (default 20 bp)
                seed = sequence[0:12]  # seed region of the sgRNA (12bp)
                rc_sequence = reverse_complement(sequence)
                rc_seed = rc_sequence[0:12]
                if candidate_strand[n] == '+':  # if NGG, keep original sequence
                    if gDNA.count(seed) == 1:
                        if gDNA.count(reverse_complement(seed)) == 1:
                            candidate_seq.update({n: []})
                            candidate_seq[n].append(sequence)
                else:  # if CCN, keep reverse complementary
                    if gDNA.count(rc_seed) == 1:
                        if gDNA.count(reverse_complement(rc_seed)) == 1:
                            candidate_seq.update({n: []})
                            candidate_seq[n].append(rc_sequence)
            t3 = timeit.default_timer()
            print("Finished screening seed regions in: ", format((t3 - t2), '.2f'), " seconds." + "\n")

            sgRNA = {}
            for m in candidate_seq:
                for t in tss_index:
                    if -up_range < candidate_start[m] - tss_coordinate[t] < down_range - 20:
                        sgRNA.update({m: candidate_seq[m]})
                        sgRNA[m].append(tss_locus[t])
                        # candidate_seq = {candidate index: [sequence, tss_coordinate] }
            # look for locus_tags and make output file
            for key in sgRNA:
                target = '; '.join(map(str, sgRNA[key][1:]))
                sg_list.append([candidate_chrom[key], candidate_strand[key], candidate_start[key], candidate_end[key],
                                sgRNA[key][0], target])
            df = pd.DataFrame(sg_list, columns=['ID', 'strand', 'SGR_start', 'SGR_end', 'SGR sequence', 'Target(s)'])
            df.to_csv(sg_outfile, sep='\t', index=False)
    t4 = timeit.default_timer()
    print("Finished recording sgRNA candidates: ", format((t4 - t3), '.2f'), " seconds." + "\n")
    print("Recorded ", len(candidate_seq), "sgRNA sequences." + "\n")
    print("sgRNA sequences saved as: ", sg_outfile, ".", "\n")


def binary_search(arr_tss, l, r, p, up, down):
    if r >= l:
        mid = l + (r - l) // 2;
        if -up <= p - arr_tss[mid] <= down - 20:
            return mid
        elif p - arr_tss[mid] < -up:
            return binary_search(arr_tss, l, mid - 1, p, up, down)
        else:
            return binary_search(arr_tss, mid + 1, r, p, up, down)
    else:
        return -1


def reverse_complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq = seq[::-1]
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)


def main():
    opts, args = options.parse_args()
    pam = opts.infile
    pam_path = os.path.abspath(pam)
    pam_dir = os.path.dirname(pam_path)
    pam_filename = os.path.basename(pam_path)
    tss = opts.tss
    reference = opts.reference
    up_range = int(opts.up_range)
    down_range = int(opts.down_range)
    sg_outfile = opts.outfile

    base_pam = os.path.splitext(pam_filename)[0]
    shortlist = pam_dir + '/' + 'PAM_shortlist_' + base_pam + '.bed'

    if not (os.path.isabs(sg_outfile)):
        sg_outfile = pam_dir + '/' + sg_outfile

    t0 = timeit.default_timer()
    sg_candidate = tss_pam(pam, tss, up_range, down_range, shortlist)
    t1 = timeit.default_timer()
    print("finished PAM shortlist in: ", format(t1 - t0, '.2f'), " seconds" + "\n")
    sgRNA_finder(sg_candidate, tss, reference, up_range, down_range, sg_outfile)


if __name__ == '__main__':
    main()
