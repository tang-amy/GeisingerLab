## Yunfei Dai
## 03/26/2020

"""

This script identifies all n-base (default is 20) protospacer sequences in the genome based on protospacer adjacent motifs (PAM),
and outputs their positions.

Packages required: Bio, pandas, optparse, re

Usage: python3 PAM_finder.py -i [input_directory] -n [PAM_length] -o [output_directory]

-i input file in fasta format, multiple-chromosome is supported

"""

# !/home/bin/python

from optparse import OptionParser
from Bio import SeqIO
import re
import pandas as pd

options = OptionParser()
options.add_option("-i", "--infile", dest="infile",
                   help="provide input directory")
options.add_option("-n", "--PAM_length", dest="PAM_length", default="20",
                   help="desired sequence length, default is 20")
options.add_option("-o", "--output_name", dest="outfile", default="pam_list",
                   help="custom output filename, default is pam_list")


def find_pam(sequence, length=20, outfile="pam_list"):
    length = int(length)
    plus_strand = "+"
    minus_strand = "-"
    with open(sequence, "r") as f:
        for record in SeqIO.parse(f, "fasta"):                  # read each sequence (chromosome)
            seq_ID = record.id
            POS_PAM_plus = {}
            POS_PAM_minus = {}
            PAM_list_plus = []
            PAM_list_minus = []
            counter1 = 0
            counter2 = 0
            for PAM_NGG in re.finditer("GG", str(record.seq)):  # search PAM in plus strand ending with "...NGG"
                if PAM_NGG.start() >= 22:
                    GG_pos = PAM_NGG.start() + 1                # position of GG
                    pam_start_plus = GG_pos - length - 2        # start position of the PAM sequence
                    pam_end_plus = GG_pos - 2                   # end position of the PAM sequence
                    POS_PAM_plus.update({pam_start_plus: [seq_ID, pam_end_plus, plus_strand]})
                    counter1 += 1
            for key in POS_PAM_plus:
                PAM_list_plus.append([POS_PAM_plus[key][0], key, POS_PAM_plus[key][1], POS_PAM_plus[key][2]])
            df = pd.DataFrame(PAM_list_plus, columns=['ID', 'start_pos', 'end_pos', 'strand'])
            print("Found " + str(counter1) + " PAM sequences on \"+\" strand")


            for PAM_CCN in re.finditer("CC", str(record.seq)):  # search PAM in minus strand following "CCN..."
                CC_pos = PAM_CCN.start() + 1
                pam_start_minus = CC_pos + 2
                pam_end_minus = CC_pos + length + 2
                POS_PAM_minus.update({pam_start_minus: [seq_ID, pam_end_minus, minus_strand]})
                counter2 += 1
            for key in POS_PAM_minus:
                PAM_list_minus.append([POS_PAM_minus[key][0], key, POS_PAM_minus[key][1], POS_PAM_minus[key][2]])
            print("Found " + str(counter2) + " PAM sequences on \"-\" strand")
            print("Found " + str(counter1+counter2) + " PAM sequences in total from both strands")
            counter1 = 0
            counter2 = 0
            to_append = pd.DataFrame(PAM_list_minus, columns=df.columns)
            df = df.append(to_append, ignore_index=True)
            df = df.sort_values(by=['start_pos'])
            df.to_csv(outfile, sep='\t', index=False, header=None)

def main():
    opts, args = options.parse_args()
    infile = opts.infile
    pam_length = opts.PAM_length
    outfile = opts.outfile
    find_pam(infile, pam_length, outfile)


if __name__ == '__main__':
    main()
