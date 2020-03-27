## Yunfei Dai
## 03/26/2020

"""


This script identifies all n-base (default is 20) protospacer sequences in the genome based on protospacer adjacent motifs (PAM),
and outputs their positions.

Packages required: Bio, pandas

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
            POS_PAM = {}
            for PAM_NGG in re.finditer("GG", str(record.seq)):  # search PAM in plus strand ending with "...NGG"
                if PAM_NGG.start() >= 23:
                    GG_pos = PAM_NGG.start() + 1                # position of GG
                    pam_start_plus = GG_pos - length - 2        # start position of the PAM sequence
                    pam_end_plus = GG_pos - 2                   # end position of the PAM sequence
                    POS_PAM.update({pam_start_plus: [seq_ID, pam_end_plus, plus_strand]})
            for PAM_CCN in re.finditer("CC", str(record.seq)):  # search PAM in minus strand following "CCN..."
                CC_pos = PAM_CCN.start() + 1
                pam_start_minus = CC_pos + 3
                pam_end_minus = CC_pos + length + 3
                POS_PAM.update({pam_start_minus: [seq_ID, pam_end_minus, minus_strand]})
            PAM_list = []
            for key in POS_PAM:
                PAM_list.append([POS_PAM[key][0], key, POS_PAM[key][1], POS_PAM[key][2]])
            df = pd.DataFrame(PAM_list, columns=['ID', 'start_pos', 'end_pos', 'strand'])
            df.to_csv(outfile, sep='\t', index=False)

def main():
    opts, args = options.parse_args()
    infile = opts.infile
    pam_length = opts.PAM_length
    outfile = opts.outfile
    find_pam(infile, pam_length, outfile)


if __name__ == '__main__':
    main()
