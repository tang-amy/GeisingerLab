## Yunfei Dai
## 04/24/2020

"""
counts number of bmfR motifs in each fasta file
"""

from optparse import OptionParser
import pandas as pd, re
from Bio import SeqIO, Seq



options = OptionParser()

options.add_option("-i", "--input", dest="infile",
                   help="provide input fasta file")
options.add_option("-m", "--motif", dest="motif",
                   default="/Users/yunfei/GeisingerLab/Chip-Seq_Analysis/BfmR_Motif_Searches.csv",
                   help="provide list of query motifs")
opts, args = options.parse_args()
infile = opts.infile
motif = opts.motif
df_motif = pd.read_csv(motif)
motifs = df_motif["Match"].tolist()

with open(infile, "r") as f:
    counter = 0
    for record in SeqIO.parse(f, "fasta"):  # read sequence(s) from fasta
        seq = str(record.seq)
        for motif in motifs:
            if motif in seq or Seq.reverse_complement(motif) in seq:
                counter += 1
    print("Find", counter, " matches in ", infile)
