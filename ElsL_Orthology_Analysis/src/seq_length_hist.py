## Yunfei Dai
## 10/06/2021

from sys import argv
from Bio import SeqIO
from matplotlib import pyplot as plt
import pandas as pd

try:
    infile = argv[1]
except Exception:
    print("Please provide input table.")

try:
    outfile = argv[2]
except Exception:
    outfile = "output_histogram.pdf"

try:
    infile_type = argv[3]
except Exception:
    infile_type = "info_table"

def get_length(seq_info, seq_type):
    if seq_type == "fasta":
        # Read fasta file that contains protein sequences
        fasta_sequences = SeqIO.parse(open(seq_info), 'fasta')
        seq_dict = {rec.id : rec.seq for rec in fasta_sequences}
        seq_lengths = []
        for seq_record in seq_dict:
            seq_lengths.append(len(seq_dict[seq_record]))      
    else:
        df_info = pd.read_csv(infile, sep='\t')
        seq_lengths = df_info["Length"].tolist()
    return seq_lengths

def main():
    seq_lengths = get_length(infile, infile_type)
    plt.figure(figsize=(8,6))
    plt.hist(seq_lengths, bins='auto', alpha=0.5, label=infile)
    plt.xlabel("Sequence Length")
    plt.ylabel("Count")
    plt.legend(loc="upper right")
    plt.savefig(outfile, format='pdf')
    
if __name__ == '__main__':
    main()
