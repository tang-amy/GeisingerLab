## Yunfei Dai
## 2021/10/21

# This script takes blastp result table as input, and outputs a fasta file containing the protein sequences of all hits.

import pandas as pd
from Bio import SeqIO
from optparse import OptionParser

options = OptionParser()
options.add_option("-i", "--infile", dest="infile", help="Please provide input table.")
options.add_option("-r", "--reference", dest="reference", help="Please provide reference fasta file.")
options.add_option("-o", "--outfile", dest="outfile", default="test_weboflife.fasta")

def get_seq_dict(reference):
    fasta_sequences = SeqIO.parse(open(reference), 'fasta')
    seq_dict = {rec.id : rec.seq for rec in fasta_sequences}
    return seq_dict

def main():
    opts, args = options.parse_args()
    infile = opts.infile
    reference = opts.reference
    outfile = opts.outfile

    seq_dict = get_seq_dict(reference)

    fasta_output = open(outfile, 'a+')
    with open(infile, 'r') as blast_result:
        for line in blast_result:
            if "#" not in line:
                line = line.strip()
                if len(line.split('\t')) > 1:
                    hit = line.split('\t')[1]
                else:
                    hit = line
                header = ">" + hit
                seq = str(seq_dict[hit])
                fasta_output.write(header + '\n')
                fasta_output.write(seq + '\n')
    fasta_output.close()
if __name__ == '__main__':
        main()

