## Yunfei Dai
## 05/06/2022

# This script translates DNA sequences in a multifasta file into protein sequences.
# Translation starts with ATG, if a DNA sequence does not start with ATG but ends with TAC, the reverse complement sequence will be translated.
# Usage:
# python translate_DNAtoProtein.py [INPUT] [OUTPUT]

from sys import argv
from Bio import SeqIO

try:
    infile = argv[1]
except Exception:
    print("Please provide input fasta file.")

try:
    outfile = argv[2]
except Exception:
    outfile = "translated.fasta"

def translate_inframe(infile, outfile):
    fasta_sequences = SeqIO.parse(open(infile), 'fasta')
    with open(outfile, "a+") as protein_fasta:
        for record in fasta_sequences:
            if record.seq[0:3]=="ATG":
                protein_fasta.write(">"+str(record.id) + "\n")
                protein_fasta.write(str(record.seq.translate()) + "\n")
            elif record.seq[-3:]=="CAT":
                rc_record = record.seq.reverse_complement()
                protein_fasta.write(">" + str(record.id) + "\n")
                protein_fasta.write(str(rc_record.translate()) + "\n")
            else:
                print(str(record.id) + ": no start codon found.")

def main():
    translate_inframe(infile, outfile)

if __name__ == '__main__':
    main()
