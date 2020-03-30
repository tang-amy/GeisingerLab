# Identification of sgRNA sequences from genome

## Requirements
+ Python3 with the following packages installed:
  + biopython (Bio)
  + pandas 
+ Reference files:
  + Genomic FASTA file (eg. NZ_CP012004.fasta)
  + List of transcription start sites (TSSs) (Kroger_TSS.txt, from Kroger _et al_, 2017)
 
## Steps
1. Identify all pam sequences from the genome

2. Identify sgRNA sequences from the list of PAM sequences in step 1
