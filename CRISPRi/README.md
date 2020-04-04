# Identification of sgRNA sequences from genome

## Requirements
+ Python3 with the following packages installed:
  + biopython (Bio)
  + pandas 
+ Reference files:
  + Genomic FASTA file (eg. NZ_CP012004.fasta)
  + List of transcription start sites (TSSs) (Kroger_TSS.txt, from Kroger _et al_, 2017)
 
## Steps
1. List all pam sequences from fasta genome file
  Use PAM_finder.py to identify all pam regions in the genome. (PAM regions could be preceeding NGG if on the + strand, or following CCN on the - strand. )
  
    * Usage:
  ```
  python3 PAM_finder.py -i [input.fasta] -o [output.bed] -l [shortguide length]
  ```
2. Identify sgRNA sequences from the list of PAM sequences in step 1
