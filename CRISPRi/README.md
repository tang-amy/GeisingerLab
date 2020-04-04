# Identification of sgRNA sequences from genome

## Requirements
+ Python3 with the following packages installed:
  + biopython (Bio)
  + pandas 
+ Reference files:
  + Genomic FASTA file (eg. NZ_CP012004.fasta)
  + List of transcription start sites (TSSs) (Kroger_TSS.txt, from Kroger _et al_, 2017)
 
## Steps
**1. Find PAM regions from genome**
  Use `PAM_finder.py` to identify all pam regions in the genome. (PAM regions could be preceeding NGG if on the + strand, or following CCN on the - strand. )
    
   - Usage:
   ```
    python3 PAM_finder.py -i [input.fasta] -o [output.bed]
   ```
    
   - Options
   ```
    -l desired short guide length, default is 20 bp.
   ```
   
   - Output:
   ```
    Column 1: ID
    Column 2: Type (NGG) or (CCN)
    Column 3: Start position (SGR)
    Column 4: End position (SGR)
    Column 5: PAM coordinate (position of "N")
   ```
    
  
**2. Find sgRNA sequences**
  Use `sgRNA_finder.py` to generate a list of sgRNA list that could be used for CRISPRi knock-down.
  
   - Usage:
   ```
   python3 sgRNA_finder.py -i [input.bed] -t [TSS list.txt] -r [genome.fasta] -g [genome.gbk] -o [output.tsv]
   ```
   
   - Options:
   ```
    -u acceptable upstream range to tss (default is 50 bp)
    -d acceptable downstream range to tss (default is 100 bp)
    -l length of desired sgRNA sequences (default is 20 bp)
   ```
   
   - Description:
   ```
   The two main functions in this script are tss_pam and  sgRNA_finder
   ```
   
   - Output:
   ```
    Column 1: ACX_60 locus tag (old as in Kroger TSS list)
    Column 2: ACX_RS60 locus tag (new)
    Column 3: protein id
    Column 4: TSS strand
    Column 5: TSS coordinate
    Column 6: pam coordinate
    Column 7: pam sequence
    Column 8: SGR start position
    Column 9: SGR end position
    Column 10: target strand (template or non-template)
    Column 11: Distance (tss to the beginning / end of SGR sequence, whichever is longer)
    Column 11: SGR sequence
    Column 12: Seed number (total occurance of the 12 bp region preceeding PAM plus "NGG" in the entire genome)
   ```
