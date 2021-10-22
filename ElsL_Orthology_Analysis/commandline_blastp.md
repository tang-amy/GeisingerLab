
### 1. Download preformatted blast database from NCBI

On bash command line terminal, use the `update_blastdb.pl` script [(user manual)](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory) in the blast+ package to download the preformatted blast database files available in the [NCBI ftp site](https://ftp.ncbi.nlm.nih.gov/blast/db/). 
```perl
perl /home/dai.yun/ncbi-blast-2.12.0+/bin/update_blastdb.pl --passive --decompress refseq_select_prot &
```
If `update_blastdb.pl` does not work (this happens when a firewall is in place), the db files can be downloaded manually using commands like `wget`. The `metadata.json` files include URLs that can be used to download the compressed binary files (.tar.gz). The downloaded files must be decompressed manually using `tar -xf`.   

### 2. Perform protein blast locally (BLASTp)
```
blastp -db refseq_select_prot -query LdcA.fasta -evalue 0.0001 -max_target_seqs 20000 -outfmt 7 -out blastp_LdcA_result_eval_1e-4_max20000_outfmt7.txt &
```
### 3. Extract protein IDs from BLASTp result
Next, save the protein IDs in a separate file. In the example below, the first record starts at line 6 (1-5 being info and header lines) - if your file is different, change the parameter following `tail` command.
```
cat blastp_LdcA_result_eval_1e-4_max20000_outfmt7.txt | tail -n +6 | cut -f2 > list_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt
```

### 4. Get protein sequences from NCBI
Get the protein sequences for each ID in the list generated in 3.
```
for line in $(cat list_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt); 
do efetch -db protein -id $line -format fasta >> seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.fasta; 
done &
```

### 5. Obtain predictions for signal peptide, transmembrane helices, and conserved domains
Use the fasta file from 4 as input, get prediction output files from the following databases:

1. [SignalP](http://www.cbs.dtu.dk/services/SignalP/): short output (no figures). Use [mirror site](https://services.healthtech.dtu.dk/service.php?SignalP) if the original site is overloaded.
2. [predisi](http://www.predisi.de/) (GN and GP results seem to be the same) Do not trim the header lines - `Filtering.py` will do that.
3. [phobius](https://phobius.sbc.su.se/): short output format (make sure to wait until the results page finishes loading, otherwise you will end up with incomplete result). The first column header must be manually changed from `SEQENCE ID` to `SEQUENCE_ID`.
4. [CDD domain](https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi) 
5. [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) Do not trim the header lines - `Filtering.py` will do that.

Note that the CDD batch search function does not allow fasta files with more than 4000 sequences. The following command splits a large fasta file into smaller ones. In the example below, the command splits `seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.fasta` (7333 sequences) into two files: `seq_blast_LdcA_results_eval_1e-4_max20000_outfmt7.00001.fasta` (4000 sequence) and `seq_blast_LdcA_results_eval_1e-4_max20000_outfmt7.04001.fasta` (3333 sequences).
```
awk -v size=4000 -v pre=seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7 -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.fasta
```

Use the `cat` command to concatenate multiple prediction outputs of smaller fasta files into one file.
```
cat file1 file2 > file_total
```
(When concatenating files, make sure all descriptor lines before the column hearder are removed.)


### 6. Exclude hits that does not meet certain structural/size criteria

Use `Filtering.py` to apply selected filters based on predicted domain/sequence features.

```
python /home/dai.yun/GeisingerLab/ElsL_Orthology_Analysis/src/Filtering.py \
-i /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.fasta \
-p /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/signalp-GP-seq_blast_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
-n /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/signalp-GN-seq_blast_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
-b /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/phobius_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
-c /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/CDD_domain_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
-t /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/TMHMM_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
-x /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/predisi_GN_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
-y /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/predisi_GP_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
-o /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/shortlist_prot_acc_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt
```
Note that `Filtering.py` the following python packages to be installed. Check package version if errors occur. 

### 7. Creat info table for the shortlist of proteins from NCBI.
For output generated in 6, use `make_shortlist_info.py` to extract taxonomy information from NCBI.
```
python /home/dai.yun/GeisingerLab/ElsL_Orthology_Analysis/src/make_shortlist_info.py \
/scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/shortlist_prot_acc_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
/scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/CDD_domain_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
/scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/shortlist_info_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt
```

### 8. Plot distribution of sequence lengths on a histogram.

If the input file is the info table generated in 7, the script will directly use the protein length listed in the table.
```
python /home/dai.yun/GeisingerLab/ElsL_Orthology_Analysis/src/seq_length_hist.py \
/scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/ldcA_blastp/predictions/shortlist_info_seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.txt \
/scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/histograms/hist_shortlist_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.pdf
```

If the input file is a fasta file, use the "fasta" option so that the script will determine the length for each protein sequence.
```
python /home/dai.yun/GeisingerLab/ElsL_Orthology_Analysis/src/seq_length_hist.py \
/scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcV_phylogenetics/seq_blastp_LdcV_result_eval_1e-4_max20000_outfmt7.fasta \
/scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/ldcA_phylogenetics/histograms/hist_unfiltered_blastp_LdcV_result_eval_1e-4_max20000_outfmt7.pdf \
fasta
```

### 9. Obtain taxon IDs from protein accessions

```
for line in $(cat /scratch/dai.yun/2021July_ElsL_LDC_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/shortlist_new/GT9TM664013_shortlist_prot_acc_new.txt); \
do efetch -db protein -id $line -format est | grep "taxon" | cut -d ":" -f2 | cut -d "\"" -f1; \
done >> taxid_ElsL_shortlist.txt &
```
