
### 1. Download preformatted blast database from NCBI

On bash command line terminal, use the `update_blastdb.pl` script [(user manual)](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory) in the blast+ package to download the preformatted blast database files available in the [NCBI ftp site](https://ftp.ncbi.nlm.nih.gov/blast/db/). 
```perl
perl /home/dai.yun/ncbi-blast-2.12.0+/bin/update_blastdb.pl --passive --decompress refseq_select_prot &
```
If `update_blastdb.pl` does not work, the db files can be downloaded manually using commands like `wget`. The `metadata.json` files include URLs that can be used to download the compressed binary files (.tar.gz). The downloaded files must be decompressed manually using `tar -xf`.   

### 2. Perform protein blast locally (BLASTp)
```
blastp -db refseq_select_prot -query LdcA.fasta -evalue 0.0001 -max_target_seqs 20000 -outfmt 7 -out blastp_LdcA_result_eval_1e-4_max20000_outfmt7.txt &
```
### 3. Extract protein IDs from BLASTp result
Next, save the protein IDs in a separate file.
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
[SignalP](http://www.cbs.dtu.dk/services/SignalP/): short output (no figures)
[predisi](http://www.predisi.de/) 
[phobius](https://phobius.sbc.su.se/): short output format
[CDD domain](https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi)
