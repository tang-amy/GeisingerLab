
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
2. [predisi](http://www.predisi.de/) (GN and GP results seem to be the same)
3. [phobius](https://phobius.sbc.su.se/): short output format (make sure to wait until the results page finishes loading, otherwise you will end up with incomplete result)
4. [CDD domain](https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi) 
5. [TMHMM](http://www.cbs.dtu.dk/services/TMHMM/)

Note that the CDD batch search function does not allow fasta files with more than 4000 sequences. The following command splits a large fasta file into smaller ones. In the example below, the command splits `seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.fasta` (7333 sequences) into two files: `seq_blast_LdcA_results_eval_1e-4_max20000_outfmt7.00001.fasta` (4000 sequence) and `seq_blast_LdcA_results_eval_1e-4_max20000_outfmt7.04001.fasta` (3333 sequences).
```
awk -v size=4000 -v pre=seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7 -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' seq_blastp_LdcA_results_eval_1e-4_max20000_outfmt7.fasta
```

Use the `cat` command to concatenate multiple prediction outputs of smaller fasta files into one file.
```
cat file1 file2 > file_total
```
