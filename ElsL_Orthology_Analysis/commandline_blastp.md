
### 1. Download preformatted blast database from NCBI

Use the `update_blastdb.pl` script in the blast+ package to download the preformatted blast database files available in the [ncbi ftp site](https://ftp.ncbi.nlm.nih.gov/blast/db/). 
```perl
perl /home/dai.yun/ncbi-blast-2.12.0+/bin/update_blastdb.pl --passive --decompress refseq_select_prot &
```
