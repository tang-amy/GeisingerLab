
### 1. Download preformatted blast database from NCBI

On bash command line terminal, use the `update_blastdb.pl` script [(user manual)](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory) in the blast+ package to download the preformatted blast database files available in the [NCBI ftp site](https://ftp.ncbi.nlm.nih.gov/blast/db/). 
```perl
perl /home/dai.yun/ncbi-blast-2.12.0+/bin/update_blastdb.pl --passive --decompress refseq_select_prot &
```
If `update_blastdb.pl` does not work, the db files can be downloaded manually using commands like `wget`. The `metadata.json` files include URLs that can be used to access the compressed files (.tar.gz). The downloaded files must be decompressed manually using `tar -xf`.   
