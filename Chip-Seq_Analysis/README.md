# Analyzing Chip-Seq Data
This is a collection of scripts for analyzing chip-seq data using various bioinformatic tools, such as `macs2` and `meme-suite`.
## Requirements
+ Python3 with the following packages installed. Follow these [instructions](https://packaging.python.org/tutorials/installing-packages/) on how to install python packages.
  + pandas
  + numpy
  + biopython (Bio)
+ [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html)
+ [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (version > 1.2.2)
+ [MACS2](https://github.com/macs3-project/MACS)
+ [MEME-ChIP](http://meme-suite.org/tools/meme-chip)
+ Sufficient disk space (>200G)
## Folder Structure Used
```bash
+-- Chip-Seq_Data
|   +-- gunzipped-FASTQ
    |   +-- BfmR-Chip-28_example1_r001.fastq.gz
    |   +-- BfmR-Chip-29_example2_r001.fastq.gz
    |   +-- BfmR-Chip-49_example3_r001.fastq.gz
|   +-- unzipped_fastq_files
    |   +-- BfmR-Chip-28_example1_r001.fastq
    |   +-- BfmR-Chip-29_example2_r001.fastq
    |   +-- BfmR-Chip-49_example3_r001.fastq
|   +-- clipped_fastq_files
    |   +-- clipped_BfmR-Chip-28_example1_r001.fastq
    |   +-- clipped_BfmR-Chip-29_example2_r001.fastq
    |   +-- clipped_BfmR-Chip-49_example3_r001.fastq
|   +-- bowtie_mapped_files
    |   +-- BfmR-Chip-28_example1_r001.bam
    |   +-- BfmR-Chip-29_example2_r001.bam
    |   +-- BfmR-Chip-49_example3_r001.bam
```
## Process Raw Data
### Access Sequencing Data from TUCF Server in Terminal
This example shows how to download data from the TUCF server. 
```bash
# Connect server with ftp command, followed by IP address of TUCF server
# Enter username and password when prompted
# Once login is successful, you will see >ftp at the beginning of a new command line
ftp 130.64.74.72
# Go to directory containing sequencing data (no need to type ftp)
>ftp cd [server directory]
# Turn off promt mode so that you don't have to enter "yes" for every file to download
>ftp prompt
# Retrieve all files in current directory using mget command
>ftp mget * > [local directory to save data]
# disconnect server
>ftp close
```
### Decompress Data
The sequencing data from TUCF are typically compressed in `.gz` format. We can unpack these files using the `gzip` command. 
```bash
# create a new folder for unzipped files
mkdir unzipped_fastq_files
gunzip -k ./gz/[file.fastq.gz]
# move unzipped file to the folder
mv ./gz/file.fastq ./unzipped_fastq_files/file.fastq
```
### Remove Barcodes with Fastx-tool
```bash
# create a new folder for clipped fastq files
mkdir unzipped_fastq_files
# discard sequences shorter than 20 bp using -l option
# use -v option to report number of sequences
# specify adapter sequence using -a
fastx_clipper -v -l 20 -a [adapter sequence] -i ./unzipped_fastq_files/file.fastq -o ./clipped_fastq_files/clipped_file.fastq
```
### Batch Unzip and Clip
To unzip and clip all the fastq.gz files in batch, use the unzip-clip_fastq.sh. 

## Sequence Alignment with Bowtie 1

## Find Peak with MACS2 peak caller
## Predict Motif using MEME-ChIP
