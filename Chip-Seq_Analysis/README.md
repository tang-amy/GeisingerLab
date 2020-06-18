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
|   +-- gunzipped-FASTQ (gz)
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
The sequencing data from TUCF are typically compressed in `.gz` format. We can unpack these files using the `gunzip` command. 
Usage for individual files:
```bash
# create a new folder for unzipped files
mkdir unzipped_fastq_files
gunzip -k ./gz/[file.fastq.gz]
# move unzipped file to the folder
mv ./gz/file.fastq ./unzipped_fastq_files/file.fastq
```
### Remove Barcodes with Fastx-tool
Once the fastq files are available, we will remove the adaptors from the reads and discard reads that are too short (less than 20 bp for example). 
Usage for individual files:
```bash
# create a new folder for clipped fastq files
mkdir unzipped_fastq_files
# discard sequences shorter than 20 bp using -l option
# use -v option to report number of sequences
# specify adapter sequence using -a
fastx_clipper -v -l 20 -a [adapter sequence] -i ./unzipped_fastq_files/file.fastq -o ./clipped_fastq_files/clipped_file.fastq
```
### Batch Unzip and Clip
To unzip and clip all the fastq.gz files in batch, use bash script [batch_unzip_clip_fastq.sh](https://github.com/tang-amy/GeisingerLab/blob/master/Chip-Seq_Analysis/src/batch_unzip_clip_fastq.sh). For a given directory (default is current directory if no input entered), this script generates two folders `unzipped_fastq_files` and `clipped_fastq_files` in the parent directory, unzips files in the given directory and move them to `unzipped_fastq_files`, then clip the adaptors from the fastq files, saving output as `clipped_file.fastq` in `clipped_fastq_files`. This process may take a while depending on your computer specs.

```bash
# the unzipped_fastq_files clipped_fastq_files folders will be generated automatically, no need to creat them in advance
# make sure the path to the script is correct
bash batch_unzip_clip_fastq.sh ./gz
# the following message will prompt, enter the sequence (no quotation marks or brackets) following the message
Please enter adaptor sequence:
# list files in the unzipped_fastq_files folder
ls unzipped_fastq_files
# list files in the check the clipped_fastq_files folder
ls clipped_fastq_files
```
Note that the script uses GNU `Parallel` for faster execution of all commands. If parallel is not installed, run the following command in terminal. 
```bash
# check if Parallel is installed
which parallel
# install Parallel
brew install parallel
```
If brew is not installed, check [Homebrew](https://brew.sh/) website for installation instructions.

## Sequence Alignment with Bowtie 1

## Find Peak with MACS2 peak caller
## Predict Motif using MEME-ChIP
