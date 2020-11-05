# Analyze Chip-Seq Data
This is a collection of scripts for analyzing chip-seq data using various bioinformatic tools, such as `macs2` and `meme-suite`.
## Requirements
+ Python3 with the following packages installed. Follow these [instructions](https://packaging.python.org/tutorials/installing-packages/) on how to install python packages.
  + pandas
  + numpy
  + biopython (Bio)
+ [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html)
+ [Samtools](http://www.htslib.org/) (see also: https://github.com/samtools/samtools)
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
|   +-- mapped_SAM_files
    |   +-- BfmR-Chip-28_example1_r001.sam
    |   +-- BfmR-Chip-29_example2_r001.sam
    |   +-- BfmR-Chip-49_example3_r001.sam
|   +-- sorted_BAM_files
    |   +-- sorted_BfmR-Chip-28_example1_r001.bam
    |   +-- sorted_BfmR-Chip-29_example2_r001.bam
    |   +-- sorted_BfmR-Chip-49_example3_r001.bam
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
# Switch to binary mode
>ftp binary
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
# brackets mean the full path to a file
# do not include the brackets when typing the commands
gunzip -k ./gz/[file.fastq.gz]
# move unzipped file to the folder
mv ./gz/file.fastq ./unzipped_fastq_files/file.fastq
```
### Remove Barcodes with Fastx-tool
Once the fastq files are available, we will remove the adaptors from the reads and discard reads that are too short (less than 20 bp for example) as quality control. 
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
# list files in the clipped_fastq_files folder
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

## Sequence Alignment
### Sequence Alignment with Bowtie1
The short reads in the fastq files need to be mapped to the genome (sequence alignment) prior to all downstream analyses. There are different sequence alignment softwares such as BWA and bowtie. Here we use bowtie 1. Check [here](http://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie) for installation instructions. Bowtie requires Bowtie index files, which are basically transformed from the genome sequence, and serve as the reference for the alignment. Bowtie index can be built with the fasta genome file using bowtie-build indexer.
Use the following command to build the index. Note that you only need to do this once for one genome.
```bash
# [genome.fasta] can be downloaded from NCBI database
# [ebwt_base] is the prefix for the index files (e.g. Ab17978)
bowtie-build [genome.fasta] [ebwt_base]
# the output is a folder (Ab17978) containing multiple .ebwt files that share the prefix Ab17978
```
To perform genome alignment for one fastq files:
```bash
# Ab17978 is the prefix for index files
# output is .sam file
bowtie -m 1 -n 1 --best -y -S Ab17978/Ab17978 [input.fastq] [output.sam]
```
### Convert to Sorted BAM file
The output from bowtie is [SAM file format](http://www.htslib.org/doc/sam.html). For the ease of downstream analyses, we need to sort the reads, and save the results as BAM file format, which is a binary version of SAM. SAM files can be sorted by `samtools sort`.
```bash
# make sure samtools is installed
which samtools
# sort the sam file and save output as bam file
samtools sort [input.sam] -o [output.bam]
```
### Batch Bowtie and Sort
For sequence alignment of multiple fastq files, use [batch_bowtie.sh](https://github.com/tang-amy/GeisingerLab/blob/master/Chip-Seq_Analysis/src/batch_bowtie.sh). For a given directory, this script generates two folders `mapped_SAM_files` and `sorted_BAM_files` in the parent directory and maps the fastq reads to the reference. It also moves the SAM output files from bowtie to `mapped_SAM_files`, sorts them with `samtools sort` and save the output as `file_sorted.bam` in `sorted_BAM_files`.

```bash
# the output folders will be generated automatically, no need to creat them in advance
bash batch_bowtie.sh [clipped_fastq] [Ab17978/Ab17978]
# list files in the mapped_SAM_files folder
ls mapped_SAM_files
# list files in the sorted_BAM_files folder
ls sorted_BAM_files
```

## Find Peak with MACS2 peak caller
### Install MACS2
Follow these instructions on [how to install MACS2](https://github.com/macs3-project/MACS/blob/master/INSTALL.md). "Install from source" is recommended.
### Downsample BAM reads

### Batch downsampling
[batch_downsampling.sh](https://github.com/tang-amy/GeisingerLab/blob/master/Chip-Seq_Analysis/src/batch_downsampling.sh)
### Run MACS2
```bash
macs2 callpeak -t $file -c $input -g 3.8e6 \
-n $oname.input_28-$i.ext200 \
--nomodel --extsize 200 >> $macs2_log.txt 2>&1
```
To call peaks from multiple bam files in a folder, use [batch_call_peak.sh](https://github.com/tang-amy/GeisingerLab/blob/master/Chip-Seq_Analysis/src/batch_call_peak.sh).

To call peaks from downsampled datasets, use [batch_call_peak_downsampled.sh](https://github.com/tang-amy/GeisingerLab/blob/master/Chip-Seq_Analysis/src/batch_call_peak_downsampled.sh)

## Predict Motif using MEME-ChIP
