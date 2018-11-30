Transposon-Sequencing Data Processing
=============================
_On: MacOS, Linux_
</br>
_Requirements: Bash, Python 2.7_

Installing FastX-Tools
-------------------------------
### Linux ###

### MacOS ###



Installing BowTie
-------------------------------
### Linux ###

### MacOS ###



Running TnSeqProcessing Script
-------------------------------
<i> Before starting, install and build FastX-Tools and Bowtie accordingly. <i>

<p> This script processes gunzipped FastQ files into Wiggle files that can be used in the Transit pipeline as well as a merged Wiggle file named "merged_wig.wig" containing all of the reads from every processed file. The merged-wig file can be used to view the read counts on IGV. </p>

1. Download the TnSeqProcessing scripts from (add link here) or by cloning this repository. 
2. Navigate to the src folder containing TnSeqProcessing scripts (TnseqProcessing/src) with the terminal. Change the scripts to be executable with the following command:
```bash
sudo chmod 777 *
```
3. To run the scripts, go to the src folder containing the scripts and run the following command in terminal:
```bash
bash tnseq-process.sh -i [input folder containing gunzipped files] -g [genbank file for the genome (.gbk) file] -r [reference indices for Bowtie]
```
The required inputs are the folder containing gunzipped FASTQ files and the Genbank genome file (.gbk). The reference indices should be provided as PATH/reference-index-folder/reference-index-name (skipping the .ebwt portion of the name). 

(If no reference indices are given, the script will ask for a FASTA file to provided in order to create the reference indices for Bowtie processing.) 

4. Folder structure/outputs:
<b> For best results, nest the folder containing the gunzipped FASTQ files inside of another folder for all of the outputs/directories made by the script to reside in. </b>
```bash
+-- TnSeq_Data
|   +-- clipped_fastq_files
    |   +-- clipped_example1_r001.fastq
    |   +-- clipped_example2_r001.fastq
    |   +-- clipped_example3_r001.fastq
|   +-- map_files
    |   +-- clipped_example1_r001.map
    |   +-- clipped_example2_r001.map
    |   +-- clipped_example3_r001.map
|   +--merged_wig.wig
|   +-- gunzipped-FASTQ
    |   +-- example1_r001.fastq.gz
    |   +-- example2_r001.fastq.gz
    |   +-- example3_r001.fastq.gz
|   +-- unzipped_fastq_files
    |   +-- example1_r001.fastq
    |   +-- example2_r001.fastq
    |   +-- example3_r001.fastq
|   +-- wig_files
    |   +-- clipped_example1_r001.wig
    |   +-- clipped_example2_r001.wig
    |   +-- clipped_example3_r001.wig
```


