# Analyzing Chip-Seq Data
This is a collection of scripts for analyzing chip-seq data using various bioinformatic tools, such as `macs2` and `meme-suite`.
## Requirements
+ Python3 with following packages installed:
  + pandas
  + numpy
  + biopython (Bio)
+ FASTX-Toolkit 
+ Bowtie (version > 1.2.2)
+ MACS2
+ MEME-ChIP
### Folder Structure Used
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
