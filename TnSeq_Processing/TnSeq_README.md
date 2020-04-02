# Tn-Seq Processing README

Includes scripts for processing FASTQ sequencing files into WIG (wiggle) files. Scripts written in Bash (shell scripting) and Python 3. For running TRANSIT (DeJesus, M. A., et al, 2015), please see their [documentation](https://transit.readthedocs.io/en/latest/tpp.html).

## Requirements
- Command Line Modules/Tools
  - FASTX_ToolKit from Hannon Lab (from [Gordon, A. & Hannon, G. J., 2010](http://hannonlab.cshl.edu/fastx_toolkit/))
  - Bowtie (from [Langmead, B., et al, 2010](http://bowtie-bio.sourceforge.net/index.shtml))
- Python 3 with biopython (Bio) package installed
- Reference files
  - Genomic FASTA file (eg. NZ_CP012004.fasta)
  - Genbank file (eg. NZ_CP012004.1.gb or NZ_CP012004.1.gbk)

## Notes
To be able to run these scripts on command line, they must be made executable. Use the following commands to change the rights for the scripts.
  ```bash
  chmod u+x script.py 
  chmod u+x script.sh
  ```

## Processing Steps
1. Unzip the FASTQ files using gunzip. Then clip FASTQ sequencing files using FASTX_ToolKit in Terminal. Use [unzip-clip_fastq.sh](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/src/unzip-clip_fastq.sh) to unzip/clip all zipped (.gz) files into FASTQ files (.fastq). The script will create 2 folders, unzipped and clipped, for the processed FASTQ files to live in.
```bash
# to run batch script 
# run in directory where zipped fastq (.gz) files are
./unzip-clip_fastq.sh

# to run on single file
gunzip -k zippedFile
fastx-clipper -l 20 -a "adapter_sequence" -i inFile -o OutFile
```
2. Build a reference genome from Bowtie if one hasn't been built yet. For Step 3, the reference directory should be enterd as GenomeName/GenomeName (what was used in this step). 
```bash
mkdir -p GenomeName
bowtie-build -f reference.fasta GenomeName/GenomeName
```
3. Map FASTQ files with Bowtie to genome. Use [run_bowtie_batch.sh](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/src/run_bowtie_batch.sh) to map all files in the given directory to the given Bowtie reference. The script will create a map_files folder for the outputted MAP files (renamed SAM files) to live in.
```bash
# to run the batch script
./run_bowtie_batch.sh clipped_fastq_dir reference_files

# to run on a single file
bowtie -m 1 -n 1 --best -y -S BowtieReference/BowtieReference yourfastqfile.fastq outputFile.map
```
4. Create the WIG files from the MAP files using either [map_to_wig.py](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/src/map_to_wig.py) or [batch_map2wig.py](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/src/batch_map2wig.py). Use full names for the pathing to the input/output files/directories or run in the appropriate directory. Note: batch_map2wig.py requries map_to_wig.py.
```bash
# to run the batch script
python ./batch_map2wig.py -g reference.gbk -i mapped_files_dir -o wig_files

# to run on a single file
python ./map_to_wig.py -g reference.gbk -i map_file.map -o wig_file.wig
```
5. Create prot table for TRANSIT analysis if necessary by running [make_prot_table.py](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/src/make_prot_table.py).
```bash
python ./make_prot_table.py -i reference.gbk -o outFile.prot_table
```
6. To merge all of the WIG files from one TnSeq library, run [merge_wigs.py](https://github.com/tang-amy/GeisingerLab/blob/master/TnSeq_Processing/src/merge_wigs.py). Once again, use full names for the pathing to the input/output files/directories or run in the appropriate directory.
```bash
python ./merge_wigs.py -i wig_file_dir -o merged_wig.wig -g reference.gbk
```

## Suggested Folder Structure
Processing the TnSeq data usually tends to generate a lot of files. This is the folder structures that many of the batch scripts will help generate to make ensure all file types in one folder are the same for organization.
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
