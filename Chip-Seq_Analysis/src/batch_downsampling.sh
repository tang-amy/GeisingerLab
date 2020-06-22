#!/bin/bash

## Yunfei Dai
## 04/28/2020

<<Comment
  downsample all .bam files in directory ./sorted_bam/ for each file,
  3 downsampled bam files containing n (default is 100,000) reads are created, using 3 seed numbers (2,5,9)
  exact read number from each downsampled files are counted and writted into downsampling_log.txt

  Usasge:
  bash batch_downsampling.sh [INPUT Directory] 
Comment

DIR=${1?Error: "requires directory of bam files"}
read -p "Please enter the target read number:" READS
OUT=$(dirname ${DIR})/$(echo "scale=1; (${READS}/1000)" | bc)k_downsampled_BAM_files

for seed in 2 5 9
do
  mkdir -p ${OUT}/seed${seed} ${OUT}/seed${seed}/input;
  for file in $DIR/*.bam
  do
    read_count=$(samtools view ${file} | wc -l);
    scale=$(echo "scale=6; (${READS}/${read_count})+${seed}" | bc);
    rt=$(echo "scale=1; (${READS}/1000000)" | bc)
    fname=$(basename ${file});
    oname=${fname/sorted.bam/s${seed}_${rt}M.sorted.bam};
    ofile=${OUT}/seed${seed}/${oname};
    samtools view -s $scale -b $file -o $ofile;
    echo ${fname} >> ${OUT}/downsampling_${rt}M_log.txt;
    echo "s${seed}_${rt}M: $(samtools view $ofile | wc -l) reads" >> ${OUT}/downsampling_${rt}M_log.txt;
    done
  mv ${OUT}/seed${seed}/*Input*.bam ${OUT}/seed${seed}/input/
  done

printf "Finished downsampling. \n"
