#!/bin/bash

## Yunfei Dai
## 11/05/2020

DIR=${1?Error: "enter directory of fastq files"}
OUT=${2?Error: "provide output directory name"}
mkdir -p ${OUT}/narrowPeak ${OUT}/xls ${OUT}/bed

for i in 1 2 3
        do for file in ${DIR}/JBA71FLAG_BEADS_${i}_*.bam
                do fname=$(basename $file); 
                oname=${OUT}/${fname%%.sorted.bam};
                input=${DIR}/JBA71_INPUT_${i}_*.bam;
                macs2 callpeak -t $file -c $input -g 3.8e6 \
                -n $oname.JBA71_INPUT_${i}.ext200 \
                --nomodel --extsize 200 >> ${OUT}/macs2_log.txt 2>&1;
                done
         done


mv ${OUT}/*.narrowPeak ${OUT}/narrowPeak
mv ${OUT}/*.xls ${OUT}/xls
mv ${OUT}/*.bed ${OUT}/bed
