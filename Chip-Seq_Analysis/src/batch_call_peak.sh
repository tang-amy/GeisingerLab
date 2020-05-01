#!/bin/bash

## Yunfei Dai
## 04/28/2020

<<Comment
  Run macs2 with downsample treatment files against input.
  Usage:
  bash batch_call_peak.sh <input directory> <output directory>
Comment

DIR=${1?Error: "requires directory of bam files"}
OUT=${2?Error: "provide output directory name"}
mkdir ${OUT}/narrowPeak ${OUT}/xls ${OUT}/bed
for seed in 2 5 9
do for i in 1 2 3
	do for file in ${DIR}/seed$seed/*ChIP-??-$i*
		do fname=$(basename $file);
		oname=${OUT}/${fname%%.sorted.bam};
		input=${DIR}/seed$seed/input/*28-$i*
		macs2 callpeak -t $file -c $input -g 3.8e6 -n $oname.input_28-$i.ext200 --nomodel --extsize 200 >> ${OUT}/macs2_log.txt 2>&1;
		done	
	done
done
mv ${OUT}/*.narrowPeak ${OUT}/narrowPeak
mv ${OUT}/*.xls ${OUT}/xls
mv ${OUT}/*.bed ${OUT}/bed