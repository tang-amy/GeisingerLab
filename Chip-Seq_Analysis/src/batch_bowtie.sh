#!/bin/bash

DIR=${1?Error: "enter directory of fastq files"}
REF=${2?Error: "enter directory of reference files (.ebwt), e.g. Ab17978/Ab17978"}
SAM_OUTPUT="$(dirname $DIR)/mapped_SAM_files"
BAM_OUTPUT="$(dirname $DIR)/sorted_BAM_files"

mkdir -p ${OUTPUT}

parallel bowtie -m 1 -n 1 --best -y -S $REF {} {.}.sam ::: ${DIR}/*.fastq

for file in ${DIR}/*.sam
do mv $file ${OUT}/${file/clipped}
done
