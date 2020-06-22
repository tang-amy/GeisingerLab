#!/bin/bash

DIR=${1?Error: "enter directory of fastq files"}
REF=${2?Error: "enter directory of reference files (.ebwt), e.g. Ab17978/Ab17978"}
SAM_OUTPUT="$(dirname $DIR)/mapped_SAM_files"
BAM_OUTPUT="$(dirname $DIR)/sorted_BAM_files"

mkdir -p ${SAM_OUTPUT} ${BAM_OUTPUT}

parallel bowtie -m 1 -n 1 --best -y -S $REF {} {.}.sam ::: ${DIR}/*.fastq
printf "Finished sequence alignment.\n"

for file in ${DIR}/*.sam
do fname=$(basename ${file})
mv $file ${SAM_OUTPUT}/${fname/clipped_}
done

for file in ${SAM_OUTPUT}/*.sam
do fname=$(basename ${file})
samtools sort $file -o ${BAM_OUTPUT}/${fname/.sam/_sorted.bam}
done
printf "Finished sorting SAM files.\n"
