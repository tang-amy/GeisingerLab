#!/bin/bash

## Yunfei Dai
## 04/28/2020

<<Comment
  Find consensus peaks from macs2 narrowPeak outputs of replicates using bedtool multiIntersectBed function.
  Usage:
  bash batch_find_consensus_bedIntersect.sh <input directory>
Comment

DIR=${1?Error: "requires directory of bam files"}
# do not include "/narrowPeak"

for seed in 2 5 9
do
  mkdir -p ${DIR}/consensus_fasta_bed_intersect/concensus_peak/seed$seed;
  for sample in 28 49
	do
    replicates=$(ls ${DIR}/narrowPeak/BfmR-ChIP-${sample}*s${seed}*.narrowPeak);
    outname=BfmR-ChIP-${sample}_seed${seed}.consensus_peak.bed;
    outfile=${DIR}/consensus_fasta_bed_intersect/concensus_peak/seed$seed/${outname};
    echo "concensus found for $replicates";
    multiIntersectBed -cluster -i $replicates | awk -F"\t" '$4>1' > $outfile &&
    python /Users/yunfei/GeisingerLab/Chip-Seq_Analysis/src/getfasta.py \
    -i $outfile \
    -g /Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.fasta \
    -s 'consensus' \
    -o ${DIR}/consensus_fasta_bed_intersect/${outname/bed/fasta};
	done
done

mkdir -p ${DIR}/consensus_fasta_bed_intersect/intermediate_bed;
mv ${DIR}/consensus_fasta_bed_intersect/*.bed ${DIR}/consensus_fasta_bed_intersect/intermediate_bed