#!/bin/bash

## Yunfei Dai
## 05/01/2020


<<Comment
  Find consensus peaks from macs2 narrowPeak outputs of replicates using bedtool multiIntersectBed function.
  Usage:
  bash batch_call_peak.sh <input directory> <output directory>
Comment


DIR=${1?Error: "requires directory of bam files"}
# do not include "/consensus_fasta_bed_intersect/"
OUT=${2?Error: "requires output directory"}

mkdir -p ${OUT}/bed_intersect_consensus_meme_chip/fimo_output ${OUT}/bed_intersect_consensus_meme_chip/html_output

for file in $DIR/consensus_fasta_bed_intersect/*.fasta
  do
  fname=$(basename $file)
  meme-chip -meme-p 4 \
  -o ${OUT}/bed_intersect_consensus_meme_chip/${fname%.fasta} \
  -db /Users/yunfei/meme/motif_databases/PROKARYOTE/collectf.meme $file >> meme-chip_log.txt 2>&1;
done

# copy all fimo outputs into a single folder
for file in ${OUT}/bed_intersect_consensus_meme_chip/*/fimo_out*/*.gff
do
  fname=$(echo $file | rev | cut -d'/' -f3 | rev);
  ext=$(echo $file | rev| cut -d'/' -f2 | rev);
  cp -a $file ${OUT}/bed_intersect_consensus_meme_chip/fimo_output/${fname}.${ext}.gff;
done

# copy all html outputs into a single folder
for file in ${OUT}/bed_intersect_consensus_meme_chip/*/*.html
do
  fname=$(echo $file | rev | cut -d'/' -f2 | rev);
  cp -a $file ${OUT}/bed_intersect_consensus_meme_chip/html_output/$fname.meme-chip.html;
done
