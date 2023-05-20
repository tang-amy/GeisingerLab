#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=100G
#SBATCH --job-name=mastertable
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/master.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/master.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dai.yun@northeastern.edu

## Yunfei Dai
## 04May2023 

INDIR=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/MACS2_output_500k_ext200_input28-2/peak_stat_20May2023
OUTDIR=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/RNAseq_Chipseq_mastertable_20May2023
script=/work/geisingerlab/GeisingerLab/2023ChipSeq/src/combine_chipseq_rnaseq_tables.py

module load anaconda3/2021.05
source activate /work/geisingerlab/Yunfei/2022_Nov_ChipSeq_Analysis/YD_chipseq_venv

mkdir -p $OUTDIR

for file in $INDIR/*.tsv;
    do fname=$(basename $file);
    outfile=$OUTDIR/${fname/.nearest_orf.tsv}.master_table.tsv;
    python $script $file $outfile;
done
