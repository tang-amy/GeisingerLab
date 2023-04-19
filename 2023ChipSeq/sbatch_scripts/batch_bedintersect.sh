#!/bin/bash

## Yunfei Dai
## 18Apr2023

#SBATCH --time=08:00:00
#SBATCH --job-name=bedintersect
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/bedintersect.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/bedintersect.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dai.yun@northeastern.edu

module load anaconda3/2021.05
source activate YD_chipseq_venv

INDIR=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/MACS2_output_500k_ext200_input28-2/narrowPeak
OUTDIR=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/MACS2_output_500k_ext200_input28-2/bed_multiinter
SCRIPT=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Scripts/batch_bedintersect.py

mkdir -p $OUTDIR
python $SCRIPT $INDIR $OUTDIR