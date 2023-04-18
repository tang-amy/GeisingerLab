#!/bin/bash

## Yunfei Dai
## 17Apr2023

#SBATCH --time=08:00:00
#SBATCH --job-name=macs2
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/slurm_macs2.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/slurm_macs2.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dai.yun@northeastern.edu

module load anaconda3/2021.05
source activate YD_chipseq_venv

SCRIPT=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Scripts/batch_run_macs2.py
INDIR=/work/geisingerlab/Yunfei/2022_Nov_ChipSeq_Analysis/BAM_Downsampled_500k
OUTDIR=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/MACS2_output_500k_ext200_input28-2

mkdir -p ${OUTDIR}/narrowPeak ${OUTDIR}/xls ${OUTDIR}/bed

python $SCRIPT $INDIR $OUTDIR

mv ${OUTDIR}/*.narrowPeak ${OUTDIR}/narrowPeak
mv ${OUTDIR}/*.xls ${OUTDIR}/xls
mv ${OUTDIR}/*.bed ${OUTDIR}/bed

