#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=100G
#SBATCH --job-name=peak
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/peak.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/Slurm_log/peak.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dai.yun@northeastern.edu

## Yunfei Dai
## 19Apr2023

#Loading modules and venv
module load anaconda3/2021.05
source activate /work/geisingerlab/Yunfei/2022_Nov_ChipSeq_Analysis/YD_chipseq_venv

script=/work/geisingerlab/GeisingerLab/2023ChipSeq/src/find_peak_summit.py
INDIR=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/MACS2_output_500k_ext200_input28-2/bed_multiinter
OUTDIR=/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/MACS2_output_500k_ext200_input28-2/average_summit_19APR2023

mkdir -p $OUTDIR

for file in $INDIR/*.intersect.bed;
    do fname=$(basename $file);
    outfile=$OUTDIR/${fname/.intersect.bed}.average_summit.bed;
    python $script $file $outfile;
done