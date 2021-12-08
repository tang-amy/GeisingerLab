#!/bin/bash
#SBATCH --array=0-11%12
#SBATCH --ntasks=13
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --job-name=WIG
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_wig.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_wig.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=y.dai@northeastern.edu

#Loading modules and venv
module load python/3.7.1
source /home/dai.yun/yunfei_python_venv/bin/activate

#Directories
SAM_DIR=/work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/SAM
OUT_DIR=/work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/WIG_exact_match
EMPTY_WIG=/work/geisingerlab/Yunfei/Ab_reference/17978-mff_empty.wig
WIG_SCRIPT=/home/dai.yun/GeisingerLab/TnSeq_Processing/src/exact-match_BC_populate_empty_wig_SAM.py

#Create output folders
mkdir -p $OUT_DIR

#Array of inputs
INPUT=( $(find $SAM_DIR -type f -name "EKL*.sam") )
func () {
    SAM=$1
    fname=$(basename $SAM)
    echo "processing ${SAM}"
    WIG_OUT=$OUT_DIR/${fname%%.sam}.wig
    python $WIG_SCRIPT -i $EMPTY_WIG -s $SAM -o $WIG_OUT
}

#Lauch job array
func "${INPUT[$SLURM_ARRAY_TASK_ID]}"
