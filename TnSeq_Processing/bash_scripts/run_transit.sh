#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=transit
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_transit.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_transit.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=y.dai@northeastern.edu


module load python/3.7.1
source /work/geisingerlab/python_venv/transit_venv/bin/activate

INPUT_DIR=/work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/WIG_exact_match
OUT_DIR=/work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis
CONTROL_DIR=/work/geisingerlab/Yunfei/Ab_reference/WT17978MarinerTnseqLibraries_ATAJB_SS83_WIGandPROT_TABLE
TRANSIT=/work/geisingerlab/software/transit/src/transit.py
PROT_TABLE=/work/geisingerlab/Yunfei/Ab_reference/WT17978MarinerTnseqLibraries_ATAJB_SS83_WIGandPROT_TABLE/NZ_CP012004.prot_table

CONTROL_WIGS=$(echo $(find $CONTROL_DIR -type f -name "*.wig") | tr " " ,)
INPUT_WIGS=$(echo $(find $INPUT_DIR -type f -name "EKL*.wig") | tr " " ,)
RESULT=$OUT_DIR/transit_result/EKL_resampling_result.txt
mkdir -p $OUT_DIR/transit_result
python $TRANSIT resampling \
    $CONTROL_WIGS \
    $INPUT_WIGS \
    $PROT_TABLE $RESULT \
    -s 10000 \
    -n TTR \
    -l \
    -iN 0 \
    -iC 10 \
    --ctrl_lib "AbC" \
    --exp_lib "EKL"

