## Yunfei Dai
## 2021/12/08

#!/bin/bash
#SBATCH --array=0-12%13
#SBATCH --ntasks=13
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --job-name=tnseq
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_tnseq.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_tnseq.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=y.dai@northeastern.edu

#Loading modules and venv
module load anaconda3
module load bowtie/1.3.0
module load gcc/4.8.5
source activate /work/geisingerlab/conda_env/tnseq_env

#Directories
FASTA_DIR=/work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/raw_data/211203-0504H_Edward_Geisinger_transposed/fastq_Lane8
OUT_DIR=/work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis
EMPTY_WIG=/work/geisingerlab/Yunfei/Ab_reference/17978-mff_empty.wig
REF=/work/geisingerlab/Yunfei/Ab_reference/Ab_all/Ab_all
WIG_SCRIPT=/work/geisingerlab/GeisingerLab/TnSeq_Processing/src/exact-match_BC_populate_empty_wig_SAM.py
CLIPPER=/work/geisingerlab/software/fastx_toolkit/bin/fastx_clipper
ADAPTOR=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

#Create output folders
mkdir -p $OUT_DIR/raw_fastq
mkdir -p $OUT_DIR/clipped_fastq
mkdir -p $OUT_DIR/SAM
mkdir -p $OUT_DIR/WIG_exact_match

#Array of inputs
INPUT=( $(find $FASTA_DIR -type f -name "*.gz") )
func () {
    GZ=$1
    fname=$(basename $GZ)
    echo "processing ${GZ}"
    FASTQ=$OUT_DIR/raw_fastq/${fname%%.gz}
    CLIPPED=$OUT_DIR/clipped_fastq/${fname%%.fastq.gz}.clipped.fastq
    SAM_OUT=$OUT_DIR/SAM/${fname%%.fastq.gz}.sam
    WIG_OUT=$OUT_DIR/WIG_exact_match/${fname%%.fastq.gz}.wig
    gzip -cd $GZ > $FASTQ
    $CLIPPER -v -l 20 -d 0 -Q 33 -a $ADAPTOR -i $FASTQ -o $CLIPPED
    bowtie -m 1 -n 1 --best -y -S $REF $CLIPPED $SAM_OUT
    python $WIG_SCRIPT -i $EMPTY_WIG -s $SAM_OUT -o $WIG_OUT
}

#Lauch job array
func "${INPUT[$SLURM_ARRAY_TASK_ID]}"
