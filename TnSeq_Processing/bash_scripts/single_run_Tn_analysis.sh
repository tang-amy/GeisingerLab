## Yunfei Dai
## 2021/12/08

#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --job-name=tnseq_single
#SBATCH --mem=200G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_test.%N.%j.out
#SBATCH -e /work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis/slurm_log/slurm_test.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=y.dai@northeastern.edu

#Load modules and packages
module load python/3.7.1
module load bowtie/1.3.0
module load gcc/4.8.5
source /home/dai.yun/yunfei_python_venv/bin/activate

GZ=${1?Error: "requires directory of a gz file"}
OUT_DIR=/work/geisingerlab/Yunfei/2021Dec_dacC_Tnseq_analysis
EMPTY_WIG=/work/geisingerlab/Yunfei/Ab_reference/17978-mff_empty.wig
REF=/work/geisingerlab/Yunfei/Ab_reference/AB17978/AB17978
WIG_SCRIPT=/home/dai.yun/GeisingerLab/TnSeq_Processing/src/exact-match_BC_populate_empty_wig_SAM.py
CLIPPER=/work/geisingerlab/software/fastx_toolkit/bin/fastx_clipper
ADAPTOR=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

#Make directories if missing.
mkdir -p $OUT_DIR/raw_fastq
mkdir -p $OUT_DIR/clipped_fastq
mkdir -p $OUT_DIR/SAM
mkdir -p $OUT_DIR/WIG_exact_match

#Execute pipeline
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

echo "Finished pipeline for ${GZ}."
