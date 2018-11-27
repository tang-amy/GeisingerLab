#!/home/bin/bash

# assuming in directory of clipped fastq files if no directory is entered
DIR=${1:-$PWD}

cd ${DIR}

ODIR="$(dirname $PWD)/map_files/"
mkdir ${ODIR}



parallel 'bowtie -m 1 -n 1 --best -y Ab17978 {} ${ODIR}/{}/%fastq/map' ::: ${DIR}/*.fastq.gz

# ${file/%fastq/map}"'