#!/home/bin/bash

# assuming in directory of clipped fastq files if no directory is entered
DIR=${1?Error: "requires directory of fastq files"}
REF=${2?Error: "requires reference"}

ODIR="$(dirname $DIR)/map_files"
mkdir -p ${ODIR}

for file in $DIR/*.fastq
do
bowtie -m 1 -n 1 --best -y $REF $file "$ODIR/$(basename ${file/%fastq/map})"
done

# parallel 'bowtie -m 1 -n 1 --best -y Ab17978 {} ${ODIR}/{}/%fastq/map' ::: ${DIR}/*.fastq.gz

# ${file/%fastq/map}"'