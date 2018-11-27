#!/home/bin/bash

# assuming in directory of zipped fastq files
# will make all directories/files in parent directory of zipped files

DIR=${1:-$PWD}

# creating var for output directories
cd ${DIR}
UNZIPPED="$(dirname $PWD)/unzipped_fastq_files/"
CLIPPED="$(dirname $PWD)/clipped_fastq_files/"

#unzipping files
#keeps original gz zipped file
mkdir '${UNZIPPED}'
parallel 'gunzip -k {}' ::: *.gz
#moves unzipped files to new directory
for file in ${DIR}/*.fastq
do mv "$file" "${UNZIPPED}"
done

echo "Finished unzipping all .gz files into .fastq files."

#clips files and adds clipped files to its own directory
mkdir "${CLIPPED}"
parallel 'fastx_clipper -i {} -o ${CLIPPED}/clipped_{}' ::: ${UNZIPPED}/*.fastq

echo "Finished clipping all .fastq files with fastx_clipper."