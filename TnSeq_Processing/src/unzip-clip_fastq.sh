#!/home/bin/bash

# assuming in directory of zipped fastq files
# will make all directories/files in parent directory of zipped files

DIR=${1:-$PWD}

# creating var for output directories
UNZIPPED="$(dirname $DIR)/unzipped_fastq_files"
CLIPPED="$(dirname $DIR)/clipped_fastq_files"

#unzipping files
#keeps original gz zipped file
mkdir -p $UNZIPPED
for file in $DIR/*.gz
do gunzip -k $file
done
for file in $DIR/*.fastq
do mv "$file" $UNZIPPED
done

echo "Finished unzipping all .gz files into .fastq files."

#clips files and adds clipped files to its own directory
mkdir -p $CLIPPED
for file in $UNZIPPED/*.fastq
do
fastx_clipper -v -l 20 -a [adapter sequence] -i $file -o $CLIPPED/clipped_$(basename $file)
# parallel 'fastx_clipper -v -l 20 -a [adapter sequence] -i {} -o ${CLIPPED}/clipped_{}' ::: $UNZIPPED/*.fastq
done

echo "Finished clipping all .fastq files with fastx_clipper."
