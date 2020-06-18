#!/bin/bash


# This script unzips all the .fastq.gz files in the given directory using the gunzip command.
# After unzipping, adaptors are removed from the reads in each file using fastx_clipper. 
# For faster processing, parallel command is utilized.


DIR=${1:-$PWD}
read -p "Please enter adaptor sequence:" ADAPTOR
UNZIPPED="$(dirname $DIR)/unzipped_fastq_files"
CLIPPED="$(dirname $DIR)/clipped_fastq_files"
mkdir -p $UNZIPPED $CLIPPED

parallel gunzip -k {} ::: $DIR/*.gz
printf "Finished unzipping files, now clipping adaptors.\n"
printf "Adapator sequence $ADAPTOR will be removed from reads.\n"
for file in $DIR/*.txt
do mv $file $UNZIPPED
done

parallel fastx_clipper -v -l 20 -a $ADAPTOR -i {} -o ${CLIPPED}/clipped_{} ::: $UNZIPPED/*.fastq
printf "Finished clipping adaptors.\n"
