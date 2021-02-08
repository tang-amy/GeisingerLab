#!/home/bin/bash

# assuming in directory of zipped fastq files
# will make all directories/files in parent directory of zipped files

DIR=${1:-$PWD}
read -p "Please enter adaptor sequence:" ADAPTOR
UNZIPPED="$(dirname $DIR)/unzipped_fastq_files"
CLIPPED="$(dirname $DIR)/clipped_fastq_files"
mkdir -p $UNZIPPED $CLIPPED

# unzip all .gz files in $DIR
<<<<<<< HEAD
parallel gunzip -c {} ::: $DIR/*.gz
=======
# removed -k option of gunzip to adapt to RC
parallel gunzip -d {} ::: $DIR/*.gz
>>>>>>> b4fd803458db7cac6d0599a0ebd5fbe3a7db9a84
printf "Finished unzipping files, now clipping adaptors.\n"
printf "Adapator sequence $ADAPTOR will be removed from reads.\n"
for file in $DIR/*.fastq
do mv $file $UNZIPPED
done

# clips adaptors from *.fastq files in #UNZIPPED
# adaptor used for Tn-seq: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
cd $UNZIPPED
parallel fastx_clipper -v -l 20 -d 0 -Q 33 -a $ADAPTOR -i {} -o ${CLIPPED}/clipped_{} ::: *.fastq
printf "Finished clipping adaptors.\n"
cd $DIR
