# Yunfei Dai
# 2021/02/05
# This script allows batch alignment of fastq files in a folder

#!/home/bin/bash

# enter directories

DIR=${1?Error: "enter directory of fastq files"}
REF=${2?Error: "enter directory of reference files (.ebwt), e.g. Ab17978/Ab17978"}
SAM_OUTPUT="$(dirname $DIR)/mapped_SAM_files"

mkdir -p ${SAM_OUTPUT}

parallel bowtie -m 1 -n 1 --best -y -S $REF {} {.}.sam ::: ${DIR}/*.fastq
printf "Finished sequence alignment.\n"

# sequence alignment for all sam files in $DIR
for file in ${DIR}/*.sam
do fname=$(basename ${file})
# move files to SAM_OUTPUT
mv $file ${SAM_OUTPUT}/${fname/clipped_}
done

printf "Finished transfering files.\n"
