## Yunfei Dai
## Feb 05 2021

# using sam_to_wig.py created by Amy Tang (original Apr 4 2020)

DIR=${1?Error: "enter directory of SAM files"}
REF=${2?Error: "enter directory of genome annotation file (.gbk)"}
OUTPUT="$(dirname $DIR)/wig_files"

mkdir -p ${OUTPUT}

parallel python ~/GeisingerLab/TnSeq_Processing/src/sam_to_wig.py -g ${REF} -i {} -o {.}.wig ::: ${DIR}/*.sam

printf "Finished converting SAM files to WIG files.\n"

for file in ${DIR}/*.wig
do fname=(basename ${file})
mv $file ${OUTPUT}/${fname}
done

