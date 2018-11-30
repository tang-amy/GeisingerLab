#!/home/bin/bash

# user input should give the directory name of the zipped fastq files
# script will create directories as necessary in parent directory for organization of files


# assumed var
DIR=${PWD}

# getting proper parameters via flags by user
while getopts ":i:g:r:h" opt; do
case ${opt} in
  i)
    DIR="${OPTARG}" >&2
    ;;
  g)
    GENBANK="${OPTARG}" >&2
    ;;
  r)
    REF="${OPTARG}" >&2
    ;;
  h)
    echo "usage: tnseq-process.sh -i [DIR] -g [GENBANK FILE (.gbk)] -r [BOWTIE REFERENCES]"
    echo "If -i flag not used, the directory containing the gunzipped files will be assumed to be the current directory."
    echo "The -g is mandatory and the proper genbank (.gbk) file must be provided for the species."
    echo "The -r flag is optional. If reference is not provided, you will be asked to supply FASTA (.fasta) file for creating the references."
    ;;
  \?)
    echo "Invalid flag: -${OPTARG}" >&2
    exit 1
esac
done

# ensuring genbank file is entered
if [ ! "${GENBANK}" ]
then
    echo "Missing genbank file."
    exit 1
fi

# asking for fasta file to create reference files later on if no REF folder for bowtie is provided
if [ ! "${REF}" ]
then
    while [ ! "${FASTA}" ]; do
        echo "Please provide a FASTA file to create the bowtie reference (.ebwt) files from and then [ENTER]: "
        read FASTA
    done
    echo "Please name the reference the FASTA file is used to create and then [ENTER]: "
    read REFERENCE
    echo "Type in location to save the bowtie reference and then [ENTER] or [ENTER] to save in parent directory: "
    read LOC

    # if missing name for reference, name it according to fasta file
    if [ ! "${REFERENCE}" ]
    then
    REFERENCE="${FASTA/%.fasta}"
    fi
    # if no location specified, then create the reference folder in parent directory
    if [ ! "${LOC}" ]
    then
    LOC="$(dirname ${DIR})/${REFERENCE}_ebwt/"
    fi

    # set the reference properly
    REF="${LOC}/${REFERENCE}"

    # build bow-tie reference file
    bowtie-build "${FASTA}" "${REF}"
fi



#make prot table if it does not exist
echo "Does a prot_table need to be created with this genbank file? Type Y/N and then [ENTER]: "
read answer
if [ "${answer}" == "Y" ]
then
    echo "Type in location and name to save the prot_table (include .prot_table) to followed by [ENTER]: "
    read PROT_TABLE
    python make_prot_table.py -i ${GENBANK} -o ${PROT_TABLE}
    echo "Finished creating prot_table for Transit analysis."
fi

# getting name of parent dir
PDIR=$(dirname ${DIR})

# unzipping and clipping fastq.gz files
bash process_zipped_fastq.sh "${DIR}"

# running bowtie on clipped files
bash run_bowtie_batch.sh "${PDIR}/clipped_fastq_files" "$REF"

# python requires absolute paths to be given, force genbank's path into an absolute path
GENBANK=$(realpath $GENBANK)
PDIR=$(realpath $PDIR)

# create the wig files
python batch_map2wig.py -i "${PDIR}/map_files/" -o "${PDIR}/wig_files/" -g "$GENBANK"

#mkdir -p "${PDIR}/wig_files"
#for file in ${PDIR}/map_files/*.map
#do
#    python bowtie_to_wig.py -g "$GENBANK" -i "$file" -o "${PDIR}/wig_files/$(basename ${file/%map/wig})"
#done
#echo "Finished processing map files into wig files."


# create the merged wig file
python merge_wigs.py -i "${PDIR}/wig_files/" -o "${PDIR}/merged_wig.wig" -g "$GENBANK"
echo "Created a merged .wig file of all reads."

echo "Done."