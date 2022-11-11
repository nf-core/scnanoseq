#!/usr/bin/env bash

########################
### INPUT PARAMETERS ###
########################

IN_DIR=""
OUT_FILE=""

while [[ $# -gt 0 ]]; do
    flag=$1

    case "${flag}" in
        -i) IN_DIR=$2; shift;;
        -o) OUT_FILE=$2; shift;;
        *) echo "Unknown option $1 ${reset}" && exit 1
    esac

    shift
done

OUT_DIR=$(dirname $OUT_FILE)

# Print out the header
echo -e "Sample\tFeature_Type\tEstimated_Cell_Number\tMean_Reads_Per_Cell\tMedian_Genes_Per_Cell" > ${OUT_FILE}.header.tmp

# The files are just input from the command line
#for in_file in $@
for in_file in ${IN_DIR}/*.tsv
do
    # Parse the values from the individual files
    # Each file only has one sample, so the parsing is easier

    sample=$(basename $in_file | cut -f1 -d'.')
    feature_type=$(basename $(dirname $in_file))
    values=$(tail -n1 $in_file | tr ',' '\t')

    # Output out the values
    echo -e "${sample}\t${feature_type}\t${values}" >> ${OUT_DIR}/${feature_type}.data.tmp
done

# cat the header file with the data
for data_file in ${OUT_DIR}/*.data.tmp
do
    feature=$(basename $data_file | cut -f1 -d'.')
    cat ${OUT_FILE}.header.tmp ${data_file} > ${OUT_FILE}
done

# Remove intermediary files
rm -f ${OUT_DIR}/*.tmp
