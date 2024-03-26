#!/usr/bin/env bash

########################
### INPUT PARAMETERS ###
########################

IN_DIR=`pwd`
OUT_FILE=""
FEATURE_TYPE=""

while [[ $# -gt 0 ]]; do
    flag=$1

    case "${flag}" in
        -i) IN_DIR=$2; shift;;
        -o) OUT_FILE=$2; shift;;
        -f) FEATURE_TYPE=$2; shift;;
        *) echo "Unknown option $1 ${reset}" && exit 1
    esac

    shift
done

###############
### COMBINE ###
###############

# Print out the header
echo -e "Sample\tFeature_Type\tEstimated_Cell_Number\tMean_Reads_Per_Cell\tMedian_Features_Per_Cell\tTotal_Number_of_Features" > .header.tmp

# The files are just input from the command line
for in_file in ${IN_DIR}/*.csv
do
    # Parse the values from the individual files
    # Each file only has one sample, so the parsing is easier

    sample=$(basename $in_file | cut -f1 -d'.')
    values=$(tail -n1 $in_file | tr ',' '\t')

    # Output out the values
    echo -e "${sample}\t${FEATURE_TYPE}\t${values}" >> ${FEATURE_TYPE}.data.tmp
done

# cat the header file with the data
cat .header.tmp ${FEATURE_TYPE}.data.tmp > ${OUT_FILE}

# Remove intermediary files
rm -f *.tmp
