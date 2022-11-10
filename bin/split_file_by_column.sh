#!/usr/bin/env bash

########################
### INPUT_PARAMETERS ###
########################

INFILE=""
SPLIT_AMT=""

while [[ $# -gt 0 ]]
do
    flag=$1

    case "${flag}" in
        -i) INFILE=$2; shift;;
        -s) SPLIT_AMT=$2; shift;;
        *) echo "Unknown option $1" && exit 1
    esac

    shift
done

######################
### SPLIT THE FILE ###
######################

# The first column will need to be added to all files, so lets isolate that one
dir_name=$(dirname $INFILE)
cut -d $'\t' -f1 $INFILE > ${INFILE}.first_col.tmp
cut -d $'\t' -f2- $INFILE > ${INFILE}.other_cols.tmp

# Determine the number of columns
num_col=$(awk 'NR==1{print NF}' "${INFILE}.other_cols.tmp")

# We will use this variable to create a window
for i in $(seq 1 $SPLIT_AMT $num_col)
do
    start_val=$i
    end_val=$(($i+$SPLIT_AMT-1))
    
    cut -f$start_val-$end_val -d $'\t' ${INFILE}.other_cols.tmp > ${INFILE}.other_cols.${i}.tmp

    # Check to make sure the file has data and is not just an empty space
    if [ -s ${INFILE}.other_cols.${i}.tmp ]
    then
        paste ${INFILE}.first_col.tmp ${INFILE}.other_cols.${i}.tmp -d $"\t" > $(basename $INFILE | sed 's#.tsv##g').${i}.tsv
    fi
done

# Clean up the temp files
#rm -rf ${INFILE}.*.tmp
