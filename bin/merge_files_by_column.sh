#!/usr/bin/env bash

########################
### INPUT PARAMETERS ###
########################
OUT_FILE=$1
shift
IN_FILES=$@

###################
### MERGE FILES ###
###################

in_files=$IN_DIR/*
local_tmp_dir=$(dirname $OUT_FILE)/tmp_dir

mkdir -p $local_tmp_dir

for file in $IN_FILES
do
    echo "Evaluating $file"
    file_prefix=$local_tmp_dir/$(basename $file)

    # Lets remove the header so we can sort
    echo "Removing header..."
    head -n1 $file > ${file_prefix}.header.tmp

    # Sort and save off the file
    echo "File sorting..."
    tail -n +2 $file | sort -k1 > ${file_prefix}.other_rows.tmp

    # Remove and save off the first column (Its repetitive)
    echo "Removing first column..."
    cut -f1 -d$'\t' $file > ${file_prefix}.first_col.tmp
    cut -f1 -d$'\t' ${file_prefix}.header.tmp > ${file_prefix}.header.first_col.tmp
    
    cut -f2- -d$'\t' ${file_prefix}.other_rows.tmp > ${file_prefix}.no_first_col.tmp
    cut -f2- -d$'\t' ${file_prefix}.header.tmp > ${file_prefix}.header.no_first_col.tmp

    # Now add the header and file back together
    echo "Adding the header and file back together..."
    cat ${file_prefix}.header.no_first_col.tmp ${file_prefix}.no_first_col.tmp > ${file_prefix}.sorted.tmp

    # Now lets add the stuff back to the growing output file
    echo "Growing the file..."
    if [[ -f ${OUT_FILE}.no_first_col.tmp ]]; then
        paste ${OUT_FILE}.no_first_col.tmp ${file_prefix}.sorted.tmp > ${OUT_FILE}.no_first_col.tmp2
        mv ${OUT_FILE}.no_first_col.tmp2 ${OUT_FILE}.no_first_col.tmp
    else
        mv ${file_prefix}.sorted.tmp ${OUT_FILE}.no_first_col.tmp
    fi

    # Lets also make sure gene columns match to be extra safe
    echo "Checking gene column"
    if [[ ! -f ${file_prefix}.first_col ]]
    then
        cp ${file_prefix}.first_col.tmp ${file_prefix}.first_col
    else
        cmp --s ${file_prefix}.first_col.tmp ${file_prefix}.first_col || echo "First columns do not match!"
    fi
done

echo "After the for loop"

paste ${file_prefix}.first_col ${OUT_FILE}.no_first_col.tmp > ${OUT_FILE}
