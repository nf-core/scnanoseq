#!/usr/bin/env bash

set -e

########################
### INPUT PARAMETERS ###
########################

IN_FASTA=""
IN_GTF=""
DELIMITER=""

while [[ $# -gt 0 ]]; do
    flag=$1
    case "${flag}" in
        -f) IN_FASTA=$2; shift;;
        -g) IN_GTF=$2; shift;;
        -d) DELIMITER=$2; shift;;
    esac
    shift
done


#####################
### MERGE REGIONS ###
#####################


# Shrink down the GTF so the greps aren't as expensive

SIMPLE_GTF="simplified.gtf"
grep transcript_id $IN_GTF | awk '{print $1"\t"$12}' | sed 's#[";]##g' | uniq > $SIMPLE_GTF


# Iterate over IN_FASTA
grep '^>' $IN_FASTA | sed 's#>##g' | while read feature_name
do
    # Try and determine delimiter if not provided
    if [ -z $DELIMITER ]; then
        max_delim_count=0
        max_delim_str=""

        for test_delim in "|" ";"
        do
            delim_counts=$(echo $feature_name | grep -o -n $test_delim | wc -l)

            if [[ "$delim_counts" -gt "$max_delim_count" ]]; then
                max_delim_count=$delim_counts
                max_delim_str=$test_delim
            fi

        done

        DELIMITER=$max_delim_str

    fi

    # Determine which chr the transcript originates from
    # NOTE: This does not need to be correct, its just for being able to split downstream
    echo $feature_name

    # Its possible there wasn't a delim found, in which case it might be space
    transcript_name=""
    if [ -z $DELIMITER ]; then
        transcript_name=$(echo $feature_name | cut -f1 -d " ")
    else
        transcript_name=$(echo $feature_name | cut -f1 -d $DELIMITER)
    fi
    echo $transcript_name

    chr=$(grep $transcript_name $SIMPLE_GTF | awk '{print $1}' | sort -u | head -n1 | cut -f2 -d':')
    echo $chr
    if [ -z "$chr" ]; then
        chr="chr_Unknown"
    fi
    echo ""

    echo $feature_name >> ${chr}.transcripts.txt
done
