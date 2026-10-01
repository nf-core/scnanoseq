#!/usr/bin/env bash

# Sum the per-chunk umi_tools dedup logs of one sample into a single
# duplication summary row.
#
# Deduplication runs once per chromosome / transcript-group chunk (and, with
# --dedup_per_gene, once more per chunk for the positional fallback), so no
# single umi_tools log carries the sample-level counts. This adds the
# "Input Reads" and "Number of reads out" figures across every log it is given.

set -euo pipefail

output=""
sample=""
logs=()

while [[ $# -gt 0 ]]
do
    case "$1" in
        --output) output=$2; shift;;
        --sample) sample=$2; shift;;
        --*) echo "Unknown option $1" >&2 && exit 1;;
        *) logs+=("$1");;
    esac
    shift
done

if [[ -z "$output" || -z "$sample" || ${#logs[@]} -eq 0 ]]
then
    echo "usage: $(basename "$0") --sample <name> --output <tsv> <dedup.log> [<dedup.log> ...]" >&2
    exit 1
fi

awk -v sample="$sample" -v OFS='\t' '
    # umi_tools writes the read tallies in "most common first" order on a single
    # line, so "Input Reads" is not necessarily the first field. Match it anywhere.
    match($0, /Input Reads: [0-9]+/) {
        split(substr($0, RSTART, RLENGTH), field, ": ")
        input_reads += field[2]
    }
    match($0, /Number of reads out: [0-9]+/) {
        split(substr($0, RSTART, RLENGTH), field, ": ")
        output_reads += field[2]
        n_runs++
    }
    END {
        if (n_runs == 0) {
            print "No umi_tools dedup read counts found in the supplied logs" > "/dev/stderr"
            exit 1
        }
        duplicate_reads = input_reads - output_reads
        duplication_rate = (input_reads > 0) ? duplicate_reads / input_reads : 0
        print "sample", "input_reads", "output_reads", "duplicate_reads", "duplication_rate", "n_dedup_runs"
        print sample, input_reads, output_reads, duplicate_reads, sprintf("%.4f", duplication_rate), n_runs
    }
' "${logs[@]}" > "$output"
