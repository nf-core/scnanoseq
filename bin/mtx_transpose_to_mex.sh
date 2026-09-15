#!/usr/bin/env bash

# Rewrite a cells x features MatrixMarket matrix as a features x cells MEX
# directory: matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz.
#
# kallisto quant-tcc, bustools count and oarfish all write cells as rows and
# features as columns; Seurat's Read10X and CellBender want the transpose. For
# a coordinate MatrixMarket file that is only a matter of swapping the first
# two fields of the dimension line and of every entry, which is done in one
# streamed pass straight into gzip.
#
# With --integer the header's field type is rewritten from real to integer,
# and every value is first checked to be a plain non-negative integer, so the
# header never promises what the values do not deliver. bustools prints its
# counts as doubles, so a fractional value (from --multimapping) or a %g
# exponent (1e+06) would otherwise be silently misread by scipy's mmread.

set -euo pipefail

mtx=""
features=""
barcodes=""
outdir=""
integer=0

while [[ $# -gt 0 ]]
do
    case "$1" in
        --mtx) mtx=$2; shift;;
        --features) features=$2; shift;;
        --barcodes) barcodes=$2; shift;;
        --outdir) outdir=$2; shift;;
        --integer) integer=1;;
        *) echo "Unknown option $1" >&2 && exit 1;;
    esac
    shift
done

if [[ -z "$mtx" || -z "$features" || -z "$barcodes" || -z "$outdir" ]]
then
    echo "usage: $(basename "$0") --mtx <cells_x_features.mtx> --features <tsv> --barcodes <tsv> --outdir <dir> [--integer]" >&2
    exit 1
fi

mkdir -p "$outdir"

awk -v integer="$integer" '
    /^%/ {
        if (integer) sub(/ real /, " integer ")
        print
        next
    }
    NF == 0 { next }
    NF != 3 {
        print FILENAME ": line " FNR " does not have three fields: " $0 > "/dev/stderr"
        exit 1
    }
    integer && $3 !~ /^[0-9]+$/ {
        print FILENAME ": line " FNR " is not a whole number, refusing to declare the matrix integer: " $0 > "/dev/stderr"
        exit 1
    }
    { print $2, $1, $3 }
' "$mtx" | gzip -c > "$outdir/matrix.mtx.gz"

gzip -c "$features" > "$outdir/features.tsv.gz"
gzip -c "$barcodes" > "$outdir/barcodes.tsv.gz"
