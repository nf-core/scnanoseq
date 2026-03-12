get_fastqc_counts()
{
    fastqc_file=$1
    counts=$(unzip -p "${fastqc_file}" "$(basename "${fastqc_file}" .zip)/fastqc_data.txt" | \
        grep 'Total Sequences' | \
        cut -f2 -d$'\t')
    echo "$counts"
}

get_nanoplot_counts()
{
    nanoplot_file=$1
    counts=$(grep 'Number of reads' "$nanoplot_file" | awk '{print $NF}' | cut -f1 -d'.' | sed 's/,//g')
    echo "$counts"
}

output=""
input=""

while [[ $# -gt 0 ]]
do
    flag=$1
    case "${flag}" in
        --input) input=$2; shift;;
        --output) output=$2; shift;;
        *) echo "Unknown option $1" && exit 1
    esac
    shift
done

header="sample,base_fastq_counts,trimmed_read_counts,extracted_read_counts,corrected_read_counts"
echo "$header" > "$output"

# Collect all sample names from both barcode file types
sample_names=$(find "$input" -type f -name "*.corrected_bc_umi.tsv" -o -name "*_known_barcodes.txt" | \
    sed -E 's|.*/||' | sed -E 's/_known_barcodes\.txt$//; s/\.corrected_bc_umi\.tsv$//' | sort -u)

for sample_name in $sample_names
do
    raw_fastqc="${sample_name}.raw_fastqc.zip"
    raw_nanoplot="${sample_name}.raw_NanoStats.txt"

    trim_fastqc="${sample_name}.trimmed_fastqc.zip"
    trim_nanoplot="${sample_name}.trimmed_NanoStats.txt"

    extract_fastqc="${sample_name}.extracted_fastqc.zip"
    extract_nanoplot="${sample_name}.extracted_NanoStats.txt"

    corrected_tsv="${sample_name}.corrected_bc_umi.tsv"
    known_barcodes="${sample_name}_known_barcodes.txt"

    data="$(basename "$sample_name")"

    ####################
    # RAW FASTQ COUNTS #
    ####################
    if [[ -s "$raw_fastqc" ]]; then
        fastqc_counts=$(get_fastqc_counts "$raw_fastqc")
        data="$data,$fastqc_counts"
    elif [[ -s "$raw_nanoplot" ]]; then
        nanoplot_counts=$(get_nanoplot_counts "$raw_nanoplot")
        data="$data,$nanoplot_counts"
    else
        data="$data,"
    fi

    ###############
    # TRIM COUNTS #
    ###############
    if [[ -s "$trim_fastqc" ]]; then
        trim_counts=$(get_fastqc_counts "$trim_fastqc")
        data="$data,$trim_counts"
    elif [[ -s "$trim_nanoplot" ]]; then
        nanoplot_counts=$(get_nanoplot_counts "$trim_nanoplot")
        data="$data,$nanoplot_counts"
    else
        data="$data,"
    fi

    #####################
    # PREEXTRACT COUNTS #
    #####################
    if [[ -s "$extract_fastqc" ]]; then
        extract_counts=$(get_fastqc_counts "$extract_fastqc")
        data="$data,$extract_counts"
    elif [[ -s "$extract_nanoplot" ]]; then
        nanoplot_counts=$(get_nanoplot_counts "$extract_nanoplot")
        data="$data,$nanoplot_counts"
    else
        data="$data,"
    fi

    ##################
    # CORRECT COUNTS #
    ##################
    if [[ -s "$known_barcodes" ]]; then
        correct_sum=$(awk -F'\t' '{if ($2 != "") sum += $2} END {print sum}' "$known_barcodes")
        data="$data,$correct_sum"
    elif [[ -s "$corrected_tsv" ]]; then
        correct_counts=$(cut -f6 "$corrected_tsv" | awk '{if ($0 != "") print $0}' | wc -l)
        data="$data,$correct_counts"
    else
        data="$data,"
    fi

    echo "$data" >> "$output"
done
