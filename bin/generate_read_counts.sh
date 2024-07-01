
get_fastqc_counts()
{
    fastqc_file=$1
    counts=$(unzip -p ${fastqc_file} $(basename ${fastqc_file} .zip)/fastqc_data.txt | \
        grep 'Total Sequences' | \
        cut -f2 -d$'\t')
    echo $counts

}

output=""
input=""

while [[ $# -gt 0 ]]
do
    flag=$1

    case "${flag}" in
        --input) input=$2; shift;;
        --output) output=$2; shift;;
        *) echo "Unknown option $1 ${reset}" && exit 1
    esac
    shift
done

header=""
data=""

header="sample,base_fastq_counts,trimmed_read_counts,extracted_read_counts,corrected_read_counts"
echo "$header" > $output

for sample_name in $(for file in $(readlink -f $input)/*.zip; do echo $file; done | cut -f1 -d'.' | sort -u)
do
    raw_fastqc="${sample_name}.raw_fastqc.zip"
    trim_fastqc="${sample_name}.trimmed_fastqc.zip"
    extract_fastqc="${sample_name}.extracted_fastqc.zip"
    correct_csv="${sample_name}.corrected_bc_umi.tsv"
    data="$(basename $sample_name)"

    # RAW FASTQ COUNTS
    
    if [[ -s "$raw_fastqc" ]]
    then
        fastqc_counts=$(get_fastqc_counts "$raw_fastqc")
        data="$data,$fastqc_counts"
    else
        data="$data,"
    fi
    
    # TRIM COUNTS
    
    if [[ -s "$trim_fastqc" ]]
    then
        trim_counts=$(get_fastqc_counts "$trim_fastqc")
        data="$data,$trim_counts"
    else
        data="$data,"
    fi
    
    # PREEXTRACT COUNTS
    
    if [ -s "$extract_fastqc" ]
    then
        extract_counts=$(get_fastqc_counts "$extract_fastqc")
        data="$data,$extract_counts"
    else
        data="$data,"
    fi
    
    # CORRECT COUNTS
    
    
    if [ -s $correct_csv ]
    then
        correct_counts=$(cut -f6 $correct_csv | awk '{if ($0 != "") {print $0}}' | wc -l)
        data="$data,$correct_counts"
    else
        data="$data,"
    fi
    echo "$data" >> $output
done
