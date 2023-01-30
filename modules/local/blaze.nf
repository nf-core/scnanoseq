process BLAZE {
    tag "$meta.id"
    label 'process_low'

    conda ("conda-forge::python=3.7 conda-forge::biopython conda-forge::pandas conda-forge::numpy conda-forge::tqdm conda-forge::matplotlib conda-forge::pip conda-forge::python-levenshtein")

    input:
    tuple val(meta), path(reads)
    val exp_cell_amount
    val in_whitelist

    output:
    tuple val(meta), path("*.putative_bc.csv") , emit: putative_bc 
    tuple val(meta), path("*.whitelist.csv")   , emit: whitelist
    tuple val(meta), path("*.bc_count.txt")    , emit: bc_count
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    blaze.py \\
        --expect-cells=${exp_cell_amount} \\
        --threads=12 \\
        --full-bc-whitelist=${in_whitelist} \\
        --out-putative-bc=${prefix}.putative_bc \\
        --out-bc-whitelist=${prefix}.whitelist \\
       \$(pwd) 

    sed -i 's#-1##g' ${prefix}.whitelist.csv
    cut -f2 -d',' ${prefix}.putative_bc.csv | sort | uniq -c | awk '{if (\$2 != "") {print \$2"\\t\\t"\$1"\\t"}}' > ${prefix}.bc_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blaze: 1.0 
    END_VERSIONS
    """
}
