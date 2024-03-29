process BLAZE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda "conda-forge::python=3.7 conda-forge::biopython conda-forge::pandas conda-forge::numpy conda-forge::tqdm conda-forge::matplotlib conda-forge::pip conda-forge::python-levenshtein"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d7d88ab042b250d0c097b70dc83ce3eb755b9be3:39cb062b9290449cd9de547d494e563502b5e0eb-0' :
        'biocontainers/mulled-v2-d7d88ab042b250d0c097b70dc83ce3eb755b9be3:39cb062b9290449cd9de547d494e563502b5e0eb-0' }"

    input:
    tuple val(meta), path(reads)
    val in_whitelist

    output:
    tuple val(meta), path("*.putative_bc.csv") , emit: putative_bc
    tuple val(meta), path("*.whitelist.csv")   , emit: whitelist
    tuple val(meta), path("*.bc_count.txt")    , emit: bc_count
    tuple val(meta), path("*.knee_plot.png")   , emit: knee_plot
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping BLAZE code
    def cell_count = "${meta.cell_count}"

    """
    blaze.py \\
        --expect-cells=${cell_count} \\
        --full-bc-whitelist=${in_whitelist} \\
        --out-putative-bc=${prefix}.putative_bc \\
        --out-bc-whitelist=${prefix}.whitelist \\
        --threads=$task.cpus \\
        ${args} \\
        \$(pwd)

    sed -i 's#-1##g' ${prefix}.whitelist.csv
    grep -f ${prefix}.whitelist.csv ${prefix}.putative_bc.csv | cut -f2 -d',' | sort | uniq -c | awk '{print \$2","\$1}'> ${prefix}.bc_count.txt
    mv knee_plot.png ${prefix}.knee_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blaze: $VERSION
    END_VERSIONS
    """
}
