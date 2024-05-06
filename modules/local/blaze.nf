process BLAZE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda "conda-forge::python=3.7 conda-forge::biopython conda-forge::pandas conda-forge::numpy conda-forge::tqdm conda-forge::matplotlib conda-forge::pip conda-forge::python-levenshtein"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/biocontainers/blaze:2.2.0_cv1' :
    //    'biocontainers/blaze:2.2.0_cv1' }"

    container "docker://biocontainers/blaze:2.2.0_cv1"

    input:
    tuple val(meta), path(reads)
    val in_whitelist

    output:
    tuple val(meta), path("*.putative_bc.mod.csv") , emit: putative_bc
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
    blaze \\
        --expect-cells ${cell_count} \\
        --full-bc-whitelist ${in_whitelist} \\
        --output-prefix ${prefix}. \\
        --threads $task.cpus \\
        ${args} \\
        \$(pwd)


    tail -n +2 ${prefix}.putative_bc.csv > ${prefix}.putative_bc.mod.csv

    cat ${prefix}.putative_bc.mod.csv | cut -f2 -d',' | sort -T \$(pwd) | uniq -c | awk '{print \$2","\$1}'> ${prefix}.bc_count.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blaze: $VERSION
    END_VERSIONS
    """
}
