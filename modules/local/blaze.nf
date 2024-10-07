process BLAZE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda "atrull314::fast_edit_distance=1.2.1 conda-forge::matplotlib=3.8.4 conda-forge::biopython=1.83 conda-forge::pandas=2.2.2 conda-forge::numpy=2.0.0rc2 conda-forge::tqdm=4.66.4"

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://agtrull314/blaze:2.2.0' :
        'docker.io/agtrull314/blaze:2.2.0'}"

    input:
    tuple val(meta), path(reads)
    path in_whitelist

    output:
    tuple val(meta), path("*.putative_bc.no_header.csv") , emit: putative_bc
    tuple val(meta), path("*.whitelist.csv")             , emit: whitelist
    tuple val(meta), path("*.bc_count.txt")              , emit: bc_count
    tuple val(meta), path("*.knee_plot.png")             , emit: knee_plot
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.2.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping BLAZE code
    def cell_count = "${meta.cell_count}"

    """
    main.py \\
        --expect-cells ${cell_count} \\
        --full-bc-whitelist ${in_whitelist} \\
        --output-prefix ${prefix}. \\
        --threads $task.cpus \\
        ${args} \\
        \$(pwd)


    tail -n +2 ${prefix}.putative_bc.csv > ${prefix}.putative_bc.no_header.csv

    cat ${prefix}.putative_bc.no_header.csv | \\
        cut -f2 -d',' | \\
        sort -T \$(pwd) | \\
        uniq -c | \\
        awk '{print \$2","\$1}'> ${prefix}.bc_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blaze: $VERSION
    END_VERSIONS
    """
}
