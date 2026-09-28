process MARK_DUPLICATES_DNA {
    tag "$meta.id"
    label 'process_high'
    // One worker per reference sequence, each holding ~165 bytes per primary read of
    // its contig: twelve large contigs at full depth (~850M primaries) need ~80 GB.
    label 'process_high_memory'

    // Same environment as SPLIT_BAM: the script needs pysam and numpy only, and numpy
    // comes in with biopython here, so no new container is pulled.
    conda "conda-forge::editdistance=0.6.0 bioconda::pysam=0.19.1 conda-forge::pygtrie=2.5.0 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea:ad24cd6a9acfe3ff51beb4c454076e18e778f7c0-0' :
        'biocontainers/mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea:ad24cd6a9acfe3ff51beb4c454076e18e778f7c0-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.bam")        , emit: bam
    tuple val(meta), path("*.metrics.txt")        , emit: metrics
    tuple val(meta), path("*.dedup_summary.tsv")  , emit: summary
    path "versions.yml"                           , emit: versions_mark_duplicates_dna, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.dedup"
    if ("${bam}" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    mark_dna_duplicates.py \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --metrics ${prefix}.metrics.txt \\
        --summary ${prefix}.dedup_summary.tsv \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.dedup"
    """
    touch ${prefix}.bam
    touch ${prefix}.metrics.txt
    touch ${prefix}.dedup_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
    END_VERSIONS
    """
}
