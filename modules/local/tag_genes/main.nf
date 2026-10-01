process TAG_GENES {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::pysam=0.19.1 conda-forge::python=3.10.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.19.1--py310hff46b53_1' :
        'biocontainers/pysam:0.19.1--py310hff46b53_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(gtf)

    output:
    tuple val(meta), path("*.gene_tagged.bam"), emit: gene_tagged_bam
    tuple val(meta), path("*.gene_assignment.tsv"), emit: summary
    path "versions.yml"                       , emit: versions_tag_genes, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    tag_genes.py \\
        -b ${bam} \\
        -g ${gtf} \\
        -o ${prefix}.gene_tagged.bam \\
        -s ${prefix}.gene_assignment.tsv \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene_tagged.bam
    touch ${prefix}.gene_assignment.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
