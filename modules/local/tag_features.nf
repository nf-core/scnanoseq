process TAG_FEATURES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pysam=0.19.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.19.1--py310hff46b53_1' :
        'quay.io/biocontainers/pysam:0.19.1--py310hff46b53_1' }"

    input:
    tuple val(meta), path(bam), path(feature_file)

    output:
    tuple val(meta), path("*.feature.bam"), emit: feature_bam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    tag_features.py \\
        -i $bam \\
        -f $feature_file \\
        -o ${prefix}.feature.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
