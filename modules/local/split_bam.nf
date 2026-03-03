process SPLIT_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::editdistance=0.6.0 bioconda::pysam=0.19.1 conda-forge::pygtrie=2.5.0 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea:ad24cd6a9acfe3ff51beb4c454076e18e778f7c0-0' :
        'biocontainers/mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea:ad24cd6a9acfe3ff51beb4c454076e18e778f7c0-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(contig_file)

    output:
    tuple val(meta), path("*.split.bam"), emit: split_bam
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    split_bam.py \\
        ${args} \\
        --input ${bam} \\
        --output ${prefix}.split.bam \\
        --contig_file ${contig_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.corrected_bc_umi.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
