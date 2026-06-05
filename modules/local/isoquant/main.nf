process ISOQUANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::isoquant=3.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.13.0--pyh106432d_0' :
        'biocontainers/isoquant:3.13.0--pyh106432d_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(gtf)
    val group_category

    output:
    tuple val(meta), path("*/*/*.gene_grouped*counts.matrix.mtx"),         emit: grouped_gene_mtx
    tuple val(meta), path("*/*/*.gene_grouped*counts.barcodes.tsv"),       emit: grouped_gene_mtx_barcodes
    tuple val(meta), path("*/*/*.gene_grouped*counts.features.tsv"),       emit: grouped_gene_mtx_features
    tuple val(meta), path("*/*/*.transcript_grouped*counts.matrix.mtx"),   emit: grouped_transcript_mtx
    tuple val(meta), path("*/*/*.transcript_grouped*counts.barcodes.tsv"), emit: grouped_transcript_mtx_barcodes
    tuple val(meta), path("*/*/*.transcript_grouped*counts.features.tsv"), emit: grouped_transcript_mtx_features
    path "versions.yml", emit: versions_isoquant, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def group_flag = group_category ? "--read_group $group_category" : ""
    def ref_flag   = fasta ? "--reference $fasta" : ""
    def gtf_flag   = gtf ? "--genedb $gtf": ""

    """
    export HOME=\$(pwd)

    isoquant ${args} \\
        --threads $task.cpus \\
        --prefix $prefix \\
        --bam ${bam} \\
        --output ${prefix} \\
        ${ref_flag} \\
        ${gtf_flag} \\
        ${group_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant -v | sed 's#IsoQuant ##')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p isoquant_out/tmp

    touch isoquant_out/tmp/${prefix}.gene_grouped.counts.matrix.mtx
    touch isoquant_out/tmp/${prefix}.gene_grouped.counts.barcodes.tsv
    touch isoquant_out/tmp/${prefix}.gene_grouped.counts.features.tsv
    touch isoquant_out/tmp/${prefix}.transcript_grouped.counts.matrix.mtx
    touch isoquant_out/tmp/${prefix}.transcript_grouped.counts.barcodes.tsv
    touch isoquant_out/tmp/${prefix}.transcript_grouped.counts.features.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant -v | sed 's#IsoQuant ##')
    END_VERSIONS
    """
}
