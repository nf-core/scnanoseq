process ISOQUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::isoquant=3.12.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.12.2--pyh106432d_0' :
        'biocontainers/isoquant:3.12.2--pyh106432d_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(gtf)
    val group_category

    output:
    tuple val(meta), path("*/*/*.gene_grouped*counts.matrix.mtx"),       emit: grouped_gene_mtx
    tuple val(meta), path("*/*/*.gene_grouped*counts.barcodes.tsv"),            emit: grouped_gene_mtx_barcodes
    tuple val(meta), path("*/*/*.gene_grouped*counts.features.tsv"),            emit: grouped_gene_mtx_features
    tuple val(meta), path("*/*/*.transcript_grouped*counts.matrix.mtx"), emit: grouped_transcript_mtx
    tuple val(meta), path("*/*/*.transcript_grouped*counts.barcodes.tsv"),      emit: grouped_transcript_mtx_barcodes
    tuple val(meta), path("*/*/*.transcript_grouped*counts.features.tsv"),      emit: grouped_transcript_mtx_features
    path "versions.yml",                                                emit: versions

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
        isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dir_prefix=\$(isoquant_output/${prefix})
    mkdir -p \$dir_prefix


    touch \$dir_prefix/${prefix}.read_assignments.tsv.gz
    touch \$dir_prefix/${prefix}.corrected_reads.bed.gz
    touch \$dir_prefix/${prefix}.transcript_tpm.tsv
    touch \$dir_prefix/${prefix}.transcript_counts.tsv
    touch \$dir_prefix/${prefix}.gene_tpm.tsv
    touch \$dir_prefix/${prefix}.gene_counts.tsv
    touch \$dir_prefix/isoquant.log
    touch \$dir_prefix/${prefix}.exon_counts.tsv
    touch \$dir_prefix/${prefix}.intron_counts.tsv
    touch \$dir_prefix/${prefix}.novel_vs_known.SQANTI-like.tsv
    touch \$dir_prefix/${prefix}.gene_grouped_tpm.tsv
    touch \$dir_prefix/${prefix}.gene_grouped_counts.tsv
    touch \$dir_prefix/${prefix}.transcript_grouped_tpm.tsv
    touch \$dir_prefix/${prefix}.transcript_grouped_counts.tsv
    touch \$dir_prefix/${prefix}.exon_grouped_counts.tsv
    touch \$dir_prefix/${prefix}.intron_grouped_counts.tsv
    touch \$dir_prefix/${prefix}.transcript_models.gtf
    touch \$dir_prefix/${prefix}.transcript_model_reads.tsv.gz
    touch \$dir_prefix/${prefix}.transcript_model_tpm.tsv
    touch \$dir_prefix/${prefix}.transcript_model_counts.tsv
    touch \$dir_prefix/${prefix}.extended_annotation.gtf
    touch \$dir_prefix/${prefix}.transcript_model_grouped_counts.tsv
    touch \$dir_prefix/${prefix}.transcript_model_grouped_tpm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
    END_VERSIONS
    """
}
