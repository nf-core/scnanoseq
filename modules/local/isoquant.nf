process ISOQUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::isoquant=3.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.6.1--hdfd78af_0' :
        'biocontainers/isoquant:3.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(gtf)
    val group_category

    output:
    tuple val(meta), path("*/*/*.read_assignments.tsv.gz"),             emit: read_assignments
    tuple val(meta), path("*/*/*.corrected_reads.bed.gz"),              emit: corrected_reads
    tuple val(meta), path("*/*/*.transcript_tpm.tsv"),                  emit: transcript_tpm
    tuple val(meta), path("*/*/*.transcript_counts.tsv"),               emit: transcript_counts
    tuple val(meta), path("*/*/*.gene_tpm.tsv"),                        emit: gene_tpm
    tuple val(meta), path("*/*/*.gene_counts.tsv"),                     emit: gene_counts
    tuple val(meta), path("*/isoquant.log"),                            emit: log
    tuple val(meta), path("*/*/*.exon_counts.tsv"),                     emit: exon_counts,                     optional: true
    tuple val(meta), path("*/*/*.intron_counts.tsv"),                   emit: intron_counts,                   optional: true
    tuple val(meta), path("*/*/*.novel_vs_known.SQANTI-like.tsv"),      emit: sqanti_output,                   optional: true
    tuple val(meta), path("*/*/*.gene_grouped_tpm.tsv"),                emit: grouped_gene_tpm,                optional: true
    tuple val(meta), path("*/*/*.gene_grouped_counts.tsv"),             emit: grouped_gene_counts,             optional: true
    tuple val(meta), path("*/*/*.transcript_grouped_tpm.tsv"),          emit: grouped_transcript_tpm,          optional: true
    tuple val(meta), path("*/*/*.transcript_grouped_counts.tsv"),       emit: grouped_transcript_counts,       optional: true
    tuple val(meta), path("*/*/*.exon_grouped_counts.tsv"),             emit: grouped_exon_counts,             optional: true
    tuple val(meta), path("*/*/*.intron_grouped_counts.tsv"),           emit: grouped_intron_counts,           optional: true
    tuple val(meta), path("*/*/*.transcript_models.gtf"),               emit: transcript_models,               optional: true
    tuple val(meta), path("*/*/*.transcript_model_reads.tsv.gz"),       emit: transcript_model_reads,          optional: true
    tuple val(meta), path("*/*/*.transcript_model_tpm.tsv"),            emit: transcript_model_tpm,            optional: true
    tuple val(meta), path("*/*/*.transcript_model_counts.tsv"),         emit: transcript_model_counts,         optional: true
    tuple val(meta), path("*/*/*.extended_annotation.gtf"),             emit: extended_gtf,                    optional: true
    tuple val(meta), path("*/*/*.transcript_model_grouped_counts.tsv"), emit: grouped_transcript_model_counts, optional: true
    tuple val(meta), path("*/*/*.transcript_model_grouped_tpm.tsv"),    emit: grouped_transcript_model_tpm,    optional: true
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

    isoquant.py ${args} \\
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
