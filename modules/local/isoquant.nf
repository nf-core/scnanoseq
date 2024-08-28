process ISOQUANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::isoquant=3.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.5.0--hdfd78af_0' :
        'biocontainers/isoquant:3.5.0--hdfd78af_0' }"

    // setting custom home mount (see issue #30)
    containerOptions "--bind ${task.workDir}:/home"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_gtf), path(gtf)
    tuple val(meta_fa), path(fasta)
    tuple val(meta_fai), path(fai)
    val group_category

    output:
    tuple val(meta), path("*.gene_counts.tsv")      , emit: gene_count_mtx
    tuple val(meta), path("*.transcript_counts.tsv"), emit: transcript_count_mtx
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ( !group_category?.trim() ){
        """
        isoquant.py ${args} \\
                    --threads $task.cpus \\
                    --datatype nanopore \\
                    --reference $fasta \\
                    --genedb $gtf \\
                    --bam $bam \\
                    -o .

        mv OUT/OUT.gene_counts.tsv ${prefix}.gene_counts.tsv
        mv OUT/OUT.transcript_counts.tsv ${prefix}.transcript_counts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
        END_VERSIONS
        """
    } else {
        """
        isoquant.py ${args} \\
                    --threads $task.cpus \\
                    --data_type nanopore \\
                    --reference $fasta \\
                    --genedb $gtf \\
                    --bam $bam \\
                    -o . \\
                    --read_group $group_category

        mv OUT/OUT.gene_grouped_counts.tsv ${prefix}.gene_counts.tsv
        mv OUT/OUT.transcript_grouped_counts.tsv ${prefix}.transcript_counts.tsv

        sed -i "1s/#//" ${prefix}.gene_counts.tsv
        sed -i "1s/#//" ${prefix}.transcript_counts.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
        END_VERSIONS
        """

    }
}
