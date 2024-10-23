process ISOQUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::isoquant=3.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.5.0--hdfd78af_0' :
        'biocontainers/isoquant:3.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(gtf)
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

    // setting custom home via export (see issue #30)
    if ( !group_category?.trim() ){
        """
        export HOME=\$(pwd)

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
        export HOME=\$(pwd)

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
