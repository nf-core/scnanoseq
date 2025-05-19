process READ_COUNTS {
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

    input:
    path raw_fastqc
    path trim_fastqc
    path preextract_fastqc
    path correct_tsv

    output:
    path "read_counts.csv" , emit: read_counts
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    generate_read_counts.sh \\
        $args \\
        --input ./ \\
        --output read_counts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | head -n2 | tail -n1 |  sed -n 's/.*(v\\([^)]*\\)).*/\\1/p')
    END_VERSIONS
    """

    stub:
    """
    touch read_counts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | head -n2 | tail -n1 |  sed -n 's/.*(v\\([^)]*\\)).*/\\1/p')
    END_VERSIONS
    """
}
