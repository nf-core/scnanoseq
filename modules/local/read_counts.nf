process READ_COUNTS {
    label 'process_low'

    
    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

   //conda "conda-forge::sed=4.7"
   //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   //    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
   //    'ubuntu:20.04' }"

    input:
    path raw_fastqc
    path trim_fastqc
    path preextract_fastqc
    path correct_tsv

    output:
    path "read_counts.csv"        , emit: read_counts

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    generate_read_counts.sh \\
        $args \\
        --input ./ \\
        --output read_counts.csv
    """
}
