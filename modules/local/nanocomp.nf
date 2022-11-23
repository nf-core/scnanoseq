process NANOCOMP {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::nanocomp=1.20.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp:1.20.0--pyhdfd78af_0':
        'quay.io/biocontainers/nanocomp:1.20.0--pyhdfd78af_0' }"

    input:
    path ont_files

    output:
    path "*.html"       , emit: html
    path "*.txt"        , emit: txt
    path "*.log"        , emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_files = ("$ont_files".contains(".fastq.gz") || "$ont_files".contains(".fq.gz")) ? "--fastq ${ont_files}" :
        ("$ont_files".contains(".bam")) ? "--bam ${ont_files}" : ''

    """
    NanoComp \\
    $args \\
    -t $task.cpus \\
    $input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocomp: \$(echo \$(NanoComp --version 2>&1) | sed 's/^.*NanoComp //; s/ .*\$//')
    END_VERSIONS
    """
}
