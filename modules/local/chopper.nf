process CHOPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::nanofilt=0.10.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.10.0--hcdda2d0_0':
        'biocontainers/chopper:0.10.0--hcdda2d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.filtered.fastq.gz"), emit: reads
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    FILE_PREFIX=${prefix}
    if [ ${params.split_amount} -gt 0 ]; then
        IDX=\$(basename ${reads} | cut -f2 -d'.')
        FILE_PREFIX=\${FILE_PREFIX}.\${IDX}
    fi

    chopper -t ${task.cpus} $args --input $reads | \\
      gzip -c > \${FILE_PREFIX}.filtered.fastq.gz 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$( chopper --version | sed -e "s/chopper //g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$( chopper --version | sed -e "s/chopper //g" )
    END_VERSIONS
    """
}
