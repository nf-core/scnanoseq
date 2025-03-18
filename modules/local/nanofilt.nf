process NANOFILT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::nanofilt=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanofilt:2.8.0--py_0':
        'biocontainers/nanofilt:2.8.0--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.filtered.fastq")   , emit: reads
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
    cat $reads | NanoFilt $args > \${FILE_PREFIX}.filtered.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanofilt: \$( NanoFilt --version | sed -e "s/NanoFilt //g" )
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.filtered.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanofilt: \$( NanoFilt --version | sed -e "s/NanoFilt //g" )
    END_VERSIONS
    """
}
