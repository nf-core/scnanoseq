process UNZIP_FILE {
    tag "$meta.id"
    label 'process_low'

    // TODO: Figure out what needs to be done for gunzip
    //conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //    'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(zipped_file)

    output:
    // TODO: Figure out how to generalize this
    tuple val(meta), path("*.fastq"), emit: unzipped_file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip $zipped_file

    cat <<-END_VERSIONSi > versions.yml
    "${task.process}":
        unzip_file: \$(echo \$( gunzip -V | grep gzip | sed 's#gzip ##g' ))
    END_VERSIONS
    """
}
