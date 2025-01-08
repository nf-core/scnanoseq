process SPLIT_FILE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(unsplit_file)
    val file_ext
    val split_amount

    output:
    tuple val(meta), path("*$file_ext"), emit: split_files
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split -a10 -l ${split_amount} -d --additional-suffix ${file_ext} ${unsplit_file} ${prefix}.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split: \$(echo \$(split --version 2>&1 | head -n1 | sed 's#split (GNU coreutils) ##g'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${file_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split: \$(echo \$(split --version 2>&1 | head -n1 | sed 's#split (GNU coreutils) ##g'))
    END_VERSIONS
    """
}
