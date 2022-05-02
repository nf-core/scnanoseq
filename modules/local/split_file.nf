process SPLIT_FILE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(unsplit_file)
    val file_ext

    output:
    // TODO: Make this more generalizable. Gunzip probably a good example
    tuple val(meta), path("*.fastq"), emit: split_files
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: Convert split amount to parameter
    """
    split -a4 -l 10 -d --additional-suffix ${file_ext} ${unsplit_file} ${prefix}.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split: \$(echo \$(split --version 2>&1) | sed 's#split (GNU coreutils) ##g')
    END_VERSIONS
    """
}
