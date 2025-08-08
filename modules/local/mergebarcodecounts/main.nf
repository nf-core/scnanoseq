process MERGEBARCODECOUNTS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04':
        'biocontainers/ubuntu:22.04' }"

    input:
    tuple val(meta), path(barcode_counts)

    output:
    tuple val(meta), path("merged_barcode_counts.txt"), emit: barcode_counts
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    awk -F'\t' '{counts[\$1]+=\$2} END {for (b in counts) print b "\t" counts[b]}' ${barcode_counts} \
    | sort -k2,2nr > merged_barcode_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU Awk: \$(awk --version | sed -n '1s/.*GNU Awk \\([0-9.]*\\).*/\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch merged_barcode_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU Awk: \$(awk --version | sed -n '1s/.*GNU Awk \\([0-9.]*\\).*/\1/p')
    END_VERSIONS
    """
}
