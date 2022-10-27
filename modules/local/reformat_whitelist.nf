process REFORMAT_WHITELIST {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(whitelist)

    output:
    tuple val(meta), path("*.bc_count.txt")        , emit: whitelist_bc_count
    tuple val(meta), path("*.bc_count.txt.bc_only"), emit: whitelist_bc_list
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    reformat_whitelist.py \\
        -i ${whitelist} \\
        -o ${prefix}.bc_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat_whitelist: v1.0 
    END_VERSIONS
    """
}
