process UMI_TOOLS_WHITELIST {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38hbff2b2d_1':
        'quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0' }"

    input:
    tuple val(meta), path(reads)
    val exp_cell_amount
    val regex_pattern

    output:
    tuple val(meta), path(reads), path("*.txt"), emit: whitelist
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        umi_tools \\
            whitelist \\
            --log2stderr \\
            --stdin=${reads} \\
            --bc-pattern ${regex_pattern} \\
            --set-cell-number ${exp_cell_amount} \\
            ${args} > ${prefix}.whitelist.txt 2> ${prefix}.err

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umi_tools_whitelist: \$(echo \$(umi_tools -v 2>&1) | sed 's/^UMI-tools version: //' ))
        END_VERSIONS
        """
        } else {
        """
        umi_tools \\
            whitelist \\
            --log2stderr \\
            --stdin=${reads[0]} \\
            --bc-pattern ${regex_pattern} \\
            --set-cell-number ${exp_cell_amount} \\
            ${args} > ${prefix}.whitelist.txt 2>${prefix}.err
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umi_tools_whitelist: \$(echo \$(umi_tools -v 2>&1) | sed 's/^UMI-tools version: //' )
        END_VERSIONS
        """
    }
}
