process NANOPLOT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.46.1--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.46.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.html")                , emit: html
    tuple val(meta), path("*.png") , optional: true, emit: png
    tuple val(meta), path("*.txt")                 , emit: txt
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_file = [ ".fastq.gz", ".fastq", ".fq", ".fq.gz" ].any { ext -> "$ontfile".endsWith(ext) } ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        $input_file

    for nanoplot_file in *.html *.png *.txt *.log
    do
        if [[ -s \$nanoplot_file ]]
        then
            mv \$nanoplot_file ${prefix}_\$nanoplot_file
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_LengthvsQualityScatterPlot_dot.html
    touch ${prefix}_LengthvsQualityScatterPlot_kde.html
    touch ${prefix}_NanoPlot-report.html
    touch ${prefix}_NanoStats.txt
    touch ${prefix}_Non_weightedHistogramReadlength.html
    touch ${prefix}_Non_weightedLogTransformed_HistogramReadlength.html
    touch ${prefix}_WeightedHistogramReadlength.html
    touch ${prefix}_WeightedLogTransformed_HistogramReadlength.html
    touch ${prefix}_Yield_By_Length.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
