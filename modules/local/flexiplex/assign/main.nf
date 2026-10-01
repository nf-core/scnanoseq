process FLEXIPLEX_ASSIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flexiplex:1.02.7--py312h9b2995c_0':
        'biocontainers/flexiplex:1.02.7--py312h9b2995c_0' }"

    input:
    tuple val(meta), path(reads), path(barcodes)

    output:
    tuple val(meta), path("*flexiplex.fastq.gz")            , emit: reads
    path "versions.yml"                                     , emit: versions_flexiplex_assign, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}${meta.part ? "_part_${meta.part}" : ''}"
    // With -a, flexiplex reports every read and writes CB/CR/UB/UR into the header
    // comment, leaving CB as "-" for a read it could not match to the known list and
    // CR as "-" for a read where it found no barcode region at all. The awk below
    // does two things:
    //
    //   1. Drops the reads with no barcode region. They belong to no droplet, called
    //      or empty, and flexiplex gives them a one-character UMI placeholder that
    //      umi_tools would refuse to bundle alongside a full-width UMI.
    //   2. Derives XB from the two barcodes, so a single tag identifies the droplet
    //      a surviving read came from:
    //
    //        XB = CB   barcode called against the known list (a cell)
    //        XB = CR   no known barcode, but a barcode region was found (a droplet
    //                  below the knee flexiplex-filter called)
    //
    // Doing it here rather than over the aligned bam keeps it in a stream that is
    // already being read and written, spares minimap2 the reads that are about to be
    // discarded, and lets minimap2 -y carry XB into the bam alongside the tags
    // flexiplex wrote itself.
    """
    # Run in assignment mode

    zcat ${reads} | flexiplex \\
        ${args} \\
        -k ${barcodes} \\
        -p ${task.cpus} \\
        | awk '
            BEGIN { FS = "\\t"; OFS = "\\t" }
            # flexiplex repeats the whole header, tags and all, on the "+" line. Both
            # copies have to be rewritten or the two captions disagree and readers that
            # check them (NanoPlot, via biopython) reject the file. keep and xb are
            # decided once, on the header, and reused on the "+" line so the two copies
            # cannot drift; discarding a record discards all four of its lines.
            NR % 4 == 1 {
                cb = ""; cr = ""; keep = 0; xb = ""
                for (i = 2; i <= NF; i++) {
                    if (substr(\$i, 1, 5) == "CB:Z:") cb = substr(\$i, 6)
                    else if (substr(\$i, 1, 5) == "CR:Z:") cr = substr(\$i, 6)
                }
                if (cr != "" && cr != "-") {
                    keep = 1
                    xb = (cb != "" && cb != "-") ? cb : cr
                }
            }
            !keep { next }
            NR % 4 == 1 || NR % 4 == 3 {
                print \$0 "\\tXB:Z:" xb
                next
            }
            { print }
        ' \\
        | gzip -c > ${prefix}.flexiplex.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help |& sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.flexiplex.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex --help 2>/dev/null | sed '1!d ; s/FLEXIPLEX //')
    END_VERSIONS
    """
}
