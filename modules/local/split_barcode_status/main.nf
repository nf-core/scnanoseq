//
// Split a bam into the alignments UMI deduplication can act on and the rest.
//
// Since flexiplex is run with -a it reports every read, including reads where no
// barcode region was found at all. Those carry XB "-" and a UMI of "-", so they
// have nothing to deduplicate against and, pooled into a single cell, their
// one-character UMI sits alongside the full-width UMIs of the uncalled droplets --
// which trips umi_tools' assertion that every UMI in a bundle is the same length.
// Splitting them off keeps them in the output without letting them reach umi_tools.
//
process SPLIT_BARCODE_STATUS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.22.1 bioconda::htslib=1.22.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val barcode_status_filter

    output:
    tuple val(meta), path("*.barcoded.bam") , path("*.barcoded.bam.bai") , optional: true, emit: barcoded_bam
    tuple val(meta), path("*.nobarcode.bam"), path("*.nobarcode.bam.bai"), optional: true, emit: nobarcode_bam
    path "versions.yml"                     , emit: versions_split_barcode_status, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // samtools view preserves the input ordering, so neither output needs a
    // re-sort before umi_tools sees it.
    """
    samtools view \\
        --threads ${task.cpus} \\
        -b \\
        -e '${barcode_status_filter}' \\
        -o ${prefix}.barcoded.bam \\
        -U ${prefix}.nobarcode.bam \\
        ${bam}

    # The two halves must add back up to the input. If they ever do not, reads would
    # be silently lost between here and the merge, which is the one failure mode of
    # this design that produces plausible-looking output.
    n_in=\$(samtools view -c -@ ${task.cpus} ${bam})
    n_barcoded=\$(samtools view -c -@ ${task.cpus} ${prefix}.barcoded.bam)
    n_nobarcode=\$(samtools view -c -@ ${task.cpus} ${prefix}.nobarcode.bam)

    if [ \$(( n_barcoded + n_nobarcode )) -ne \$n_in ]; then
        echo "ERROR: barcoded (\$n_barcoded) + no-barcode (\$n_nobarcode) != input (\$n_in)" >&2
        exit 1
    fi

    # umi_tools aborts on an empty bam, and a contig can legitimately hold nothing but
    # unbarcoded reads. Drop an empty half instead of emitting it; both outputs are
    # optional and the caller mixes whatever arrives.
    if [ "\$n_barcoded" -eq 0 ]; then
        rm -f ${prefix}.barcoded.bam
    else
        samtools index -@ ${task.cpus} ${prefix}.barcoded.bam
    fi

    if [ "\$n_nobarcode" -eq 0 ]; then
        rm -f ${prefix}.nobarcode.bam
    else
        samtools index -@ ${task.cpus} ${prefix}.nobarcode.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.barcoded.bam
    touch ${prefix}.barcoded.bam.bai
    touch ${prefix}.nobarcode.bam
    touch ${prefix}.nobarcode.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
