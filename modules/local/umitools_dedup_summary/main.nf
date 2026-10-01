//
// Sum the per-chunk umi_tools dedup logs of a sample into one duplication summary.
//
// Deduplication runs per chromosome / transcript-group chunk, so the individual
// umi_tools logs only ever hold chunk-level counts. This is the one place the
// sample-level input reads, output reads and duplication rate are reported.
//
process UMITOOLS_DEDUP_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(logs)

    output:
    tuple val(meta), path("*.umi_dedup_summary.tsv"), emit: summary
    path "versions.yml"                              , emit: versions_umitools_dedup_summary, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    summarize_umitools_dedup.sh \\
        ${args} \\
        --sample ${prefix} \\
        --output ${prefix}.umi_dedup_summary.tsv \\
        ${logs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/^GNU Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tinput_reads\\toutput_reads\\tduplicate_reads\\tduplication_rate\\tn_dedup_runs\\n' > ${prefix}.umi_dedup_summary.tsv
    printf '${prefix}\\t0\\t0\\t0\\t0.0000\\t0\\n' >> ${prefix}.umi_dedup_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/^GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}
