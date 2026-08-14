//
// Split a bam into the alignments per-gene dedup can act on and the rest.
//
// umi_tools --per-gene silently drops any read whose gene tag is missing or
// whose assigned-status tag matches --skip-tags-regex. Splitting explicitly
// instead of relying on that means the two record counts must add back up to
// the input, so the positional fallback is verifiable rather than implicit.
//
process SPLIT_GENE_STATUS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.22.1 bioconda::htslib=1.22.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val gene_status_filter

    output:
    tuple val(meta), path("*.gene.bam"), path("*.gene.bam.bai"), emit: gene_bam
    tuple val(meta), path("*.pos.bam") , path("*.pos.bam.bai") , emit: pos_bam
    path "versions.yml"                , emit: versions_split_gene_status, topic: versions

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
        -e '${gene_status_filter}' \\
        -o ${prefix}.gene.bam \\
        -U ${prefix}.pos.bam \\
        ${bam}

    samtools index -@ ${task.cpus} ${prefix}.gene.bam
    samtools index -@ ${task.cpus} ${prefix}.pos.bam

    # The two halves must add back up to the input. If they ever do not, reads
    # would be silently lost between here and the merge, which is the one
    # failure mode of this design that produces plausible-looking output.
    n_in=\$(samtools view -c -@ ${task.cpus} ${bam})
    n_gene=\$(samtools view -c -@ ${task.cpus} ${prefix}.gene.bam)
    n_pos=\$(samtools view -c -@ ${task.cpus} ${prefix}.pos.bam)

    if [ \$(( n_gene + n_pos )) -ne \$n_in ]; then
        echo "ERROR: gene (\$n_gene) + positional (\$n_pos) != input (\$n_in)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene.bam
    touch ${prefix}.gene.bam.bai
    touch ${prefix}.pos.bam
    touch ${prefix}.pos.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
