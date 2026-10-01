process SPLIT_BC_UMI {
    tag "$meta.id"
    // A single-threaded pass over the fastq: ~20k reads/s, ~30 MB resident.
    label 'process_single'

    conda "conda-forge::python=3.10.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(reads), path(barcode_counts)
    val bc_length
    val umi_length
    val barcode_tag
    val min_reads

    output:
    tuple val(meta), path("*.bc_umi.fastq.gz"), path("*.cdna.fastq.gz"), emit: reads
    tuple val(meta), path("*.flagstat")                                , emit: flagstat
    path "versions.yml"                                                , emit: versions_split_bc_umi, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}${meta.type ? ".${meta.type}" : ""}"
    // An empty tag leaves the barcode where it always came from, the read name.
    def tag_arg = barcode_tag ? "-t ${barcode_tag}" : ""
    // No barcode counts staged leaves the barcode space unbounded, which is the
    // right thing for the called-cells pass: its barcodes are already the knee
    // called list. The all-droplet pass keys on XB, the raw uncorrected barcode,
    // where sequencing error alone mints a droplet per read -- there the counts
    // file plus a minimum bounds the column count the lr-kallisto EM has to
    // cover. The threshold is applied inside split_bc_umi.py rather than by an
    // awk pre-step: this container carries no guaranteed awk.
    def allowlist_arg = barcode_counts ? "-l ${barcode_counts} -m ${min_reads}" : ""
    """
    split_bc_umi.py \\
        -i ${reads} \\
        -o ${prefix} \\
        -b ${bc_length} \\
        -u ${umi_length} \\
        ${tag_arg} \\
        ${allowlist_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}${meta.type ? ".${meta.type}" : ""}"
    """
    echo "" | gzip -c > ${prefix}.bc_umi.fastq.gz
    echo "" | gzip -c > ${prefix}.cdna.fastq.gz
    touch ${prefix}.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
