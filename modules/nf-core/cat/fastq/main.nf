process CAT_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads
    tuple val("${task.process}"), val("python"), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_list_str = reads instanceof List ? reads.join(' ') : reads.toString()
    def single_end_flag = meta.single_end ? '--single_end' : ''

    """
    cat_fastq.py \\
        --prefix ${prefix} \\
        ${single_end_flag} \\
        --reads ${read_list_str}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect { item -> item.toString() } : [reads.toString()]
    if (meta.single_end) {
        if (readList.size >= 1) {
            """
            echo '' | gzip > ${prefix}.merged.fastq.gz
            """
        } else {
            error("Could not find any FASTQ files to concatenate in the process input")
        }
    }
    else {
        if (readList.size >= 2) {
            """
            echo '' | gzip > ${prefix}_1.merged.fastq.gz
            echo '' | gzip > ${prefix}_2.merged.fastq.gz
            """
        } else {
            error("Could not find any FASTQ file pairs to concatenate in the process input")
        }
    }
}
