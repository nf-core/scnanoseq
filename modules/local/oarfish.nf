process OARFISH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::oarfish=0.6.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/oarfish:0.6.5--h43eeafb_0' :
        'biocontainers/oarfish:0.6.5--h43eeafb_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*features.tsv.gz") , emit: features
    tuple val(meta), path("*barcodes.tsv.gz") , emit: barcodes
    tuple val(meta), path("*matrix.mtx.gz")    , emit: mtx
    tuple val(meta), path("*meta_info.json")  , emit: meta_info
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    oarfish \\
        --output ${prefix} \\
        --alignments $bam \\
        --threads ${task.cpus} \\
        ${args}

    mv *features.txt features.tsv
    mv *barcodes.txt barcodes.tsv

    grep '^%' *count.mtx > matrix.mtx
    grep -v '^%' *count.mtx | awk '{print \$2" "\$1" "\$3}' >> matrix.mtx

    for tsv_file in *features.tsv *barcodes.tsv *matrix.mtx
    do
        gzip \$tsv_file
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version | sed 's#oarfish ##g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.features.tsv.gz
    touch ${prefix}.barcodes.tsv.gz
    touch ${prefix}.matrix.mtx.gz
    touch ${prefix}.meta_info.json
    """
}
