process SPLIT_SEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::seqkit=2.10.0 conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:0.10.0--1' :
        'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(unsplit_file)
    val file_ext
    val split_amount

    output:
    tuple val(meta), path("output/*$file_ext"), emit: split_files
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Split the file by number of reads
    seqkit -j ${task.cpus} split2 ${args} \\
      -s ${split_amount} --out-dir output --force ${unsplit_file}

    # rename files to have the correct extension
    for file in ./output/*.part_*; do
        if [[ -f "\$file" ]]; then
            base_name=\$(basename "\$file")
            # Remove .gz suffix if present
            if [[ "\$base_name" == *.gz ]]; then
                base_name_no_gz="\${base_name%.gz}"
            else
                base_name_no_gz="\$base_name"
            fi
            # Remove the remaining extension (after the last dot)
            base_name_final="\${base_name_no_gz%.*}"
            # Remove up to .part_
            base_name_final="\${base_name_final#*.part_}"
            # Rename the file
            mv "\$file" "output/${prefix}.\${base_name_final}${file_ext}"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit version | head -n1 | sed 's/seqkit version //g'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.part_001${file_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit version | head -n1 | sed 's/seqkit version //g'))
    END_VERSIONS
    """
}
