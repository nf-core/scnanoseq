process SPLIT_FASTA {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)

    output:
    path "*.split.fa"   , emit: split_fasta
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    awk '/^>/{chrom=(split(substr(\$0,2), a, " ")); filename=( a[1] ".split.fa"); print > filename; next}{print >> filename}' $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version) | sed 's/^.*cat (GNU coreutils) //; s/ .*//')
    END_VERSIONS
    """
}
