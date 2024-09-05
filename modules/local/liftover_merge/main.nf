process LIFTOVER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::multiqc=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cnv_main)
    tuple val(meta), path(cnv_main_index)
    tuple val(meta), path(cnv_end)
    tuple val(meta), path(cnv_end_index)

    output:
    tuple val(meta), path("*.lifted.vcf")  , emit: vcf_lifted
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    liftover_merge.py \
        -i ${cnv_main} \
        -i2 ${cnv_end} \
        -o ${prefix}.copy_number.lifted.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
