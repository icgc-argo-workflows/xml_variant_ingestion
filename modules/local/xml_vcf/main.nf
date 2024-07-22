process XML_VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/pandas' :
        'docker.io/gfeng2023/pandas-pyfaidx:latest' }"
        // 'biocontainers/pandas:2.2.1' }"


    input:
    tuple val(meta), path(xml)
    path (hg19_fa)
    path (hg19_fai)

    output:
    tuple val(meta), path ("*.short_variant.vcf"), emit: short_variant_vcf
    tuple val(meta), path ("*.rearrangement.vcf"), emit: rearrangement_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    shortvariant.py \
        -i ${xml} \
        -r ${hg19_fai} \
        -o ${prefix}.short_variant.vcf

    rearrangement.py \
        -i ${xml} \
        -r ${hg19_fa} \
        -r2 ${hg19_fai} \
        -o ${prefix}.rearrangement.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
//
    // rearrangement.py \\
    //     -i ${xml} \\
    //     -r ${hg19_fa} \\
    //     -r2 ${hg19_fai} \\
        // -o ${prefix}.rearrangement.vcf
