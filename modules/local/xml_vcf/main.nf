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
    tuple val(meta), path ("*.short_variant.vcf"), optional:true, emit: short_variant_vcf
    tuple val(meta), path ("*.rearrangement.vcf"), optional:true, emit: rearrangement_vcf
    tuple val(meta), path ("*.copy_number.vcf"), optional:true, emit: copy_number_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Check for 'short-variant' in the XML file and capture output
    shortVariantOutput=\$(grep -m 1 'short-variant ' ${xml} 2>&1) || true

    # Check for 'rearrangement' in the XML file and capture output
    rearrangementOutput=\$(grep -m 1 'rearrangement ' ${xml} 2>&1) || true

    # Check for 'copy-number-alteration' in the XML file and capture output
    copyNumberOutput=\$(grep -m 1 'copy-number-alteration ' ${xml} 2>&1) || true

    # Run shortvariant.py if 'short-variant ' is found in the output
    if [ -n "\$shortVariantOutput" ]; then
        shortvariant.py \\
            -i ${xml} \\
            -r ${hg19_fai} \\
            -o ${prefix}.short_variant.vcf
    fi

    # Run rearrangement.py if 'rearrangement ' is found in the output
    if [ -n "\$rearrangementOutput" ]; then
        rearrangement.py \\
            -i ${xml} \\
            -r ${hg19_fa} \\
            -r2 ${hg19_fai} \\
            -o ${prefix}.rearrangement.vcf
    fi

    # Run copynumber.py if 'copy-number-alteration ' is found in the output
    if [ -n "\$copyNumberOutput" ]; then
        copynumber.py \\
            -i ${xml} \\
            -r ${hg19_fa} \\
            -r2 ${hg19_fai} \\
            -o ${prefix}.copy_number.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
