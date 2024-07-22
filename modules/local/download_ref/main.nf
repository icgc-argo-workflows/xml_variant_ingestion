process DOWNLOAD_REF {
    // tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas' :
        'quay.io/biocontainers/pandas' }"

    input:
    val fasta_url
    val fasta_checksum
    val fai_url
    val fai_checksum
    val chain_url
    val chain_checksum

    output:
    path "data/*.fa.gz", emit: fasta_file
    path "data/*.fa.fai", emit: fai_file
    path "data/*.chain.gz", emit: chain_file

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p data
    wget -O data/hg19ToHg38.over.chain.gz ${chain_url}
    wget -O data/GRCh38_hla_decoy_ebv.fa.fai ${fai_url}
    wget -O data/GRCh38_hla_decoy_ebv.fa.gz ${fasta_url}

    echo "${chain_checksum}  data/hg19ToHg38.over.chain.gz" | sha256sum -c -

    if [ \$? -ne 0 ]; then
        echo "LiftOver chain file checksum verification failed. Exiting."
        exit 1
    fi



    echo "${fai_checksum}  data/GRCh38_hla_decoy_ebv.fa.fai" | sha256sum -c -

    if [ \$? -ne 0 ]; then
        echo "fai file checksum verification failed. Exiting."
        exit 1
    fi


    echo "${fasta_checksum}  data/GRCh38_hla_decoy_ebv.fa.gz" | sha256sum -c -

    if [ \$? -ne 0 ]; then
        echo "fasta file checksum verification failed. Exiting."
        exit 1
    fi



    """
}


