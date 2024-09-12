/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XML_VCF } from '../modules/local/xml_vcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_SV } from '../modules/nf-core/picard/liftovervcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_RA } from '../modules/nf-core/picard/liftovervcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_CN } from '../modules/nf-core/picard/liftovervcf/main'
include { PREP_META } from '../modules/local/prep/metadata/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_SV } from '../modules/local/payload/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_RA } from '../modules/local/payload/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_CN } from '../modules/local/payload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_SV } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_RA } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_CN } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCALL {

    // take:
    // xml

    // main:

    ch_versions = Channel.empty()

    // Prepare metadata, adapted from SanityCheck
    PREP_META(
        file(params.experiment_info_tsv, checkIfExists: true),
        params.api_token,
        params.song_url,
        params.clinical_url,
        params.skip_duplicate_check
    )
    ch_versions = ch_versions.mix(PREP_META.out.versions)

    PREP_META.out.updated_experiment_info_tsv
    .collectFile(keepHeader: true, name: 'updated_sample.tsv')
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        [
            study_id: row.program_id,
            id: row.sample_id,
            donor_id: row.donor_id,
            specimen_id: row.specimen_id,
            experimental_strategy: row.experimental_strategy
        ]
    }.set{meta_ch}

    // XML to VCF conversion

    meta_ch.combine(Channel.fromPath(params.xml)).set{xml_ch}

    // xml_ch = Channel.fromPath(params.xml)
    //             .map { path -> [ [id: '2001205343'], path ] }

    XML_VCF (
        xml_ch,
        Channel.fromPath(params.hg19_ref_fa, checkIfExists: true),
        Channel.fromPath(params.hg19_ref_fai, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(XML_VCF.out.versions)

    // VCF lift over
    hg38_ref_ch = Channel.fromPath(params.hg38_ref_fa, checkIfExists: true)
                            .map{ path -> [ [id: 'fasta'], path ] }

    hg38_ref_dict = Channel.fromPath(params.hg38_ref_dict, checkIfExists: true)
                            .map{ path -> [ [id: 'dict'], path ] }

    hg19_to_hg38_chain_ch = Channel.fromPath(params.hg19_to_hg38_chain, checkIfExists: true)
                            .map{ path -> [ [id: 'chain'], path ] }

    // Short Vraint \\
    // lift over
    PICARD_LIFTOVERVCF_SV (
        XML_VCF.out.short_variant_vcf,
        hg38_ref_dict,
        hg38_ref_ch,
        hg19_to_hg38_chain_ch
    )
    ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_SV.out.versions)

    //Payload Generation
    PICARD_LIFTOVERVCF_SV.out.vcf_lifted
    .combine(PICARD_LIFTOVERVCF_SV.out.vcf_lifted_index)
    .combine(meta_ch)
    .combine(PREP_META.out.updated_experiment_info_tsv)
    .map{
        metaA, vcf, metaB, index, meta, metadata_analysis ->
        [
            meta, [vcf, index], metadata_analysis
        ]
    }.set{vcf_and_index_sv}

    PAYLOAD_VARIANT_CALL_SV (
        vcf_and_index_sv,
        Channel.empty()
        .mix(XML_VCF.out.versions)
        .mix(PICARD_LIFTOVERVCF_SV.out.versions)
        .collectFile(name: 'collated_versions.yml')
    )
    ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_SV.out.versions)

    // PAYLOAD_VARIANT_CALL_SV.out.payload_files.subscribe { println("payload output: ${it}") }

    // Upload
    SONG_SCORE_UPLOAD_SV(PAYLOAD_VARIANT_CALL_SV.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
    ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_SV.out.versions)

    // Rearrangement \\
    // lift over
    PICARD_LIFTOVERVCF_RA (
        XML_VCF.out.rearrangement_vcf,
        hg38_ref_dict,
        hg38_ref_ch,
        hg19_to_hg38_chain_ch
    )
    ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_RA.out.versions)

    //Payload Generation
    PICARD_LIFTOVERVCF_RA.out.vcf_lifted
    .combine(PICARD_LIFTOVERVCF_RA.out.vcf_lifted_index)
    .combine(meta_ch)
    .combine(PREP_META.out.updated_experiment_info_tsv)
    .map{
        metaA, vcf, metaB, index, meta, metadata_analysis ->
        [
            meta, [vcf, index], metadata_analysis
        ]
    }.set{vcf_and_index_ra}

    PAYLOAD_VARIANT_CALL_RA (
        vcf_and_index_ra,
        Channel.empty()
        .mix(XML_VCF.out.versions)
        .mix(PICARD_LIFTOVERVCF_RA.out.versions)
        .collectFile(name: 'collated_versions.yml')
    )
    ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_RA.out.versions)

    // Upload
    SONG_SCORE_UPLOAD_RA(PAYLOAD_VARIANT_CALL_RA.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
    ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_RA.out.versions)


    // Copy Number \\
    // lift over
    PICARD_LIFTOVERVCF_CN (
        XML_VCF.out.copy_number_vcf,
        hg38_ref_dict,
        hg38_ref_ch,
        hg19_to_hg38_chain_ch
    )
    ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_CN.out.versions)

    //Payload Generation
    PICARD_LIFTOVERVCF_CN.out.vcf_lifted
    .combine(PICARD_LIFTOVERVCF_CN.out.vcf_lifted_index)
    .combine(meta_ch)
    .combine(PREP_META.out.updated_experiment_info_tsv)
    .map{
        metaA, vcf, metaB, index, meta, metadata_analysis ->
        [
            meta, [vcf, index], metadata_analysis
        ]
    }.set{vcf_and_index_cn}

    PAYLOAD_VARIANT_CALL_CN (
        vcf_and_index_cn,
        Channel.empty()
        .mix(XML_VCF.out.versions)
        .mix(PICARD_LIFTOVERVCF_CN.out.versions)
        .collectFile(name: 'collated_versions.yml')
    )
    ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_CN.out.versions)

    // Upload
    SONG_SCORE_UPLOAD_CN(PAYLOAD_VARIANT_CALL_CN.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
    ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_CN.out.versions)

    emit:

     versions = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
