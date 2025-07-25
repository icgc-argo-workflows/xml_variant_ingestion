/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XML_VCF } from '../modules/local/xml_vcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_SNV } from '../modules/nf-core/picard/liftovervcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_INDEL } from '../modules/nf-core/picard/liftovervcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_RA } from '../modules/nf-core/picard/liftovervcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_CN } from '../modules/nf-core/picard/liftovervcf/main'
include { PREP_META } from '../modules/local/prep/metadata/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_SNV } from '../modules/local/payload/vcf/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_INDEL } from '../modules/local/payload/vcf/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_RA } from '../modules/local/payload/vcf/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_CN } from '../modules/local/payload/vcf/main'
include { PAYLOAD_XML_SUPPLEMENT } from '../modules/local/payload/supplement/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_XML } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_SNV } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_INDEL } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_RA } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UPLOAD_CN } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCALL {

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

    // Submit original XML

    meta_ch.combine(Channel.fromPath(params.xml)).set{xml_ch}

    PAYLOAD_XML_SUPPLEMENT (
        xml_ch.combine(PREP_META.out.updated_experiment_info_tsv)
    )

    // Upload
    if (params.test_analysis_id) {
        input_analysis_id = [[],params.test_analysis_id]
    } else {
       SONG_SCORE_UPLOAD_XML( PAYLOAD_XML_SUPPLEMENT.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
       ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_XML.out.versions)
       input_analysis_id= SONG_SCORE_UPLOAD_XML.out.analysis_id
    }

    // XML to VCF conversion
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

    
    // SNV \\
    if ("SNV" in params.vcf){
        // lift over
        PICARD_LIFTOVERVCF_SNV (
            XML_VCF.out.snv_vcf,
            hg38_ref_dict,
            hg38_ref_ch,
            hg19_to_hg38_chain_ch
        )
        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_SNV.out.versions)

        //Payload Generation
        PICARD_LIFTOVERVCF_SNV.out.vcf_lifted
        .combine(PICARD_LIFTOVERVCF_SNV.out.vcf_lifted_index)
        .combine(meta_ch)
        .combine(PREP_META.out.updated_experiment_info_tsv)
        .map{
            metaA, vcf, metaB, index, meta, metadata_analysis ->
            [
                meta, [vcf, index], metadata_analysis
            ]
        }.set{vcf_and_index_snv}

        PAYLOAD_VARIANT_CALL_SNV (
            vcf_and_index_snv,
            input_analysis_id,
            Channel.empty()
            .mix(XML_VCF.out.versions)
            .mix(PICARD_LIFTOVERVCF_SNV.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_SNV.out.versions)

        // Upload
        if ("UPLOAD" in params.vcf) {
            SONG_SCORE_UPLOAD_SNV(PAYLOAD_VARIANT_CALL_SNV.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
            ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_SNV.out.versions)
        }
    }
    // INDEL \\
    if ("INDEL" in params.vcf){
        // lift over
        PICARD_LIFTOVERVCF_INDEL (
            XML_VCF.out.indel_vcf,
            hg38_ref_dict,
            hg38_ref_ch,
            hg19_to_hg38_chain_ch
        )
        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_INDEL.out.versions)

        //Payload Generation
        PICARD_LIFTOVERVCF_INDEL.out.vcf_lifted
        .combine(PICARD_LIFTOVERVCF_INDEL.out.vcf_lifted_index)
        .combine(meta_ch)
        .combine(PREP_META.out.updated_experiment_info_tsv)
        .map{
            metaA, vcf, metaB, index, meta, metadata_analysis ->
            [
                meta, [vcf, index], metadata_analysis
            ]
        }.set{vcf_and_index_indel}

        PAYLOAD_VARIANT_CALL_INDEL (
            vcf_and_index_indel,
            input_analysis_id,
            Channel.empty()
            .mix(XML_VCF.out.versions)
            .mix(PICARD_LIFTOVERVCF_INDEL.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_INDEL.out.versions)

        // Upload
        if ("UPLOAD" in params.vcf)  {
        SONG_SCORE_UPLOAD_INDEL(PAYLOAD_VARIANT_CALL_INDEL.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_INDEL.out.versions)
        }
    }
    // Rearrangement \\
    if ("RA" in params.vcf){
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
            input_analysis_id,
            Channel.empty()
            .mix(XML_VCF.out.versions)
            .mix(PICARD_LIFTOVERVCF_RA.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_RA.out.versions)

        // // Upload
        if ("UPLOAD" in params.vcf)  {
        SONG_SCORE_UPLOAD_RA(PAYLOAD_VARIANT_CALL_RA.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_RA.out.versions)
        }
    }
    // Copy Number \\
    if ("CNV" in params.vcf){
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
            input_analysis_id,
            Channel.empty()
            .mix(XML_VCF.out.versions)
            .mix(PICARD_LIFTOVERVCF_CN.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL_CN.out.versions)

        // Upload
        if ("UPLOAD" in params.vcf)  {
        SONG_SCORE_UPLOAD_CN(PAYLOAD_VARIANT_CALL_CN.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD_CN.out.versions)
        }
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
