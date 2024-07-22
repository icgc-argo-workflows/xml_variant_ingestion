#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfcore/variantcall
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfcore/variantcall
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.study_id                    = WorkflowMain.getGenomeAttribute(params, 'study_id')
params.analysis_id                 = WorkflowMain.getGenomeAttribute(params, 'analysis_id')
params.samplesheet                 = WorkflowMain.getGenomeAttribute(params, 'samplesheet')

params.fasta_url                   = WorkflowMain.getGenomeAttribute(params, 'fasta_url')
params.fai_url                     = WorkflowMain.getGenomeAttribute(params, 'fai_url')
params.chain_url                   = WorkflowMain.getGenomeAttribute(params, 'chain_url')

params.api_token                   = WorkflowMain.getGenomeAttribute(params, 'api_token')
params.score_url_upload            = WorkflowMain.getGenomeAttribute(params, 'score_url_upload')
params.song_url_upload             = WorkflowMain.getGenomeAttribute(params, 'song_url_upload')
params.score_url_download          = WorkflowMain.getGenomeAttribute(params, 'score_url_download')
params.song_url_download           = WorkflowMain.getGenomeAttribute(params, 'song_url_download')
params.score_url                   = WorkflowMain.getGenomeAttribute(params, 'score_url')
params.song_url                    = WorkflowMain.getGenomeAttribute(params, 'song_url')

params.tools                       = WorkflowMain.getGenomeAttribute(params, 'tools')
params.outdir                      = WorkflowMain.getGenomeAttribute(params, 'outdir')

params.local_mode                  = WorkflowMain.getGenomeAttribute(params, 'local_mode')
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTCALL } from './workflows/variantcall'

//
// WORKFLOW: Run main nfcore/variantcall analysis pipeline
//
workflow NFCORE_VARIANTCALL {
    VARIANTCALL ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    NFCORE_VARIANTCALL ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
