#!/usr/bin/env nextflow

/*
========================================================================================
                         smRNA sequencing pipeline
========================================================================================
Ben Pastore
pastore.28@osu.edu
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the smRNA sequencing pipeline is as follows:

    nextflow main.nf -profile cluster --reads /path/to/reads --results /path/to/results

    Mandatory arguments:
    --design                            CSV file with design
    --results                           Path to results directory
    --outprefix                         Outprefix of files 
    -profile [str]                      Configuration profile to use. Can use local / cluster

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

params.bin = "${params.base}/../../bin"
params.index = "${params.base}/../../index"

/*
////////////////////////////////////////////////////////////////////
Check bowtie index is made
////////////////////////////////////////////////////////////////////
*/
genome_fasta = file("${params.genome}")
genome_name = "${genome_fasta.baseName}"
params.bowtie_index_path = "${params.index}/bowtie/${genome_name}"
params.bowtie_index = "${params.bowtie_index_path}/${genome_name}"
bowtie_exists = file(params.bowtie_index_path).exists()

/*
////////////////////////////////////////////////////////////////////
Check tailor index is made
////////////////////////////////////////////////////////////////////
*/
params.tailor_index_path = "${params.index}/tailor/${genome_name}"
params.tailor_index = "${params.tailor_index_path}/${genome_name}"
tailor_exists = file(params.tailor_index_path).exists()

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { BOWTIE_INDEX } from '../../modules/bowtie/main.nf'
include { BOWTIE_ALIGN_GENOME } from '../../modules/bowtie/main.nf'
include { REMOVE_CONTAMINANT } from '../../modules/bowtie/main.nf'
include { RBIND_ALIGNMENT_LOG } from '../../modules/misc/main.nf'
include { COUNT_FEATURES } from '../../modules/counter/main.nf'
include { TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { RBIND_COUNTS } from '../../modules/misc/main.nf'
include { DGE } from '../../modules/misc/main.nf'
include { MASTER_TABLE } from '../../modules/misc/main.nf'
include { TAILOR_INDEX } from '../../modules/tailor/main.nf'
include { TAILOR_MAP } from '../../modules/tailor/main.nf'

/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    /*
     * Parse design file
     */
    Channel
        .fromPath( "${params.design}")
        .splitCsv( header: ['reads', 'condition'], sep: "\t", skip: 1)
        .map{ row -> [ file(row.reads).simpleName, row.reads ] }
        .tap{ reads_ch }
        
    /*
     * Trimming & QC
     */
    if (params.format == 'fastq' & params.trimming) {
        
        TRIM_GALORE( reads_ch )
        TRIM_GALORE.out.collapsed_fa.tap{ fasta_ch }

    } else if (params.format == 'fasta') {

        fasta_ch = reads_ch

    } else {

        exit 1, "Directory with fastq or fasta files must be specified as well as the format of the input reads (See nexflow.config file)"

    }
    
    /*
     * Remove reads mapping to contaminating RNAs
     */
    if (params.contaminant_fa){

        REMOVE_CONTAMINANT( params.contaminant_fa, fasta_ch )
        REMOVE_CONTAMINANT.out.xk_fasta.tap{ fasta_xc_ch }

    } else {

        fasta_xc_ch = fasta_ch

    }
    
    /*
     * Bowtie Index
     */
    if ( bowtie_exists ){
    
        bowtie_index_ch = params.bowtie_index
    
    } else {

        BOWTIE_INDEX( params.genome, params.junction )
        BOWTIE_INDEX.out.tap{ bowtie_index_ch }
    
    }

    /*
     * Bowtie align
     */
    BOWTIE_ALIGN_GENOME( bowtie_index_ch, fasta_xc_ch )
    
    BOWTIE_ALIGN_GENOME
        .out
        .alignment_logs
        .collect()
        .tap{ bowtie_alignment_logs_ch }

    RBIND_ALIGNMENT_LOG( params.outprefix, bowtie_alignment_logs_ch )


    /*
     * Counter
     */
    if (params.run_counter ) {

        BOWTIE_ALIGN_GENOME
            .out
            .bowtie_alignment
            .tap{ counter_input_ch }

        COUNT_FEATURES(
            params.features,
            params.bed, 
            counter_input_ch
            )

        COUNT_FEATURES
            .out
            .normalization_constants
            .tap{ norm_consts_ch }
        
        COUNT_FEATURES
            .out
            .counts
            .tap{ genomic_counts_ch }

    } else {

        BOWTIE_ALIGN_GENOME
            .out
            .normalization_constants
            .tap{ norm_consts_ch }

        Channel
            .empty() 
            .tap{ genomic_counts_ch }
    }

    /*
     * Transcripts
     */
    if (params.run_transcripts) {
        
        fasta_xc_ch
            .join(norm_consts_ch)
            .tap{ transcripts_input_ch }

        TRANSCRIPTS(params.transcripts, transcripts_input_ch)

        TRANSCRIPTS
            .out
            .counts
            .tap{ transcript_counts_ch }
        
    } else {
        Channel
            .empty() 
            .tap{ transcript_counts_ch }
    }
    
    /*
     * Combine genomic and transcripts counts
     */
    
    if (params.run_counter || params.run_transcripts ) {

        genomic_counts_ch
            .mix( transcript_counts_ch )
            .groupTuple()
            .tap{ counts_ch }

        RBIND_COUNTS( counts_ch )

        RBIND_COUNTS
            .out
            .tables
            .tap{ master_table_input }
        
        MASTER_TABLE( params.outprefix, master_table_input.collect() )

        MASTER_TABLE
            .out
            .tables
            .flatten() 
            .tap{ master_table_ch }

    }
    
    /*
     * Differential gene expression analysis
     */
    if (params.comparisons) {

        DGE( params.comparisons, master_table_ch )

    }

    /*
     * Tailor pipeline
     */
    if (params.tailor){

        /*
         * Check if tailor index is build
         */
        if ( tailor_exists ){
    
            tailor_index_ch = params.tailor_index
    
        } else {

            TAILOR_INDEX( params.genome )
            TAILOR_INDEX.out.tap{ tailor_index_ch }
        
        }

        /*
         * Set tailor process input channel and preform alignment 
         */
        
        BOWTIE_ALIGN_GENOME
            .out
            .tailor_input
            .join(norm_consts_ch)
            .tap{ tailor_input_ch }

        TAILOR_MAP(
            params.genome, 
            tailor_index_ch, 
            params.features, 
            params.bed, 
            tailor_input_ch
            )
    }

}






