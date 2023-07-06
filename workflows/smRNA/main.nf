#!/usr/bin/env nextflow

/*
========================================================================================
                         smRNA sequencing pipeline (can be used for riboseq as well)
========================================================================================
Ben Pastore
pastore.28@osu.edu
----------------------------------------------------------------------------------------
HELLO WORLD
To do: 
Add option to do alignment with STAR, might be useful for analyzing riboseq data
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

/*
////////////////////////////////////////////////////////////////////
set path to bin and index
////////////////////////////////////////////////////////////////////
*/
params.bin = "${params.base}/../../bin"
params.index = "${params.base}/../../index"

/*
////////////////////////////////////////////////////////////////////
validate results
////////////////////////////////////////////////////////////////////
*/
if (params.results)   { ; } else { exit 1, 'Results path not specified!' }

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { ARTIFACTS_FILTER } from '../../modules/filter/main.nf'
include { TRIM_POLYA } from '../../modules/filter/main.nf'
include { TRIM_UMI } from '../../modules/filter/main.nf'
include { BOWTIE_INDEX } from '../../modules/bowtie/main.nf'
include { BOWTIE_ALIGN_GENOME } from '../../modules/bowtie/main.nf'
include { BOWTIE_ALIGNMENT_MASTER_TABLE } from '../../modules/bowtie/main.nf'
include { REMOVE_CONTAMINANT } from '../../modules/bowtie/main.nf'
include { RBIND_ALIGNMENT_LOG } from '../../modules/misc/main.nf'
include { COUNT_FEATURES } from '../../modules/counter/main.nf'
include { INDEX_TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { RBIND_COUNTS } from '../../modules/misc/main.nf'
include { DGE } from '../../modules/misc/main.nf'
include { MASTER_TABLE } from '../../modules/misc/main.nf'
include { TAILOR_INDEX } from '../../modules/tailor/main.nf'
include { TAILOR_MAP } from '../../modules/tailor/main.nf'
include { DESIGN_INPUT } from '../../modules/misc/main.nf'
include { TRIMGALORE_INPUT } from '../../modules/misc/main.nf'
include { ALIGN_SPIKEIN } from '../../modules/spikein/main.nf'
include { NORMALIZE_SPIKEIN } from '../../modules/spikein/main.nf'
include { MERGE_BW } from '../../modules/deeptools/main.nf'

/*
////////////////////////////////////////////////////////////////////
Subworkflow
////////////////////////////////////////////////////////////////////
*/
workflow TRIMGALORE {

    if (params.fastq)   { ; } else { exit 1, 'Comma separaterd list of fastq files not specified. Use --fastq fq1,fq2 to specify this paramter!' }
    
    TRIMGALORE_INPUT( params.fastq )

    reads_ch = TRIMGALORE_INPUT
        .out
        .fastq_ch
        .splitCsv( header: ['condition', 'reads'], sep: ",", skip: 1)
        .map{ row -> [ row.condition, row.reads ] }
    
    if (params.format == 'fastq' & params.trimming) {
        
        TRIM_GALORE( reads_ch )
        fasta_ch = TRIM_GALORE.out.collapsed_fa

    } else {
        
        fasta_ch = reads_ch

    }

    if (params.trim_umi) {

        TRIM_UMI( fasta_ch )
        
        fasta_ch = TRIM_UMI.out.fasta
    
    }
    
    if (params.artifacts_filter) {

        ARTIFACTS_FILTER( fasta_ch )

        fasta_ch = ARTIFACTS_FILTER.out.fasta
        
    }

    if (params.trim_polyA) {

        TRIM_POLYA( fasta_ch )

        fasta_ch = TRIM_POLYA.out.fasta

    }

    if (params.contaminant){

        REMOVE_CONTAMINANT( params.contaminant, fasta_ch )

        fasta_xc_ch = REMOVE_CONTAMINANT.out.xk_fasta

    } else {

        fasta_xc_ch = fasta_ch

    }


}

/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    /*
    ////////////////////////////////////////////////////////////////////
    Validate mandatory inputs (design, genome, junctions, results, outprefix)
    ////////////////////////////////////////////////////////////////////
    */
    if (params.design)    { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
    if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
    if (params.junctions) { ch_junction = file(params.junctions, checkIfExists: true) } else { exit 1, 'Junction fasta file not specified!' }
    if (params.outprefix) { ; } else {'Outprefix not specified! Defaulting to smRNA_analysis'; params.outprefix = 'smRNA_analysis' }

    /*
    ////////////////////////////////////////////////////////////////////
    Check bowtie index is made
    ////////////////////////////////////////////////////////////////////
    */
    genome_fasta = file("${params.genome}")
    genome_name = "${genome_fasta.baseName}"
    bowtie_index = "${params.index}/bowtie/${genome_name}"
    bowtie_exists = file(bowtie_index).exists()

    /*
    ////////////////////////////////////////////////////////////////////
    Check tailor index is made
    ////////////////////////////////////////////////////////////////////
    */
    tailor_index = "${params.index}/tailor/${genome_name}"
    tailor_exists = file(tailor_index).exists()

    /*
     * Parse design file
     */
    DESIGN_INPUT( params.design )

    reads_ch = DESIGN_INPUT
        .out
        .fastq_ch
        .splitCsv( header: ['condition', 'reads'], sep: ",", skip: 1)
        .map{ row -> [ row.condition, row.reads ] }

    replicates_ch = DESIGN_INPUT
        .out
        .condition_ch
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.simple_name, row.group ] }
        
    /*
     * Trimming & QC
     */
    if (params.format == 'fastq' & params.trimming) {
        
        TRIM_GALORE( reads_ch )
        fasta_ch = TRIM_GALORE.out.collapsed_fa

    } else {
        
        fasta_ch = reads_ch

    }

    if (params.trim_umi) {

        TRIM_UMI( fasta_ch )
        
        fasta_ch = TRIM_UMI.out.fasta
    
    }
    
    if (params.artifacts_filter) {

        ARTIFACTS_FILTER( fasta_ch )

        fasta_ch = ARTIFACTS_FILTER.out.fasta

    }

    if (params.trim_polyA) {

        TRIM_POLYA( fasta_ch )

        fasta_ch = TRIM_POLYA.out.fasta

    }
    
    /*
     * Remove reads mapping to contaminating RNAs
     */
    if (params.contaminant){

        REMOVE_CONTAMINANT( params.contaminant, fasta_ch )
        fasta_xc_ch = REMOVE_CONTAMINANT.out.xk_fasta

    } else {

        fasta_xc_ch = fasta_ch

    }
    
    /*
     * Bowtie Index
     */
    if ( bowtie_exists ){
    
        bowtie_index_ch = "${bowtie_index}/${genome_name}"
        bowtie_chrom_sizes = "${bowtie_index}/${genome_name}_chrom_sizes"
    
    } else {

        BOWTIE_INDEX( params.genome, params.junctions, "${bowtie_index}/${genome_name}", "${bowtie_index}/${genome_name}_chrom_sizes", genome_name)
        bowtie_index_ch = BOWTIE_INDEX.out.bowtie_index
        bowtie_chrom_sizes = BOWTIE_INDEX.out.bowtie_chrom_sizes
    }

    /*
     * Bowtie align
     */
    BOWTIE_ALIGN_GENOME( bowtie_index_ch, fasta_xc_ch, bowtie_chrom_sizes )
    
    bowtie_alignment_logs_ch = BOWTIE_ALIGN_GENOME
                                .out
                                .alignment_logs
                                .collect()

    RBIND_ALIGNMENT_LOG( params.outprefix, bowtie_alignment_logs_ch )

    /* 
     * make table with chrom, start, end, seq, strand, sample1_rpm, sample2_rpm......
     */
    BOWTIE_ALIGNMENT_MASTER_TABLE( params.outprefix, BOWTIE_ALIGN_GENOME.out.bowtie_rpm.collect() )
    
    /*
     * Counter
     */
    if ( params.features ) {

        counter_input_ch = BOWTIE_ALIGN_GENOME.out.bowtie_alignment

        COUNT_FEATURES(
            params.features,
            params.bed, 
            counter_input_ch
            )

        norm_consts_ch = COUNT_FEATURES.out.normalization_constants
        
        genomic_counts_ch = COUNT_FEATURES.out.counts

    } else {

        norm_consts_ch = BOWTIE_ALIGN_GENOME.out.normalization_constants
        genomic_counts_ch = Channel.empty()
    }

    /*
     * Transcripts
     */
    if (params.transcripts) {

        INDEX_TRANSCRIPTS( params.transcripts )

        transcripts_input_ch = fasta_xc_ch.join(norm_consts_ch)

        TRANSCRIPTS(params.transcripts, transcripts_input_ch, INDEX_TRANSCRIPTS.out.transcript_index_path_ch)

        transcript_counts_ch = TRANSCRIPTS.out.counts
        
    } else {

        transcript_counts_ch = Channel.empty()
    
    }
    
    /*
     * Combine genomic and transcripts counts
     */
    if (params.features || params.transcripts ) {

        counts_ch = genomic_counts_ch
            .mix( transcript_counts_ch )
            .groupTuple()
            .unique()

        RBIND_COUNTS( counts_ch )

        master_table_input = RBIND_COUNTS.out.tables
        
        MASTER_TABLE( params.outprefix, master_table_input.collect() )

        /*
        * Spike-ins
        */
        if (params.spikein) {

            ALIGN_SPIKEIN( params.spikein,  fasta_xc_ch)

            spikein_quant_ch = ALIGN_SPIKEIN.out.spikein_quant_ch

            unnormalized_counts = MASTER_TABLE.out.unnormalized_master_table_ch

            spikein_fasta = file("${params.spikein}")
            spikein_name = "${spikein_fasta.simpleName}"

            NORMALIZE_SPIKEIN( unnormalized_counts, spikein_quant_ch.collect(), spikein_name )

            master_table_ch = MASTER_TABLE
                .out
                .tables
                .concat( NORMALIZE_SPIKEIN.out.counts_normalized_spikein_ch )
                .flatten()

        } else {

            master_table_ch = MASTER_TABLE
                .out
                .tables
                .flatten()
        }
    }
    
    /*
     * Differential gene expression analysis
     */
    if (params.dge) {

        /*
        * Parse comparisons file
        */
        conditions = DESIGN_INPUT
            .out
            .condition_ch
            .splitCsv( header: ['sample', 'condition'], sep: ",", skip: 1)
            .map{ row -> [ row.condition, row.sample ] }
            .groupTuple()

        comparisons = Channel
            .fromPath( "${params.dge}")
            .splitCsv( header: ['x', 'y'], sep: "\t", skip: 1)
            .map{ row -> [ row.x, row.y ] }

        comparisons_ch = comparisons
            .combine(conditions, by : 0)
            .map{ it -> it[1,0,2]}
            .combine(conditions, by : 0)
            .map{ it -> it[1,0,2,3]}
        
        dge_input_ch = comparisons_ch.combine(master_table_ch)

        DGE( dge_input_ch )

    }

    /*
     * Tailor pipeline
     */
    if (params.tailor && params.features && params.bed){

        /*
         * Check if tailor index is build
         */
        if ( tailor_exists ){
    
            tailor_index_ch = "${tailor_index}/${genome_name}"

        } else {

            TAILOR_INDEX( params.genome,  "${tailor_index}/${genome_name}", genome_name )
            tailor_index_ch = TAILOR_INDEX.out.tailor_index
        
        }

        /*
         * Set tailor process input channel and preform alignment 
         */
        
        tailor_input_ch = BOWTIE_ALIGN_GENOME
            .out
            .tailor_input
            .join(norm_consts_ch)

        TAILOR_MAP(
            params.genome, 
            tailor_index_ch, 
            params.features, 
            params.bed, 
            tailor_input_ch
            )
    }
}