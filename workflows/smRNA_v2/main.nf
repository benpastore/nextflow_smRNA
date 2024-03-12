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
Add option to do alignment with STAR, might be useful for analyzing riboseq data and tRNA alignment
Add option to do alignment with bowtie2
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

if (params.outdir)   { ; } else { exit 1, 'Results path not specified!' }
params.results = "${params.outdir}/${params.primary_aligner}"

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

// Preprocessing
include { DESIGN_INPUT } from '../../modules/misc/main.nf'
include { TRIMGALORE_INPUT } from '../../modules/misc/main.nf'
include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { ARTIFACTS_FILTER } from '../../modules/filter/main.nf'
include { TRIM_POLYA } from '../../modules/filter/main.nf'
include { TRIM_UMI } from '../../modules/filter/main.nf'
include { HARDTRIM } from '../../modules/filter/main.nf'

// Bowtie
include { BOWTIE_INDEX_GENOME } from '../../modules/bowtie/main.nf'
include { BOWTIE_INDEX_JUNCTION } from '../../modules/bowtie/main.nf'
include { BOWTIE_ALIGN_GENOME } from '../../modules/bowtie/main.nf'
include { BOWTIE_ALIGN_JUNCTION } from '../../modules/bowtie/main.nf'
include { COMBINE_GENOME_JUNC_BED } from '../../modules/bowtie/main.nf'
include { REMOVE_CONTAMINANTS_BED } from '../../modules/bowtie/main.nf'
include { PROCESS_ALIGNMENT } from '../../modules/bowtie/main.nf'
include { BOWTIE_ALIGNMENT_MASTER_TABLE } from '../../modules/bowtie/main.nf'
include { REMOVE_CONTAMINANT } from '../../modules/bowtie/main.nf'
include { RBIND_ALIGNMENT_LOG } from '../../modules/misc/main.nf'
include { GET_TOTAL_READS } from '../../modules/misc/main.nf'

/*
// Bowtie2
include { BOWTIE2_INDEX_GENOME } from '../../modules/bowtie2/main.nf'
include { BOWTIE2_INDEX_JUNCTION } from '../../modules/bowtie2/main.nf'
include { BOWTIE2_ALIGN_GENOME } from '../../modules/bowtie2/main.nf'
include { BOWTIE2_ALIGN_JUNCTION } from '../../modules/bowtie2/main.nf'
*/

// Feature counting 
include { COUNT_FEATURES } from '../../modules/counter/main.nf'
include { TAILOR_COUNT_FEATURES } from '../../modules/counter/main.nf'

// Transcript alignment
include { TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { TAILOR_TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { RBIND_COUNTS } from '../../modules/misc/main.nf'
include { INDEX_TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { INDEX_TRANSCRIPTS_TAILOR } from '../../modules/transcripts/main.nf'

// Read aggregation
include { MASTER_TABLE } from '../../modules/misc/main.nf'

// Differential gene expression analysis
include { DGE } from '../../modules/misc/main.nf'

// Tailor
include { TAILOR_INDEX } from '../../modules/tailor/main.nf'
include { TAILOR_INDEX_JUNCTIONS } from '../../modules/tailor/main.nf'
include { TAILOR_ALIGN } from '../../modules/tailor/main.nf'
include { TAILOR_ALIGN_JUNCTIONS } from '../../modules/tailor/main.nf'
include { TAILOR_PROCESS_ALIGNMENT } from '../../modules/tailor/main.nf'
include { TAILOR_RUN } from '../../modules/tailor/main.nf'

// Spike-in normalization
include { ALIGN_SPIKEIN } from '../../modules/spikein/main.nf'
include { NORMALIZE_SPIKEIN } from '../../modules/spikein/main.nf'

// merge bigwigs
include { MERGE_BW } from '../../modules/deeptools/main.nf'

// Align using custom tRNA pipeline
include { INDEX_TRNA } from '../../modules/tRNA_alignment/main.nf'
include { ALIGN_TRNA } from '../../modules/tRNA_alignment/main.nf'

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

    if (params.hardtrim) {

        HARDTRIM( fasta_ch )

        fasta_ch = HARDTRIM.out.fasta

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
    Check bowtie genome index is made
    ////////////////////////////////////////////////////////////////////
    */
    genome_fasta = file("${params.genome}")
    genome_name = "${genome_fasta.baseName}"
    bowtie_index = "${params.index}/bowtie/${genome_name}"
    bowtie_exists = file(bowtie_index).exists()

    if ( params.align_junction ) {

        junction_fasta = file("${params.junctions}")
        junction_name = "${junction_fasta.baseName}"
        bowtie_junction_index = "${params.index}/bowtie/${junction_name}"
        bowtie_junction_exists = file(bowtie_junction_index).exists()

        tailor_junction_index = "${params.index}/tailor/${junction_name}"
        tailor_junction_exists = file(tailor_junction_index).exists()

    }

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

    if (params.hardtrim) {

        HARDTRIM( fasta_ch )
        fasta_ch = HARDTRIM.out.fasta

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
     * Alignment using bowtie or tailor
     */
    if ( params.align_genome ) {

        /*
        * Bowtie 
        */
        if ( params.primary_aligner == 'bowtie' ) {

            // Check to make bowtie genome index exists
            if ( bowtie_exists ){
                index_ch = "${bowtie_index}/${genome_name}"
                chrom_sizes = "${bowtie_index}/${genome_name}_chrom_sizes"

            } else {
                BOWTIE_INDEX_GENOME( 
                    params.genome, 
                    "${bowtie_index}/${genome_name}", 
                    "${bowtie_index}/${genome_name}_chrom_sizes", 
                    genome_name)

                index_ch = BOWTIE_INDEX_GENOME.out.bowtie_index
                chrom_sizes = BOWTIE_INDEX_GENOME.out.bowtie_chrom_sizes
            }

            // Align genome
            BOWTIE_ALIGN_GENOME( index_ch, fasta_xc_ch )

            // option to align to junctions
            if ( params.align_junction ) {

                // check to see if junction index exists
                if ( bowtie_junction_exists ){

                    junc_index_ch = "${bowtie_junction_index}/${junction_name}"
                
                } else {
                    BOWTIE_INDEX_JUNCTION( 
                        params.junctions, 
                        "${bowtie_junction_index}/${junction_name}", 
                        junction_name)

                    junc_index_ch = BOWTIE_INDEX_JUNCTION.out.bowtie_index
                }
                
                // align to junctions
                BOWTIE_ALIGN_JUNCTION(junc_index_ch, BOWTIE_ALIGN_GENOME.out.unmapped_fa)

                // combine bed files
                combine_bed_ch = BOWTIE_ALIGN_GENOME
                        .out
                        .genome_aligned_bed
                        .join( BOWTIE_ALIGN_JUNCTION.out.junction_aligned_bed )
                    
                COMBINE_GENOME_JUNC_BED( combine_bed_ch )

                combined_bed_ch = COMBINE_GENOME_JUNC_BED
                    .out
                    .bed
                    .join( fasta_xc_ch )
            } else {
                // if not align junctions take bowtie align genome as the combined bed channel
                combined_bed_ch = BOWTIE_ALIGN_GENOME
                    .out
                    .genome_aligned_bed
                    .join ( fasta_xc_ch )

            }

        /*
        * Tailor
        */
        } else {
                
            // check to see if tailor index exists 
            if ( tailor_exists ){

                index_ch = "${tailor_index}/${genome_name}"
                chrom_sizes = "${tailor_index}/${genome_name}_chrom_sizes"

            } else {

                TAILOR_INDEX( 
                    params.genome, 
                    "${tailor_index}/${genome_name}", 
                    "${tailor_index}/${genome_name}_chrom_sizes", 
                    genome_name)

                index_ch = TAILOR_INDEX.out.tailor_index
                chrom_sizes = TAILOR_INDEX.out.chrom_sizes

            }

            // align reads using tailor
            TAILOR_ALIGN( index_ch, fasta_xc_ch )

            // if align junction
            if ( params.align_junction ) {

                // check to see if junction exists
                if ( tailor_junction_exists ){

                    junc_index_ch = "${tailor_junction_index}/${junction_name}"
                
                } else {

                    /// index junctions with tailor
                    TAILOR_INDEX_JUNCTIONS( 
                        params.junctions, 
                        "${bowtie_junction_index}/${junction_name}", 
                        junction_name)

                    junc_index_ch = TAILOR_INDEX_JUNCTIONS.out.tailor_index
                }

                // tailor align junctions
                TAILOR_ALIGN_JUNCTIONS( junc_index_ch, TAILOR_ALIGN.out.unmapped_fq_ch )

                combine_bed_ch = TAILOR_ALIGN
                    .out
                    .ntm_ch
                    .join( TAILOR_ALIGN_JUNCTIONS.out.ntm_ch )
                    
                COMBINE_GENOME_JUNC_BED( combine_bed_ch )

                combined_bed_ch = COMBINE_GENOME_JUNC_BED
                    .out
                    .bed
                    .join( fasta_xc_ch )

            } else {

                combined_bed_ch = TAILOR_ALIGN
                    .out
                    .ntm_ch
                    .join ( fasta_xc_ch )
            }
        }

        /*
        * Filter bed
        */
        if ( params.filter_bed ) {
            
            REMOVE_CONTAMINANTS_BED( combined_bed_ch, params.filter_bed )
            process_bed_ch = REMOVE_CONTAMINANTS_BED.out.bed
            
        } else {

            process_bed_ch = combined_bed_ch

        }

        /* 
        * Process alignment 
        */
        if ( params.primary_aligner == 'bowtie' ) {

            PROCESS_ALIGNMENT( process_bed_ch, chrom_sizes )

        } else {

            TAILOR_PROCESS_ALIGNMENT( process_bed_ch, chrom_sizes )

        }

        /* 
        * combine alignment logs 
        */
        if ( params.primary_aligner == 'bowtie' ){

            alignment_logs_ch = PROCESS_ALIGNMENT
                        .out
                        .alignment_logs
                        .collect()

        } else {

            alignment_logs_ch = TAILOR_PROCESS_ALIGNMENT
                        .out
                        .alignment_logs
                        .collect()

        }

        RBIND_ALIGNMENT_LOG( params.outprefix, alignment_logs_ch )

        /*
        * Counter
        */
        if ( params.features ) {

            if ( params.primary_aligner == 'bowtie' ) {

                counter_input_ch = PROCESS_ALIGNMENT.out.channel_ntm

                COUNT_FEATURES(
                    params.features,
                    params.bed, 
                    counter_input_ch
                    )

                norm_consts_ch = COUNT_FEATURES.out.normalization_constants
            
                genomic_counts_ch = COUNT_FEATURES.out.counts

            } else {

                counter_input_ch = TAILOR_PROCESS_ALIGNMENT.out.channel_ntm

                TAILOR_COUNT_FEATURES(
                    params.features,
                    params.bed, 
                    counter_input_ch,
                    params.genome
                    )

                norm_consts_ch = TAILOR_COUNT_FEATURES.out.normalization_constants
            
                genomic_counts_ch = TAILOR_COUNT_FEATURES.out.counts

            }

        } else {

            norm_consts_ch = PROCESS_ALIGNMENT.out.normalization_constants
            genomic_counts_ch = Channel.empty()
        }
    } else {

        GET_TOTAL_READS(fasta_xc_ch)
        
        norm_consts_ch = GET_TOTAL_READS.out.normalization_constants
        
        genomic_counts_ch = Channel.empty()

    }

    /*
     * define post-alignment input ch with sample, fasta, norm_consts
     */
    post_alignment_input_ch = fasta_xc_ch.join(norm_consts_ch)

    /*
    * Transcripts
    */
    if (params.transcripts) {

        if ( params.primary_aligner == 'bowtie' ){

            INDEX_TRANSCRIPTS( params.transcripts )

            transcripts_input_ch = fasta_xc_ch.join(norm_consts_ch)

            TRANSCRIPTS(params.transcripts, post_alignment_input_ch, INDEX_TRANSCRIPTS.out.transcript_index_path_ch)

            transcript_counts_ch = TRANSCRIPTS.out.counts
        
        } else {

            INDEX_TRANSCRIPTS_TAILOR( params.transcripts )

            transcripts_input_ch = fasta_xc_ch.join(norm_consts_ch)

            TAILOR_TRANSCRIPTS(params.transcripts, post_alignment_input_ch, INDEX_TRANSCRIPTS_TAILOR.out.transcript_index_path_ch)

            transcript_counts_ch = TAILOR_TRANSCRIPTS.out.counts

        }
        
    } else {

        transcript_counts_ch = Channel.empty()
    
    }

    /*
     * tRNA alignment pipelines
     */
    if ( params.tRNA_pipeline ) {

        tRNA_fasta = file("${params.tRNA_reference}")

        tRNA_name = "${tRNA_fasta.baseName}"

        tRNA_index = "${params.index}/transcripts/${tRNA_name}"

        INDEX_TRNA( params.tRNA_reference, tRNA_index )

        ALIGN_TRNA( INDEX_TRNA.out.tRNA_sequence_index_path_ch, post_alignment_input_ch)

        tRNA_counts_ch = ALIGN_TRNA.out.counts

    } else {

        tRNA_counts_ch = Channel.empty()

    }



    /*
    * Combine genomic and transcripts counts
    */
    if (params.features || params.transcripts || params.tRNA_pipeline ) {

        counts_ch = genomic_counts_ch
            .mix( transcript_counts_ch )
            .mix( tRNA_counts_ch )
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
     * will only run if tailor == true, features and bed == 2 and primary aligner is bowtie
     */
    if (params.tailor && params.features && params.bed && params.align_genome && params.primary_aligner == 'bowtie'){

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
        if (params.align_junction) {
            tailor_input_ch = BOWTIE_ALIGN_JUNCTION
                .out
                .tailor_input
                .join(norm_consts_ch)    
        } else {
            tailor_input_ch = BOWTIE_ALIGN_GENOME
                .out
                .tailor_input
                .join(norm_consts_ch)
        }

        TAILOR_RUN(
            params.genome, 
            tailor_index_ch, 
            params.features, 
            params.bed, 
            tailor_input_ch
            )
    }
}
