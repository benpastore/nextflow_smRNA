#!/usr/bin/env nextflow

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

/*
////////////////////////////////////////////////////////////////////
set path to bin and index
////////////////////////////////////////////////////////////////////
*/
params.bin = "${params.base}/../../bin"
params.index = "${params.base}/../../index"
if (params.outdir)   { ; } else { exit 1, 'Results path not specified!' }
params.date = new Date().format( 'yyyyMMdd' )
params.results = "${params.outdir}/${params.date}"
params.outprefix = "${params.date}"

/*
////////////////////////////////////////////////////////////////////
Included modules
////////////////////////////////////////////////////////////////////
*/
include { DESIGN_INPUT } from '../../modules/misc/main.nf'
include { TRIMGALORE_INPUT } from '../../modules/misc/main.nf'
include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { ARTIFACTS_FILTER } from '../../modules/filter/main.nf'
include { TRIM_POLYA } from '../../modules/filter/main.nf'
include { TRIM_UMI } from '../../modules/filter/main.nf'
include { HARDTRIM } from '../../modules/filter/main.nf'
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
include { COUNT_FEATURES } from '../../modules/counter/main.nf'
include { TAILOR_COUNT_FEATURES } from '../../modules/counter/main.nf'
include { TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { RBIND_COUNTS } from '../../modules/misc/main.nf'
include { INDEX_TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { MASTER_TABLE } from '../../modules/misc/main.nf'
include { DGE } from '../../modules/misc/main.nf'
include { TAILOR_INDEX } from '../../modules/tailor/main.nf'
include { TAILOR_INDEX_JUNCTIONS } from '../../modules/tailor/main.nf'
include { TAILOR_ALIGN } from '../../modules/tailor/main.nf'
include { TAILOR_ALIGN_JUNCTIONS } from '../../modules/tailor/main.nf'
include { TAILOR_PROCESS_ALIGNMENT } from '../../modules/tailor/main.nf'
include { TAILOR_RUN } from '../../modules/tailor/main.nf'
include { ALIGN_SPIKEIN } from '../../modules/spikein/main.nf'
include { NORMALIZE_SPIKEIN } from '../../modules/spikein/main.nf'
include { MERGE_BW } from '../../modules/deeptools/main.nf'
include { INDEX_TRNA } from '../../modules/tRNA_alignment/main.nf'
include { ALIGN_TRNA } from '../../modules/tRNA_alignment/main.nf'
include { INDEX_TRANSCRIPTS_TAILOR } from '../../modules/transcripts/main.nf' 
include { TAILOR_TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { RBIND_TAILOR_TRANSCRIPTS } from '../../modules/transcripts/main.nf'
include { DOWNLOAD_SRR } from '../../modules/fastqdump/main.nf'

/*
////////////////////////////////////////////////////////////////////
Workflows
////////////////////////////////////////////////////////////////////
*/

workflow download_fastq {

    take : 
        data
    
    main : 
        DOWNLOAD_SRR( data )
    
    emit : 
        fastq = DOWNLOAD_SRR.out.fastq

}

workflow parse_design {

    take : 
        data 
    
    main : 
        DESIGN_INPUT( data )

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
    
    emit : 
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
}

workflow trim_galore {

    take : 
        data 

    main : 
        TRIM_GALORE( data )
    
    emit : 
        trim_galore_output = TRIM_GALORE.out.collapsed_fa

}

workflow hardtrim {

    take : 
        data 
    
    main : 
        HARDTRIM( data )
    
    emit : 
        fasta = HARDTRIM.out.fasta
}

workflow trim_umi {

    take : 
        data 
    
    main : 
        TRIM_UMI( data )
    
    emit : 
        fasta = TRIM_GALORE.out.fasta
    
}

workflow filter_artifacts {

    take : 
        data 
    
    main : 
        ARTIFACTS_FILTER( data )
    
    emit : 
        fasta = ARTIFACTS_FILTER.out.fasta

}

workflow trim_polya {

    take : 
        data 
    
    main : 
        TRIM_POLYA( data )
    
    emit : 
        fasta = TRIM_POLYA.out.fasta

}

workflow remove_contaminant {

    take : 
        data
        contaminant

    main : 
        REMOVE_CONTAMINANT( contaminant, data )
    
    emit : 
        fasta = REMOVE_CONTAMINANT.out.xk_fasta

}

workflow bowtie_index {

    take : 
        data 
        index

    main : 
        genome_fasta = file("${data}")
        genome_name = "${genome_fasta.baseName}"
        bowtie_index = "${index}/bowtie/${genome_name}"
        bowtie_exists = file(bowtie_index).exists()

        if ( bowtie_exists ){
                index_ch = "${bowtie_index}/${genome_name}"
                chrom_sizes = "${bowtie_index}/${genome_name}_chrom_sizes"

        } else {
            BOWTIE_INDEX_GENOME( 
                data, 
                "${bowtie_index}/${genome_name}", 
                "${bowtie_index}/${genome_name}_chrom_sizes", 
                genome_name)

            index_ch = BOWTIE_INDEX_GENOME.out.bowtie_index
            chrom_sizes = BOWTIE_INDEX_GENOME.out.bowtie_chrom_sizes
        }

    emit : 
        index_ch = index_ch
        chrom_sizes = chrom_sizes

}

workflow bowtie_index_junction { 

    take : 
        data
        index

    main : 
        junction_fasta = file("${data}")
        junction_name = "${junction_fasta.baseName}"
        bowtie_junction_index = "${index}/bowtie/${junction_name}"
        bowtie_junction_exists = file(bowtie_junction_index).exists()

        if ( bowtie_junction_exists ){

            junc_index_ch = "${bowtie_junction_index}/${junction_name}"
        
        } else {
            BOWTIE_INDEX_JUNCTION( 
                params.junctions, 
                "${bowtie_junction_index}/${junction_name}", 
                junction_name)

            junc_index_ch = BOWTIE_INDEX_JUNCTION.out.bowtie_index
        }
    
    emit : 
        junc_index_ch = junc_index_ch

}

workflow bowtie_align {

    take : 
        index
        data
    
    main : 
        BOWTIE_ALIGN_GENOME(index, data)
    
    emit : 
        aligned = BOWTIE_ALIGN_GENOME.out.genome_aligned_bed
        unaligned = BOWTIE_ALIGN_GENOME.out.unmapped_fa
        tailor_input = BOWTIE_ALIGN_GENOME.out.tailor_input

}

workflow bowtie_junction_align {

    take : 
        index 
        data 
    
    main : 
        BOWTIE_ALIGN_JUNCTION( index, data )
    
    emit : 
        aligned = BOWTIE_ALIGN_JUNCTION.out.junction_aligned_bed
        tailor_input = BOWTIE_ALIGN_JUNCTION.out.tailor_input

}

workflow combine_genome_junction_bed {

    take : 
        data 
    
    main : 
        COMBINE_GENOME_JUNC_BED( data )
    
    emit : 
        combined = COMBINE_GENOME_JUNC_BED.out.bed

}

workflow filter_contaminant_bed { 

    take : 
        data
        bed 
    
    main : 
        REMOVE_CONTAMINANTS_BED( data, bed )
    
    emit : 
        bed = REMOVE_CONTAMINANTS_BED.out.bed

}

workflow process_bowtie_alignment { 

    take : 
        data 
        chrom_sizes 
    
    main : 
        PROCESS_ALIGNMENT(data, chrom_sizes)
    
    emit : 
        bed = PROCESS_ALIGNMENT.out.channel_ntm
        logs = PROCESS_ALIGNMENT.out.alignment_logs
        norm = PROCESS_ALIGNMENT.out.normalization_constants

}

workflow make_alignment_info_table {

    take : 
        name
        data

    
    main : 
        RBIND_ALIGNMENT_LOG(name, data )

}

workflow count_features {

    take : 
        features
        annotation
        bed
    
    main : 
        COUNT_FEATURES(features, annotation, bed)
    
    emit : 
        norm_consts = COUNT_FEATURES.out.normalization_constants
        counts = COUNT_FEATURES.out.counts

}

workflow get_total_reads {
    
    take : 
        data 
    
    main : 
        GET_TOTAL_READS(data)
    
    emit : 
        norm_consts = GET_TOTAL_READS.out.normalization_constants

}

workflow bowtie_index_transcripts { 

    take : 
        data
    
    main : 
        INDEX_TRANSCRIPTS( data )
    
    emit : 
        transcript_index = INDEX_TRANSCRIPTS.out.transcript_index_path_ch

}

workflow bowtie_align_transcripts {

    take : 
        transcript_table
        transcript_input
        transcript_index
    
    main : 
        TRANSCRIPTS( transcript_table, transcript_input, transcript_index )

    emit : 
        bed = TRANSCRIPTS.out.counts 

}


workflow index_tRNA_ref { 

    take : 
        index
        tRNA_reference
    
    main : 
        tRNA_fasta = file("${tRNA_reference}")
        tRNA_name = "${tRNA_fasta.baseName}"
        tRNA_index = "${index}/transcripts/${tRNA_name}"
        INDEX_TRNA( tRNA_reference, tRNA_index )
    
    emit : 
        index = INDEX_TRNA.out.tRNA_sequence_index_path_ch
    
}

workflow align_tRNA {

    take : 
        index
        data 
    
    main : 
        ALIGN_TRNA( index, data )
    
    emit : 
        counts = ALIGN_TRNA.out.counts
}

workflow rbind_counts {

    take : 
        data 
    
    main : 
        RBIND_COUNTS( data )
    
    emit : 
        counts = RBIND_COUNTS.out.tables
}

workflow master_table { 

    take : 
        data
        name
    
    main : 
        RBIND_COUNTS( data )
        MASTER_TABLE( name, RBIND_COUNTS.out.tables.collect() )

    emit : 
        master_tables = MASTER_TABLE.out.tables.flatten()
        unnormalized_counts = MASTER_TABLE.out.unnormalized_master_table_ch

}

workflow align_spikein {

    take : 
        data 
        spikeins 
    
    main : 
        ALIGN_SPIKEIN(spikeins, data)
    
    emit : 
        spikein_quant_ch = ALIGN_SPIKEIN.out.spikein_quant_ch

}

workflow quantify_spikein {

    take :
        counts 
        spikeins 
        name
    
    main : 
        NORMALIZE_SPIKEIN( counts, spikeins, name )
    
    emit : 
        counts= NORMALIZE_SPIKEIN.out.counts_normalized_spikein_ch

}

workflow differential_expression_analysis {

    take : 
        data 
    
    main : 
        DGE( data )

}

workflow tailor_index_transcripts { 

    take : 
        data
    
    main : 
        INDEX_TRANSCRIPTS_TAILOR( data )
    
    emit : 
        index = INDEX_TRANSCRIPTS_TAILOR.out.transcript_index_path_ch

}

workflow tailor_align_transcripts { 

    take : 
        transcripts
        info
        index

    main : 
        TAILOR_TRANSCRIPTS( transcripts, info, index )

    emit : 
        counts = TAILOR_TRANSCRIPTS.out.counts

}

workflow tailor_index {

    take : 
        data 
        index
    
    main : 
        genome_fasta = file("${data}")
        genome_name = "${genome_fasta.baseName}"
        tailor_index = "${index}/tailor/${genome_name}"
        tailor_exists = file(tailor_index).exists()

        if ( tailor_exists ){
                tailor_index_ch = "${tailor_index}/${genome_name}"
        } else {
            TAILOR_INDEX( data,  "${tailor_index}/${genome_name}", genome_name )
            tailor_index_ch = TAILOR_INDEX.out.tailor_index
        }
    
    emit : 
        index = tailor_index_ch
}

workflow tailor_align {

    take : 
        genome
        index
        features
        bed
        input


    main : 
        TAILOR_RUN(genome, index, features, bed, input)

}

workflow rbind_tailor_transcripts { 

    take : 
        name 
        data
    
    main : 
        RBIND_TAILOR_TRANSCRIPTS( name, data )
    
    emit : 
        table = RBIND_TAILOR_TRANSCRIPTS.out.tables

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


    // parse design 
    parse_design( params.design ) 

    // trimming and qc 
    if ( params.format == 'fastq' & params.trimming ) {
        trim_galore( parse_design.out.reads_ch )
        fastas = trim_galore.out.trim_galore_output
    } else {
        fastas = parse_design.out.reads_ch
    }

    if (params.hardtrim) { 
        hardtrim( fastas )
        fastas = hardtrim.out.fasta
    }
    
    if (params.trim_umi) {
        trim_umi(fastas)
        fastas = trim_umi.out.fasta
    }

    if (params.filter_artifacts) {
        filter_artifacts( fastas )
        fastas = filter_artifacts.out.fasta
    }

    if (params.trim_polyA) {
        trim_polya( fastas )
        fastas = trim_polya.out.fasta
    }

    if (params.contaminant) {
        remove_contaminant( params.contaminant, fastas )
    }


    // bowtie index and alignment 
    if ( params.align_genome ) {
        bowtie_index(params.genome, params.index)
        bowtie_align( bowtie_index.out.index_ch, fastas)

        if (params.align_junction) {
            bowtie_index_junction(params.junctions, params.index)
            bowtie_junction_align(bowtie_index_junction.out.junc_index_ch, bowtie_align.out.unaligned)
            combine_genome_junction_bed( bowtie_align.out.aligned.join(bowtie_junction_align.out.aligned ) )
        }

        if ( params.filter_bed ) {
            bed = filter_contaminant_bed( combine_genome_junction_bed.out.combined, params.filter_bed )
        } else {
            bed = bowtie_align.out.aligned
        }

        process_bowtie_alignment(bed.join(fastas), bowtie_index.out.chrom_sizes)

        make_alignment_info_table( params.outprefix, process_bowtie_alignment.out.logs.collect() )

        // count features 
        if (params.features) {
            count_features( params.features, params.bed, process_bowtie_alignment.out.bed )
            genomic_counts_ch = count_features.out.counts
            norm_constants = count_features.out.norm_consts
        } else { 
            genomic_counts_ch = Channel.empty()
            norm_constants = process_bowtie_alignment.out.norm
        }

    } else {
        get_total_reads(fastas)
        norm_constants = get_total_reads.out.norm_consts
        genomic_counts_ch = Channel.empty()
    }

    if (params.transcripts) {
        bowtie_index_transcripts( params.transcripts )
        bowtie_align_transcripts( params.transcripts,  fastas.join(norm_constants), bowtie_index_transcripts.out.transcript_index )
        transcript_counts = bowtie_align_transcripts.out.bed
    } else {
        transcript_counts = Channel.empty()
    }

    if (params.tRNA_pipeline) {
        index_tRNA_ref(params.index, params.tRNA_reference)
        align_tRNA(index_tRNA_ref.out.index, fastas.join(norm_constants))
        tRNA_counts = align_tRNA.out.counts
    } else {
        tRNA_counts = Channel.empty()
    }

    if (params.features || params.transcripts || params.tRNA_pipeline) {

        counts = genomic_counts_ch
            .mix(transcript_counts)
            .mix(tRNA_counts)
            .groupTuple()
            .unique()

        //rbind_counts( counts )
        master_table( counts, params.outprefix )
    }

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
        
        dge_input_ch = comparisons_ch.combine( master_table.out.master_tables )
        differential_expression_analysis( dge_input_ch )

    }

    if (params.tailor && params.features && params.bed && params.align_genome) {

        tailor_index(params.genome, params.index)

        if (params.align_junction) {
            tailor_input_ch = bowtie_junction_align
                .out
                .tailor_input
                .join(norm_constants)    
        } else {
            tailor_input_ch = bowtie_align
                .out
                .tailor_input
                .join(norm_constants)
        }

        tailor_align(params.genome, tailor_index.out.index, params.features, params.bed, tailor_input_ch)

    }

    if (params.tailor_transcripts) {

        tailor_index_transcripts( params.tailor_transcripts )

        tailor_align_transcripts( params.tailor_transcripts, fastas.join(norm_constants), tailor_index_transcripts.out.index )
        
        rbind_tailor_transcripts(params.outprefix, tailor_align_transcripts.out.counts.collect())

    }

}

