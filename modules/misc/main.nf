
// process design input
process DESIGN_INPUT {

    label 'local'

    container 'docker://benpasto/smrnaseq:latest'

    publishDir "$params.results/samples", mode : 'copy', pattern : "*csv"

    input : 
        val(design)

    output : 
        path("fastq.csv"), emit : fastq_ch
        path("replicates.csv"), emit : condition_ch

    script : 
    """
    #!/bin/bash
        
    source activate smrnaseq
    
    python3 ${params.bin}/process_design_input_smRNA.py -input ${design}
    
    """

}

// Combine Alignment Logs
process RBIND_ALIGNMENT_LOG {

    label 'low'

    publishDir "$params.results", mode : 'copy', pattern : "*.log"

    input : 
        val(projectName)
        val(files)
    
    output : 
        path("*.log")

    script : 
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    """
    #!/bin/bash

    source activate smrnaseq
    
    python3 ${params.bin}/rbind_tables.py -f "${files}" -o ${projectName}.alignment.v${mismatch}.m${multimap}.log
    """
}

process RBIND_COUNTS {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*.tsv"

    input : 
        tuple val(sampleID), val(files)
    
    output : 
        path("*.tsv"), emit : tables

    script : 
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    """
    #!/bin/bash

    source activate smrnaseq
    
    python3 ${params.bin}/rbind_tables.py -f "${files}" -o ${sampleID}.v${mismatch}.m${multimap}.counts.tsv

    """
}

// Make Counts Master Table
process MASTER_TABLE {

    label 'low'

    publishDir "$params.results/master_tables", mode : 'copy', pattern : "*.tsv"

    input :
        val project_name
        val counts

    output : 
        path("*count*tsv"), emit : tables
        path("*.count.tsv"), emit : unnormalized_master_table_ch
    
    script :
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    counts=${project_name}.aligned.v${mismatch}.m${multimap}

    time python3 ${params.bin}/make_master.py -f "${counts}" -o \$counts
    
    time python3 ${params.bin}/geometric_normalization.py -i \$counts.count.tsv -o \$counts

    """
}

process DGE {

    label 'medium'

    publishDir "$params.results/DGE", mode : 'copy', pattern : "*.tsv"

    input : 
        tuple val(x), val(y), val(samples_x), val(samples_y), val(counts)
    
    output : 
        path("*.tsv")
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    time python3 ${params.bin}/compare/compare.py \\
        -nx ${x} \\
        -ny ${y} \\
        -x "${samples_x}" \\
        -y "${samples_y}" \\
        -c ${counts} \\
        -o \$PWD 
    """
}

// TRIMGALORE Input
process TRIMGALORE_INPUT {

    label 'local'

    publishDir "$params.results/samples", mode : 'copy', pattern : "*csv"

    input : 
        val(design)

    output : 
        path("trimgalore_subworkflow_input.csv"), emit : fastq_ch

    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    
    python3 ${params.bin}/process_trimgalore_subworkflow_input.py -input ${design} > trimgalore_subworkflow_input.csv
    """
}

process GET_TOTAL_READS {

    label 'low'

    publishDir "$params.results/total_reads", mode : 'copy', pattern : "*.tsv"

    input :
        tuple val(sampleID), val(fasta)
    
    output : 
        path("${sampleID}*_total_reads.tsv")
        tuple val(sampleID), path("${sampleID}_total_reads.tsv"), emit : normalization_constants

    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    fa=${fasta}
    total_reads=${sampleID}_total_reads.tsv

    python3 ${params.bin}/sum_reads.py \$fa > \$total_reads

    """
}