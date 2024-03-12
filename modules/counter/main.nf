process COUNT_FEATURES {

    label 'high'

    publishDir "$params.results/bed", mode : 'copy', pattern : "*.bed.tsv"
    publishDir "$params.results/normalization", mode : 'copy', pattern : "*.normalization.constants.tsv"

    input :
        val features
        val reference_annotation
        tuple val(sampleID), val(alignment)
    
    output : 
        tuple val(sampleID), path("*.normalization.constants.tsv"), emit : normalization_constants
        tuple val(sampleID), path("*counts.tsv"), emit : counts
        path("*counts.tsv"), emit : master_table_input
        path("*bed.tsv"), emit : bed_counts
        tuple val(sampleID), path("*counts.tsv"), emit : counts_ch

    script : 
    normalize_command = params.features_norm ? "-n ${params.features_norm}" : ''
    rpkm_command = params.rpkm ? "-rpkm" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${alignment} .ntm)

    time python3 ${params.bin}/count.py \\
        -f ${features} \\
        -a ${reference_annotation} \\
        -i ${alignment} \\
        -o \$name \\
        ${normalize_command} ${rpkm_command}
        
    """

}

process TAILOR_COUNT_FEATURES {

    label 'high'

    publishDir "$params.results/bed", mode : 'copy', pattern : "*.bed.tsv"
    publishDir "$params.results/normalization", mode : 'copy', pattern : "*.normalization.constants.tsv"

    input :
        val features
        val reference_annotation
        tuple val(sampleID), val(alignment)
        val genome
    
    output : 
        tuple val(sampleID), path("*.normalization.constants.tsv"), emit : normalization_constants
        tuple val(sampleID), path("*counts.tsv"), emit : counts
        path("*counts.tsv"), emit : master_table_input
        path("*bed.tsv"), emit : bed_counts
        tuple val(sampleID), path("*counts.tsv"), emit : counts_ch

    script : 
    normalize_command = params.features_norm ? "-n ${params.features_norm}" : ''
    rpkm_command = params.rpkm ? "-rpkm" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${alignment} .ntm)

    python3 ${params.bin}/tailor_count_v2.py \\
        -i ${alignment} \\
        -a ${reference_annotation} \\
        -f ${features} \\
        -o \$name \\
        -fasta ${genome} \\
        ${normalize_command} \\
        ${rpkm_command}
        
    """

}
