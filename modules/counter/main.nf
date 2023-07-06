process COUNT_FEATURES {

    label 'low'

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

    script : 
    normalize_command = params.features_norm ? "-n ${params.features_norm}" : ''
<<<<<<< HEAD
    rpkm_command = params.rpkm ? "-rpkm" : ''
=======
>>>>>>> 3f1103a9195e2e904318679d1ea7e24fe48260ca
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${alignment} .ntm)

    time python3 ${params.bin}/count.py \\
        -f ${features} \\
        -a ${reference_annotation} \\
        -i ${alignment} \\
        -o \$name \\
<<<<<<< HEAD
        ${normalize_command} ${rpkm_command}
=======
        ${normalize_command}
>>>>>>> 3f1103a9195e2e904318679d1ea7e24fe48260ca
        
    """

}

