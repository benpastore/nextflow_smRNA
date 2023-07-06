
process ARTIFACTS_FILTER {

    label 'low'

    publishDir "$params.results/trim_galore/x_artifacts", mode : 'copy', pattern : "*.xartifacts.fa"

    input :
        tuple val(sampleID), val(fasta)
    
    output : 
        tuple val(sampleID), path("${sampleID}*.xartifacts.fa"), emit : fasta
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${fasta} .fa)

    ou_fasta=\$name.xartifacts.fa

    python3 ${params.bin}/artifacts_filter.py -i ${fasta} -o \$ou_fasta -t 0.65
    """
}

process TRIM_POLYA {

    label 'low'

    publishDir "$params.results/trim_galore/x_polyA", mode : 'copy', pattern : "*.xpolyA.fa"

    input :
        tuple val(sampleID), val(fasta)
    
    output : 
        tuple val(sampleID), path("${sampleID}*.xpolyA.fa"), emit : fasta
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${fasta} .fa)

    ou_fasta=\$name.xpolyA.fa

    python3 ${params.bin}/trim_3p_A.py -i ${fasta} -o \$ou_fasta
    """
}

process TRIM_UMI {

    label 'low'

    publishDir "$params.results/trim_galore/x_umi", mode : 'copy', pattern : "*xumi*"

    input :
        tuple val(sampleID), val(fasta)
    
    output : 
        tuple val(sampleID), path("${sampleID}*.xumi.uniq.fa"), emit : fasta
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${fasta} .fa)

    ou_fasta=\$name.xumi.uniq.fa

    python3 ${params.bin}/trim_3p_A.py -i ${fasta} -o \$ou_fasta -p5 ${params.umi_p5} -p3 ${params.umi_p3}
    """
}