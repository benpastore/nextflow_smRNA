process TRANSCRIPTS {

    errorStrategy 'retry'
    maxRetries 3
    label 'high' 

    publishDir "$params.results/transcripts", mode : 'copy', pattern : "*.tsv"

    input : 
        val transcripts
        tuple val(sampleID), val(fasta), val(normalization)
    
    output : 
        tuple val(sampleID), path("*transcripts.counts.tsv"), emit : counts
        path("*.transcripts.counts.tsv"), emit : master_table_input
        path("*.transcripts.bed.tsv"), emit : transcript_bed_ch

    script : 
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename $fasta .fa)
    out=\$name.v${mismatch}.m${multimap}.transcripts

    python3 ${params.bin}/align_transcripts.py -f ${fasta} \\
        -t ${params.transcripts} \\
        -o \$out \\
        -n ${normalization} \\
        -v ${mismatch} \\
        -m ${multimap}
    """
}