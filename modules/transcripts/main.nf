process INDEX_TRANSCRIPTS {

    label 'medium'

    publishDir "$params.index/transcripts", mode : 'copy', pattern : "*bwt"
    publishDir "$params.index/transcripts", mode : 'copy', pattern : "*sizes"


    input :
        val transcripts

    output : 
        path("*")
        val("$params.index/transcripts"), emit : transcript_index_path_ch
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    python3 ${params.bin}/index_transcripts.py --transcripts ${transcripts}

    """
}

process INDEX_TRANSCRIPTS_TAILOR {

    label 'medium'

    publishDir "$params.index/tailor_transcripts", mode : 'copy', pattern : "*"

    input :
        val transcripts

    output : 
        path("*")
        val("$params.index/tailor_transcripts"), emit : transcript_index_path_ch
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    python3 ${params.bin}/index_transcripts_tailor.py --transcripts ${transcripts} --tailor ${params.bin}/tailor_v11

    """
}

process TRANSCRIPTS {

    errorStrategy 'ignore'
    //maxRetries 3
    label 'high' 

    publishDir "$params.results/transcripts/bed", mode : 'copy', pattern : "*bed.tsv"
    //publishDir "$params.results/transcripts/bw", mode : 'copy', pattern : "*bw"
    publishDir "$params.results/transcripts/counts", mode : 'copy', pattern : "*counts.tsv"

    input : 
        val transcripts
        tuple val(sampleID), val(fasta), val(normalization)
        val transcript_index_path
    
    output : 
        tuple val(sampleID), path("*transcripts.counts.tsv"), emit : counts
        path("*.transcripts.counts.tsv"), emit : master_table_input
        path("*.transcripts.bed.tsv"), emit : transcript_bed_ch
        //path("*.bw")

    script : 
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    rpkm_command = params.rpkm ? "-rpkm" : ''
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
        -m ${multimap} \\
        -idx ${transcript_index_path} \\
        ${rpkm_command}

    #cat \$out.bed.tsv |\\
    #    grep -v "^gene" |\\
    #    awk '{OFS="\\t"; print \$1,\$2,\$3,\$4,\$8,\$6}' |\\
    #    sort-bed - > tmp.sorted

    #bedops --partition tmp.sorted |\\
    #    bedmap --echo --sum --delim '\\t' - tmp.sorted |\\
    #    awk -F '\\t' '{OFS="\\t"; print \$1,\$2,\$3,1*\$4}' > tmp.bg

    #${params.bin}/bedGraphToBigWig tmp.bg chrom_sizes_combined \$out.bw
        
    """
}

process TAILOR_TRANSCRIPTS {

    label 'high' 

    publishDir "$params.results/tailor_transcripts/files", mode : 'copy', pattern : "*.tsv"

    input : 
        val transcripts
        tuple val(sampleID), val(fasta), val(normalization)
        val transcript_index_path
    
    output : 
        path("*transcripts.bed.tsv"), emit : bed
        path("*.transcripts.counts.tsv"), emit : master_table_input
        path("*.transcripts.bed.tsv"), emit : transcript_bed_ch

    script : 
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    tailor_mismatch = params.tailor_mismatch ? "-allow_mismatch" : ''
    rpkm_command = params.rpkm ? "-rpkm" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    python3 ${params.bin}/uniq_fasta_to_uniq_fastq.py ${fasta} > fastq

    name=\$(basename $fasta .fa)
    out=\$name.v${mismatch}.m${multimap}.transcripts

    python3 ${params.bin}/align_transcripts_tailor.py \\
        -f fastq \\
        -t ${transcripts} \\
        -o \$out \\
        -n ${normalization} \\
        -v ${mismatch} \\
        -m ${multimap} --tailor "${params.bin}" \\
        ${tailor_mismatch} \\
        -idx ${transcript_index_path} \\
        ${rpkm_command}
        
    """
}

process RBIND_TAILOR_TRANSCRIPTS {

    label 'low'

    publishDir "$params.results/tailor_transcripts", mode : 'copy', pattern : "*.tsv"

    input :
        val name
        val files
    
    output : 
        path("*.tsv"), emit : tables

    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    
    python3 ${params.bin}/rbind_tailor_transcripts.py -f "${files}" -o ${name}.counts.tsv

    """
}