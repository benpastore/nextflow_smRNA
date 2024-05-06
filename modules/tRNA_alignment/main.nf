process INDEX_TRNA {

    publishDir "$params.index/transcripts", mode : 'copy', pattern : "*bwt"
    publishDir "$params.index/transcripts", mode : 'copy', pattern : "*sizes"


    input :
        val tRNA_sequences
        val tRNA_index_path

    output : 
        path("*")
        val("${tRNA_index_path}"), emit : tRNA_sequence_index_path_ch
        path("*.chrom_sizes"), emit : tRNA_sizes
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    python3 ${params.bin}/index_tRNA.py --reference ${tRNA_sequences}
    """
}

process ALIGN_TRNA {

    errorStrategy 'ignore'
    label 'low'

    publishDir "$params.results/tRNA_alignment/bed", mode : 'copy', pattern : "*.bed.tsv"
    publishDir "$params.results/tRNA_alignment/counts", mode : 'copy', pattern : "*.counts.tsv"
    publishDir "$params.results/tRNA_alignment/bw", mode : 'copy', pattern : "*.bw"

    input :
        val tRNA_index
        tuple val(sampleID), val(fasta), val(normalization)

    output : 
        path("*.tsv")
        tuple val(sampleID), path("*.counts.tsv"), emit : counts
        path("*bw")
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    # align tRNA with 0,1,2,3 mismatches
    id=\$(basename ${fasta} .fa)

    python3 ${params.bin}/align_tRNAs.py \\
        --index ${tRNA_index} \\
        --fasta ${fasta} \\
        --basename \$id \\
        --norm_consts ${normalization}

    # convert bed file to bw
    cat \$id.tRNA.aligned.bed.tsv |\\
        grep -v "^gene" |\\
        awk '{OFS="\\t"; print \$1,\$2,\$3,\$4,\$10,\$6}' |\\
        sort-bed - > tmp.sorted

    bedops --partition tmp.sorted |\\
        bedmap --echo --sum --delim '\\t' - tmp.sorted |\\
        awk -F '\\t' '{OFS="\\t"; print \$1,\$2,\$3,1*\$4}' > tmp.bg

    ${params.bin}/bedGraphToBigWig tmp.bg ${tRNA_index}.chrom_sizes \$id.tRNA.bw

    """
}