

process TAILOR_MAP {

    label 'low'

    publishDir "$params.results/alignment/tailor/bed", mode : 'copy', pattern : "*tailor*bed"
    publishDir "$params.results/alignment/tailor/bam", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results/alignment/tailor/bam", mode : 'copy', pattern : "*.bai"
    publishDir "$params.results/alignment/tailor/bed", mode : 'copy', pattern : "*bed.tsv"
    publishDir "$params.results/alignment/tailor/counts", mode : 'copy', pattern : "*counts.tsv"

    input :
        val genome
        val talor_index
        val tailor_junction_index
        val features
        val reference_annotation
        tuple val(sampleID), val(fasta)

    output :
        path("*bai")
        path("*bam")
        path("*bed")
        path("*tsv")
        //val (sampleID), path("*.counts.tsv"), emit : tailor_counts_ch


    script :
    normalize_command = params.features_norm ? "-n ${params.features_norm}" : ''
    rpkm_command = params.rpkm ? "-rpkm" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    # Bed file format
    # 1-3 -> chrom, start, end
    # 4 -> seq
    # 5 -> ntm
    # 6 -> strand (all should be +)
    # 7 -> tail
    # 8 -> tail length
    # 9 -> number of locations mapped

    # MUST map sense and be unique mapper

    name=${sampleID}  #\$(basename $fastq .fa)
    sam=\$name.tailor.sam
    bam=\$name.tailor.bam
    bed=\$name.aligned.v0.m1.tailor.bed
    counts=\$name.aligned.v0.m1.tailor
    fastq=\$name.fa
    unmapped=\$name.unmapped.fastq

    # convert fa to fq
    python3 ${params.bin}/uniq_fasta_to_uniq_fastq.py ${fasta} > \$fastq
    
    # analyze 3' nontemplated nucleotides
    ${params.bin}/tailor_v11 map \\
        -i \$fastq \\
        -p ${talor_index} \\
        -n ${task.cpus} \\
        2> tailor.genome.log | \\
    tee \$sam | \\
    ${params.bin}/tailor_sam_to_bed | \\
    awk -v num=\$nTag -F'\\t' -v OFS="\\t" '{
        split(\$4,a,":")
        print \$1,\$2,\$3,a[1],a[2]/\$5,\$6,\$8,\$5
    }' > \$bed

    # sort sam file -> bam file
    samtools sort -m 1G -@ ${task.cpus} -o \$bam \$sam

    # index bam file -> bam.bai 
    samtools index -@ ${task.cpus} \$bam

    # get reads that dont align, align to junction, convert to genome coordinate bed
    python3 ${params.bed}/filter_reads.py -i \$fastq -o \$unmapped -b \$bed
    
    # Count reads
    python3 ${params.bin}/tailor_count_v2.py \\
        -i \$bed \\
        -a ${reference_annotation} \\
        -f ${features} \\
        -o \$counts \\
        -fasta ${genome} \\
        ${normalize_command} \\
        ${rpkm_command}

    ############################################################################
    ############################################################################
    
    """
}