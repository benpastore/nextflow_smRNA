process TAILOR_INDEX {

    label 'low'
    publishDir "$params.tailor_index_path", mode : 'copy'

    input :
        val fasta
    
    output :
        path("*")
        val("${params.tailor_index}"), emit : tailor_index

    script : 
    f = file("${fasta}")
    name = "${f.baseName}"
    """
    #!/bin/bash

    source activate smrnaseq

    ${params.bin}/tailor_v11 build -i ${fasta} -p ${name}
    """

}

process TAILOR_MAP {

    label 'low'

    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*tailor.bed"
    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*.bai"
    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*.sam"
    publishDir "$params.results/tailor/counts", mode : 'copy', pattern : "*bed.tsv"

    input :
        val genome
        val talor_index
        val features
        val reference_annotation
        tuple val(sampleID), val(fastq), val(normalization_constants)

    output :
        path("*tailor.bed")
        path("*tsv")


    script :
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

    name=\$(basename $fastq .fq)
    sam=\$name.sam
    bam=\$name.bam
    bed=\$name.aligned.v0.m1.tailor.bed
    counts=\$name.aligned.v0.m1.tailor
    
    ${params.bin}/tailor_v11 map \\
        -i ${fastq} \\
        -p ${talor_index} \\
        -n ${task.cpus} \\
        2> tailor.log | \\
    tee \$sam | \\
    ${params.bin}/tailor_sam_to_bed | \\
    awk -v num=\$nTag -F'\\t' -v OFS="\\t" '{
        if (\$8!="*")
        {
            split(\$4,a,":")

            print \$1,\$2,\$3,a[1],a[2]/\$5,\$6,\$8,\$5

        }
    }' > \$bed
    
    # sort sam file -> bam file
    samtools sort -m 1G -@ ${task.cpus} -o \$bam \$sam

    # index bam file -> bam.bai 
    samtools index -@ ${task.cpus} \$bam

    # Count reads
    python3 ${params.bin}/tailor_count.py \\
        -i \$bed \\
        -a ${reference_annotation} \\
        -f ${features} \\
        -n ${normalization_constants} \\
        -o \$counts \\
        -g ${genome}

    """
}

process TAILOR_FILTER {

    label 'medium'

    publishDir "$params.results/tailor", mode : 'copy', pattern : "*tsv"

    input :
        val directory
        val condition

    output :
        path("*tailor.bed")
        path("*tsv")        

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    # MUST map sense and be unique mapper

    name=\$(basename $fastq .fq)
    bed=\$name.aligned.v0.m1.tailor.bed
    counts=\$name.aligned.v0.m1.tailor
    
    ${params.bin}/tailor_v11 map \\
        -i ${fastq} \\
        -p ${talor_index} \\
        -n ${task.cpus} \\
        2> tailor.log | \\
    tee aligned.sam | \\
    ${params.bin}/tailor_sam_to_bed | \\
    awk -v num=\$nTag -F'\\t' -v OFS="\\t" '{
        if (\$8!="*")
        {
            split(\$4,a,":")

            print \$1,\$2,\$3,a[1],a[2]/\$5,\$6,\$8,\$5

        }
    }' > \$bed

    # Count reads
    python3 ${params.bin}/tailor_count.py \\
        -i \$bed \\
        -a ${reference_annotation} \\
        -f ${features} \\
        -n ${normalization_constants} \\
        -o \$counts \\
        -g ${genome}

    """
}