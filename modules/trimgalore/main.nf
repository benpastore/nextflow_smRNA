process TRIM_GALORE {

    label 'low'

    publishDir "$params.results/trim_galore/fastq", mode : 'copy', pattern : "*_trimmed.fq.gz"
    publishDir "$params.results/trim_galore/logs", mode : 'copy', pattern : "*_trimming_report.txt"
    publishDir "$params.results/trim_galore/collapsed", mode : 'copy', pattern : "*.trimmed.uniq.fa"

    input :
        tuple val(sampleID), val(fastq)
    
    output : 
        path("${sampleID}*_trimming_report.txt")
        tuple val(sampleID), path("${sampleID}_trimmed.fq.gz")
        tuple val(sampleID), path("${sampleID}.trimmed.uniq.fa"), emit : collapsed_fa

    script : 
    adapter = params.adapter ? "-a ${params.adapter}" : ''
    min_length = params.min_length ? "--length ${params.min_length}" : ''
    max_length = params.max_length ? "--max_length ${params.max_length}" : ''
    """
    #!/bin/bash 

    source activate smrnaseq

    fq=${fastq}
    t_fq=${sampleID}_trimmed.fq.gz
    fa=${sampleID}.trimmed.uniq.fa
    
    trim_galore -j ${task.cpus} ${adapter} -q 30 -e 0.1 --gzip ${min_length} ${max_length} --fastqc --basename ${sampleID} \$fq

    sh ${params.bin}/fastq_to_uniq_fasta.sh \$t_fq \$fa
    
    """
}