process DOWNLOAD_SRR {

    label 'fastqdump'

    tag "${sample}_download"

    publishDir "$params.results/fastq", mode : 'copy', pattern : '*fastq.gz'

    input :
        tuple val(sra), val(sample)
    
    output : 
        path("*")
        tuple val(sample), path("*.fastq.gz"), emit : fastq
    
    script : 
    """
    #!/bin/bash

    wget -q https://sra-pub-run-odp.s3.amazonaws.com/sra/${sra}/${sra} 

    fasterq-dump --threads ${task.cpus} --outdir \$PWD --split-files ${sra} --quiet

    #rename "${sra}" "${sample}" *
    rename 's/${sra}/${sample}/' *
    
    gzip *fastq
        
    """
}