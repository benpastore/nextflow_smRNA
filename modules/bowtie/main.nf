process BOWTIE_INDEX {

    label 'low'
    publishDir "$params.bowtie_index_path", mode : 'copy'
    
    input :
        val genome
        val juncs

    output :
        path("*.ebwt")
        path("*chrom_sizes")
        val("${params.bowtie_index}"), emit : bowtie_index

    script : 
    f = file("${genome}")
    name = "${f.baseName}"
    """
    #!/bin/bash

    source activate smrnaseq

    # bowtie
    bowtie-build --quiet ${genome} ${name} --threads ${task.cpus}
    bowtie-build --quiet ${juncs} ${name}_juncs --threads ${task.cpus}

    samtools faidx ${genome}
    cat ${genome}.fai | cut -f1,2 > ${name}_chrom_sizes

    """
}

process REMOVE_CONTAMINANT {
    
    label 'low'
    
    publishDir "$params.results/trim_galore/collapsed_xc", mode : 'copy', pattern : "*.xc.fa"
    publishDir "$params.results/trim_galore/xc_aligned", mode : 'copy', pattern : "*.sam"

    input : 
        val idx 
        tuple val(sampleID), path(fasta)
    
    output : 
        tuple val(sampleID), path("*.xc.fa"), emit : xk_fasta
        path("*.sam")
        
    script :
    """
    #!/bin/bash
    
    source activate smrnaseq 
    
    name=\$(basename ${fasta} .fa)
    
    #build index 
    bowtie-build --quiet ${idx} reference --threads ${task.cpus}
    
    # align 
    bowtie \\
        -x reference \\
        -f ${fasta} \\
        --un \${name}.xc.fa \\
        -S > \${name}.aligned.contaminant.sam
    """

}

process BOWTIE_ALIGN_GENOME {

    label 'low'

    publishDir "$params.results/alignment", mode : 'copy', pattern : "*.ntm"
    publishDir "$params.results/alignment", mode : 'copy', pattern : "*.bai"
    publishDir "$params.results/alignment", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results/logs", mode : 'copy', pattern : "*.log"

    input : 
        val idx
        tuple val(sampleID), path(fasta)

    output : 
        tuple val(sampleID), path("*.unmapped.v0.fq"), emit : tailor_input
        tuple val(sampleID), path("*.ntm"), emit : bowtie_alignment
        tuple val(sampleID), path("*.depth"), emit : normalization_constants
        path("*aligned*.log"), emit : alignment_logs
        path("*.sorted.bam")
        path("*.sorted.bam.bai")

    script : 
    mismatch_command = params.mismatch || params.mismatch == 0 ? "-v ${params.mismatch}" : ''
    multimap_command = params.multimap ? "-m ${params.multimap}" : ''
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${fasta} .fa)
    id=\$name.aligned.v${mismatch}.m${multimap}
    bam=\$id.genome.sorted.bam
    bai=\$id.genome.sorted.bam.bai
    bed=\$id.bed
    ntm=\$id.ntm
    log=\$id.log
    dep=\$id.depth

    #############################################################################
    # Do Alignment with -a --best --strata, get all aligned reads in best stratum
    #############################################################################
    bowtie \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -a \\
        --un \$name.tmp \\
        --best \\
        --strata \\
        ${mismatch_command} \\
        ${multimap_command} \\
        -S > aligned.genome.sam 2> genome.log
    
    bowtie \\
        -x ${idx}_juncs \\
        -f \$name.tmp \\
        -p ${task.cpus} \\
        -a \\
        --un \$name.unmapped.uni.fa \\
        --best \\
        --strata \\
        ${mismatch_command} \\
        ${multimap_command} \\
        -S > aligned.junc.sam 2> junc.log

    ##############################################################################
    # Do Alignment with -v 0, get reads that do not map perfectly, send to Tailor
    ##############################################################################
    bowtie \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -v 0 \\
        -a \\
        --un unmapped.genome.v0.tmp \\
        --best \\
        --strata \\
        -S > mapped.v0.sam

    bowtie \\
        -x ${idx}_juncs \\
        -f unmapped.genome.v0.tmp \\
        -p ${task.cpus} \\
        -v 0 \\
        -a \\
        --un unmapped.genome.junc.v0.tmp \\
        --best \\
        --strata \\
        -S > mapped.junc.v0.sam
    
    python3 ${params.bin}/uniq_fasta_to_uniq_fastq.py unmapped.genome.junc.v0.tmp > \$name.unmapped.v0.fq

    ########################################################################
    # Process alignment
    ########################################################################
    
    # sort sam file -> bam file
    samtools sort -m 1G -@ ${task.cpus} -o \$bam aligned.genome.sam

    # index bam file -> bam.bai 
    samtools index -@ ${task.cpus} \$bam
    
    # convert bam -> bed
    bam2bed < \$bam | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > aligned.genome.sorted.bed

    # process alignment to junciton. bed -> junction bed (split the read across the junction)
    sam2bed < aligned.junc.sam  | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > tmp
    python3 ${params.bin}/junc_bed2bed.py -i tmp -o aligned.junc.bed

    # combine bed and junction bed -> bed
    cat aligned.genome.sorted.bed aligned.junc.bed > \$bed

    # bed -> ntm (split read count / number of locations mapped)
    python3 ${params.bin}/bed_to_ntm.py \$bed \$ntm

    # find depth from ntm file
    depth=\$(python3 ${params.bin}/add_columns.py \$ntm 4)
    echo \$depth > depth

    # find total reads sequenced
    total_reads=\$(cat ${fasta} | grep "^>" | sed -e 's/.*://g' | ${params.bin}/addCols stdin)
    
    # compute total aligned reads (decimal percent)
    percent_aligned=\$(python3 -c "print(\$depth/\$total_reads)")
    
    # Alignment summary
    echo -e "sample\ttotal_reads\tdepth\tpercent_aligned" > \$log
    echo -e \$id'\t'\$total_reads'\t'\$depth'\t'\$percent_aligned >> \$log

    # Normalization constants 
    echo -e 'total\t'\$depth > \$dep

    # fin
    """
}
