

// Trim Galore
process TRIM_GALORE {

    label 'low'

    publishDir "$params.results/trim_galore/fastq", mode : 'copy', pattern : "*_trimmed.fq.gz"
    publishDir "$params.results/trim_galore/logs", mode : 'copy', pattern : "*_trimming_report.txt"
    publishDir "$params.results/trim_galore/collapsed", mode : 'copy', pattern : "*.trimmed.uniq.fa"

    input : 
        val min_length
        val max_length
        tuple val(sampleID), val(fastq)
    
    output : 
        path("${sampleID}*_trimming_report.txt")
        tuple val(sampleID), path("${sampleID}_trimmed.fq.gz")
        tuple val(sampleID), path("${sampleID}.trimmed.uniq.fa"), emit : collapsed_fa

    script : 
    """
    #!/bin/bash 

    source activate smrnaseq

    fq=${fastq}
    t_fq=${sampleID}_trimmed.fq.gz
    fa=${sampleID}.trimmed.uniq.fa
    
    trim_galore -j ${task.cpus} ${params.trimgalore.a} -q 30 -e 0.1 --gzip --length ${min_length} --max_length ${max_length} --fastqc --basename ${sampleID} \$fq

    sh ${params.bin}/fastq_to_uniq_fasta.sh \$t_fq \$fa
    
    """
}

// Bowtie Index
process BOWTIE_INDEX {

    label 'low'
    publishDir "$params.bowtie.directory", mode : 'copy'
    
    input :
        val genome
        val juncs

    output :
        path("*.ebwt")
        path("*chrom_sizes")
        path("${name}*"), emit : bowtie_index

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

// Remove Contaminant RNA
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


// Bowtie Align Genome
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
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename $fasta .fa)
    id=\$name.aligned.v${params.bowtie.mismatch}.m${params.bowtie.multimap}
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
        --strata ${params.bowtie.v} ${params.bowtie.m} \\
        -S > aligned.genome.sam 2> genome.log
    
    bowtie \\
        -x ${idx}_juncs \\
        -f \$name.tmp \\
        -p ${task.cpus} \\
        -a \\
        --un \$name.unmapped.uni.fa \\
        --best \\
        --strata ${params.bowtie.v} ${params.bowtie.m} \\
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

// Combine Alignment Logs
process RBIND_ALIGNMENT_LOG {

    label 'low'

    publishDir "$params.results", mode : 'copy', pattern : "*.log"

    input : 
        val(projectName)
        val(files)
    
    output : 
        path("*.log")

    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    
    python3 ${params.bin}/rbind_tables.py -f "${files}" -o ${projectName}.aligned.v${params.bowtie.mismatch}.m${params.bowtie.multimap}.log

    """
}

// Assign Mapped Reads to Genomic Features
process COUNT_FEATURES {

    label 'low'

    publishDir "$params.results/bed", mode : 'copy', pattern : "*.bed.tsv"
    publishDir "$params.results/normalization", mode : 'copy', pattern : "*.normalization.constants.tsv"

    input :
        val features
        val normalize_to
        val reference_annotation
        tuple val(sampleID), val(alignment)
    
    output : 
        tuple val(sampleID), path("*.normalization.constants.tsv"), emit : normalization_constants
        tuple val(sampleID), path("*counts.tsv"), emit : counts
        path("*counts.tsv"), emit : master_table_input
        path("*bed.tsv"), emit : bed_counts

    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename $alignment .ntm)

    time python3 ${params.bin}/count.py \\
        -f ${features} \\
        -a ${reference_annotation} \\
        -i ${alignment} \\
        -o \$name \\
        ${normalize_to}
        
    """

}

// Align Reads to Transcriptome(s)
process TRANSCRIPTS {

    errorStrategy 'retry'
    maxRetries 3
    label 'medium' 

    publishDir "$params.results/transcripts", mode : 'copy', pattern : "*.tsv"

    input : 
        val transcripts
        tuple val(sampleID), val(fasta), val(normalization)
    
    output : 
        tuple val(sampleID), path("*transcripts.counts.tsv"), emit : transcripts_output
        path("*.transcripts.counts.tsv"), emit : master_table_input
        path("*.transcripts.bed.tsv"), emit : transcript_bed_ch

    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename $fasta .fa)
    out=\$name.v${params.transcripts.mismatch}.m${params.transcripts.multimap}.transcripts

    python3 ${params.bin}/align_transcripts.py -f ${fasta} \\
        -t ${params.transcripts.file} \\
        -o \$out \\
        -n ${normalization} \\
        -v ${params.transcripts.mismatch} \\
        -m ${params.transcripts.multimap}
    """
}

// Combine Transcriptome and Genomic Assignments
process RBIND_COUNTS {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*.tsv"

    input : 
        tuple val(sampleID), val(files)
    
    output : 
        path("*.tsv"), emit : tables

    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    
    python3 ${params.bin}/rbind_tables.py -f "${files}" -o ${sampleID}.v${params.bowtie.mismatch}.m${params.bowtie.multimap}.counts.tsv

    """
}

// Make Counts Master Table
process MASTER_TABLE {

    label 'low'

    publishDir "$params.results/master_tables", mode : 'copy', pattern : "*.tsv"

    input :
        val project_name
        val counts

    output : 
        path("*norm*tsv"), emit : tables
        path("*count.tsv")
    
    script :
    """
    #!/bin/bash

    source activate smrnaseq

    counts=${project_name}.aligned.v${params.bowtie.mismatch}.m${params.bowtie.multimap}

    time python3 ${params.bin}/make_master.py -f "${counts}" -o \$counts
    
    time python3 ${params.bin}/geometric_normalization.py -i \$counts.count.tsv -o \$counts

    """
}

// Preform Differential Gene Expression (DGE) Analysis
process DGE {

    label 'medium'

    publishDir "$params.results/DGE", mode : 'copy', pattern : "*.tsv"

    input : 
        val comparisons
        val counts
    
    output : 
        path("*.tsv")
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    time python3 ${params.bin}/compare/compare_samples.py -c ${comparisons} -f ${counts}
    """
}

// Index Genome Using Tailor Index Builder
process TAILOR_INDEX {

    label 'low'
    publishDir "$params.tailor.directory", mode : 'copy'

    input :
        val fasta
    
    output :
        path("*")
        path("${name}*"), emit : tailor_index

    script : 
    f = file("${fasta}")
    name = "${f.baseName}"
    """
    #!/bin/bash

    source activate smrnaseq

    ${params.bin}/tailor_v11 build -i ${fasta} -p ${name}
    """

}

// Map reads to genome with Tailor to identify nontemplated nucleotide additions
process TAILOR_MAP {

    label 'low'

    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*tailor.bed"
    publishDir "$params.results/tailor/counts", mode : 'copy', pattern : "*.tsv"

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

/*
// MiRDeep2. Identify / quantify known and novel miRNAs in deep sequencing data
process MIRDEEP2 {

    publishDir "$params.results/mirdeep", mode : 'copy'

    input : 

    output : 

    script:
    """
    #!/bin/bash

    module load python
    source activate rnaseq_basic

    # re format fasta file to meet miRDeep format
    

    """




}
*/
