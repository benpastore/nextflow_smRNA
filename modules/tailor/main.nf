process TAILOR_INDEX {

    label 'high'
    
    publishDir "$params.index/tailor", mode : 'copy'

    input :
        val genome
        val index
        val name

    output :
        val(index), emit : tailor_index
        val("*chrom_sizes"), emit : chrom_sizes
        path("${name}/*")
        val("*chrLen")
        val("*chrStart")
        val("*NposLen.z")
        val("*_bwt.bwt")
        val("*_seq.bwt")
        val("*_table.bwt")

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    [ ! -d ${name} ] && mkdir -p ${name}

    # bowtie
    ${params.bin}/tailor_v11 build -i ${genome} -p ${name}/${name}

    samtools faidx ${genome}
    cat ${genome}.fai | cut -f1,2 > ${name}/${name}_chrom_sizes
    rm ${genome}.fai
    """
}

process TAILOR_INDEX_JUNCTIONS {

    label 'low'
    publishDir "$params.index/tailor", mode : 'copy'

    input :
        val fasta
        val index
        val name
    
    output :
        path("${name}/*")
        val(index), emit : tailor_index

    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    [ ! -d ${name} ] && mkdir -p ${name}

    ${params.bin}/tailor_v11 build -i ${fasta} -p ${name}/${name}
    """
}

process TAILOR_ALIGN {

    label 'high'

    publishDir "$params.results/alignment/ntm", mode : 'copy', pattern : "*ntm"
    publishDir "$params.results/alignment/bam", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results/alignment/bam", mode : 'copy', pattern : "*.bai"
    //publishDir "$params.results/alignment/unmapped", mode : 'copy', pattern : "*.unmapped.fastq"
    //publishDir "$params.results/alignment/tailor/bed", mode : 'copy', pattern : "*bed.tsv"
    //publishDir "$params.results/alignment/tailor/counts", mode : 'copy', pattern : "*counts.tsv"

    input :
        val talor_index
        tuple val(sampleID), val(fasta)

    output :
        path("*bai")
        path("*ntm")
        path("*bam")
        //path("*.fastq")
        //path("*bed")
        //val (sampleID), path("*.counts.tsv"), emit : tailor_counts_ch
        tuple val(sampleID), path("*.ntm"), emit : ntm_ch
        tuple val(sampleID), path("*.fq"), emit : unmapped_fq_ch


    script :
    //normalize_command = params.features_norm ? "-n ${params.features_norm}" : ''
    //rpkm_command = params.rpkm ? "-rpkm" : ''
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

    name=${sampleID}  #\$(basename ${fasta} .fa)
    sam=\$name.tailor.sam
    bam=\$name.tailor.bam
    bed=\$name.aligned.v0.m${params.multimap}.tailor.bed
    ntm=\$name.aligned.v0.m${params.multimap}.tailor.ntm
    counts=\$name.aligned.v0.m${params.multimap}.tailor
    fastq=\$name.fq
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
    awk '(\$5 <= ${params.multimap})' | \\
    awk -v num=\$nTag -F'\\t' -v OFS="\\t" '{
        split(\$4,a,":")
        print \$1,\$2,\$3,a[1],a[2]/\$5,\$6,\$8,\$5
    }' > \$ntm

    # sort sam file -> bam file
    samtools sort -m 1G -@ ${task.cpus} -o \$bam \$sam

    # index bam file -> bam.bai 
    samtools index -@ ${task.cpus} \$bam

    # get reads that dont align, align to junction, convert to genome coordinate bed
    #python3 ${params.bin}/filter_reads.py -i \$fastq -o \$unmapped -b \$ntm
    """
}

process TAILOR_ALIGN_JUNCTIONS {

    label 'high'

    publishDir "$params.results/alignment/junction_ntm", mode : 'copy', pattern : "*ntm"
    publishDir "$params.results/alignment/junction_unmapped", mode : 'copy', pattern : "*fastq"

    input :
        val talor_index
        tuple val(sampleID), val(fastq)

    output :
        path("*ntm")
        //val (sampleID), path("*.counts.tsv"), emit : tailor_counts_ch
        tuple val(sampleID), path("*.ntm"), emit : ntm_ch

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    name=${sampleID}
    ntm=\$name.aligned.v0.m1.tailor.junction.ntm

    ${params.bin}/tailor_v11 map \\
        -i ${fastq} \\
        -p ${talor_index} \\
        -n ${task.cpus} \\
        2> tailor.genome.log | \\
    tee tmp.sam | \\
    ${params.bin}/tailor_sam_to_bed | \\
    awk -v num=\$nTag -F'\\t' -v OFS="\\t" '{
        split(\$4,a,":")
        print \$1,\$2,\$3,a[1],a[2]/\$5,\$6,\$8,\$5
    }' > tmp2

    python3 ${params.bin}/tailor_junc_bed2bed.py -i tmp2 -o \$ntm
    """
}

process TAILOR_PROCESS_ALIGNMENT {

    label 'high'
    
    publishDir "$params.results/alignment/rpm", mode : 'copy', pattern : "*.rpm"
    publishDir "$params.results/alignment/bw", mode : 'copy', pattern : "*.bw"
    publishDir "$params.results/alignment/rpkm", mode : 'copy', pattern : "*.rpkm"
    publishDir "$params.results/logs", mode : 'copy', pattern : "*.log"

    input : 
        tuple val(sampleID), path(ntm), path(fasta)
        val(chrom_sizes)

    output : 
        tuple val(sampleID), path("*.depth"), emit : normalization_constants
        tuple val(sampleID), path(ntm), emit : channel_ntm
        path("*aligned*.log"), emit : alignment_logs
        path("*.rpm"), emit : channel_rpm
        path("*.bw")
        path("*.rpkm")


    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${ntm} .ntm)
    id=\$name
    rpm=\$id.rpm
    rpkm=\$id.rpkm
    bw=\$id.bw
    log=\$id.log
    dep=\$id.depth


    # find depth from ntm file
    depth=\$(python3 ${params.bin}/add_columns.py ${ntm} 4)
    echo \$depth > depth

    # normalize ntm file to rpm
    echo -e "chrom\tstart\tend\tseq\tcount_rpm\tstrand\ttail" > header 
    cat ${ntm} | awk -F'\\t' -v OFS='\\t' -v nTag=\$depth '{print \$1,\$2,\$3,\$4,1000000*(\$5/nTag),\$6,\$7}' > rpm.tmp
    cat header rpm.tmp > \$rpm

    # normalize ntm file to rpkm, remove multimappers and non-perfect alignments
    echo -e "chrom\tstart\tend\tseq\tcount_rpkm\tstrand\ttail" > header 
    cat \$rpm | awk -F'\\t' -v OFS='\\t' '{print \$1,\$2,\$3,\$4,\$5/length(\$4),\$6,\$7}' > rpkm.tmp
    cat header rpkm.tmp > \$rpkm

    # make bw file of all alignments
    cat rpm.tmp | awk -F'\\t' '{OFS="\\t"; if (\$2>=\$3) print \$1,\$2,\$3+1,\$4,\$5; else print \$1,\$2,\$3,\$4,\$5 }' | sort-bed - > tmp.sorted
    bedops --partition tmp.sorted | bedmap --echo --sum --delim '\\t' - tmp.sorted | awk -F '\\t' '{OFS="\\t"; print \$1,\$2,\$3,1*\$4}' > tmp.bg
    ${params.bin}/bedGraphToBigWig tmp.bg ${chrom_sizes} \$bw

    # find total reads sequenced
    total_reads=\$(cat ${fasta} | grep "^>" | sed -e 's/.*://g' | ${params.bin}/addCols stdin)
    
    # compute total aligned reads (decimal percent)
    percent_aligned=\$(python3 -c "print(\$depth/\$total_reads)")
    
    # Alignment summary
    echo -e "sample\ttotal_reads\tdepth\tpercent_aligned" > \$log
    echo -e \$id'\t'\$total_reads'\t'\$depth'\t'\$percent_aligned >> \$log

    # Normalization constants 
    echo -e 'total\t'\$depth > \$dep
    """
}

process TAILOR_RUN {

    label 'low'

    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*tailor*bed"
    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*.bai"
    publishDir "$params.results/tailor/alignment", mode : 'copy', pattern : "*.sam"
    publishDir "$params.results/tailor/bed", mode : 'copy', pattern : "*bed.tsv"
    publishDir "$params.results/tailor/counts", mode : 'copy', pattern : "*counts.tsv"

    input :
        val genome
        val talor_index
        val features
        val reference_annotation
        tuple val(sampleID), val(fastq), val(normalization_constants)

    output :
        path("*tailor*bed")
        path("*tsv")
        //val (sampleID), path("*.counts.tsv"), emit : tailor_counts_ch


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

    name=${sampleID}  #\$(basename $fastq .fq)
    sam=\$name.sam
    bam=\$name.bam
    bed=\$name.aligned.v0.m1.tailor.bed
    counts=\$name.aligned.v0.m1.tailor
    
    # analyze 3' nontemplated nucleotides
    ${params.bin}/tailor_v11 map \\
        -i ${fastq} \\
        -p ${talor_index} \\
        -n ${task.cpus} \\
        2> tailor.log | \\
    tee \$sam | \\
    ${params.bin}/tailor_sam_to_bed | \\
    awk '(\$5<=1)' | \\
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

process COMBINE_TAILOR_COUNTS {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*tsv"

    input :
        tuple val(sampleID), val(counts), val(tailed_counts)

    output :
        path("*tsv")

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${count} .tsv)
    outname=\$name.counts.plus.tailed.tsv

    python3 ${params.bin}/combine_counts.py \\
        -i ${counts},${tailed_counts} \\
        -o \$outname
    """

}