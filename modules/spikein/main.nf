process INDEX_SPIKEIN {

    publishDir "$params.index/spikein", mode : 'copy', pattern : "*bwt"
    publishDir "$params.index/spikein", mode : 'copy', pattern : "*sizes"

    input : 
        val spikein

    output : 
        path("*")
        val("$params.index/spikein"), emit : spikein_index_path_ch
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    python3 ${params.bin}/index_transcripts.py --transcripts ${spikein}

    """
}

process ALIGN_SPIKEIN {
    
    label 'low'
    
    publishDir "$params.results/spikein", mode : 'copy', pattern : "*spikein*"

    input : 
        val idx
        tuple val(sampleID), path(fasta)
    
    output : 
        path("*spikein"), emit : spikein_quant_ch
        path("*spikein.bed")
        path("*spikein.count")

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
        -m 1 \\
        -v 0 \\
        -f ${fasta} \\
        --un \${name}.xc.fa \\
        -S > \${name}.aligned.spikein.sam 

    sam2bed < \${name}.aligned.spikein.sam | \\
        awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' | \\
        awk '(\$6=="+")' | \\
        awk -F'\\t' -v OFS='\\t' '{ print \$1,\$2,\$3,\$4,\$5,\$6 }' > \${name}.aligned.spikein.bed

    
    # get sum per spike-in
    cat \${name}.aligned.spikein.bed | awk -F'\\t' -v OFS='\\t' '{ print \$1,\$2,\$3,\$4,\$5,\$6 }' | sort -k1,1 | bedtools groupby -g 1 -c 5 -o sum > \${name}.aligned.spikein.count
    
    # get mean of all spike-in
    norm_factor=\$(cat \${name}.aligned.spikein.count | awk -F'\\t' -v OFS='\\t' '{ print ".",\$2}' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o mean | cut -f2)

    # write to file for normalization
    echo -e ${sampleID}'\t'\$norm_factor > \${name}.aligned.spikein

    ###########################################################################################
    #norm_factor=\$(cat \${name}.aligned.spikein.bed | cut -f4 | ${params.bin}/addCols stdin)##
    ###########################################################################################
    """
}

process NORMALIZE_SPIKEIN {
    
    label 'low'
    
    publishDir "$params.results/master_tables", mode : 'copy', pattern : "*count_*.tsv"

    input :
        val unnormalized_counts
        val spikeins
        val spikein_name
    
    output : 
        path("*count_${spikein_name}.tsv"), emit : counts_normalized_spikein_ch
        
    script :
    """
    #!/bin/bash
    
    source activate smrnaseq

    name=\$(basename ${unnormalized_counts} .tsv)

    python3 ${params.bin}/normalize_spikein.py -c ${unnormalized_counts} -s "${spikeins}" -o \$name.count_${spikein_name}.tsv

    """

}