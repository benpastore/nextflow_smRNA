process MERGE_BW {

    tag "${condition}_merge_bw"

    label 'high'

    publishDir "$params.results/merge_BigWig", mode: 'copy', pattern : "*.bw"

    input : 
        tuple val(sampleIDs), path(bws), val(condition)
        val chrom_sizes

    output : 
        tuple val(condition), path("*.bw")
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    # make array with bigwigs 
    bws=(${bws.join(' ')})

    # get length of array
    N_bw=\${#bws[@]}

    echo \$N_bw

    # merge bigwig(s) --> bedgraph
    bigWigMerge ${bws.join(' ')} ${condition}.tmp

    # divide counts column by N samples to average
    cat ${condition}.tmp | awk -F'\\t' -v OFS='\\t' -v nsamp=\$N_bw '{ print \$1,\$2,\$3,\$4/nsamp }' > ${condition}.bedGraph

    # sort the bed file 
    bedSort ${condition}.bedGraph ${condition}.bedGraph.sorted

    # convert bedGraph --> bigwig
    bedGraphToBigWig ${condition}.bedGraph.sorted ${chrom_sizes} ${condition}.bw
    """

}