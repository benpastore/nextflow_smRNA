#!/bin/bash

tailor="/fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/bin/tailor_v11"
sam2bed="/fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/bin/tailor_sam_to_bed"
filter_tailor=""
reference=$1
fastq=$2

reference_name=$(basename ${reference} .fa)
sampleID=$(basename ${fastq} .fq)

# build tailor reference
[ -d "./tailor_index" ] || mkdir -p "./tailor_index"
$tailor build -i ${reference} -p ${reference_name} -f 
mv ${reference_name}* "./tailor_index"

# tailor align
sam=${sampleID}.tailor.sam
ntm=${sampleID}.tailor.ntm
${tailor} map \
        -i ${fastq} \
        -p "./tailor_index/${reference_name}" \
        -n 8 \
    tee ${sam} | \
    ${sam2bed} | \
    awk '($5 <= 1)' | \
    awk -F'\t' -v OFS="\t" '{
        split($4,a,":")
        print $1,$2,$3,a[1],a[2]/$5,$6,$8,$5
    }' > ${ntm}

cat ${ntm} | awl '($7=="*")' > ${reference_name}.perfect.ntm
cat ${ntm} | awl '($7!="*")' > ${reference_name}.tailed.ntm

python3 
