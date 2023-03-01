#!/bin/bash

while getopts ":fastq:trimmed:fasta:adapter:quality:error:minL:maxL:basename:threads" i; do
        case "${i}" in
        fastq)
                fastq=$OPTARG
        ;;
        trimmed)
                t_fastq=$OPTARG
        ;;
        fasta)
                fasta=$OPTARG
        ;;
        adapter)
                adapter=$OPTARG
        ;;
        quality)
                quality=$OPTARG
        ;;
        error)
                error=$OPTARG
        ;;
        minL)
                minL=$OPTARG
        ;;
        maxL)
                maxL=$OPTARG
        ;;
        basename)
                basename=$OPTARG
        ;;
        threads)
                threads=$OPTARG
        ;;
        esac
done

trim_galore -j ${threads} ${adapter} -q 30 -e 0.1 --gzip ${min_length} ${max_length} --fastqc --basename ${sampleID} \$fq
