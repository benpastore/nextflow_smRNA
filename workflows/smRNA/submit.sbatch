#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf -profile cluster --design /fs/ess/PAS1473/deep_sequencing/germano_riboseq/design_riboseq.txt --results /fs/ess/PAS1473/deep_sequencing/germano_riboseq/results --outprefix results --features /fs/ess/PAS1473/deep_sequencing/germano_riboseq/features --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.slim.bed --contaminant /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.xk.fa --tailor false --multimap 1000 --min_length 10 --format fastq --trimming true --adapter TGGAATTCTCGGGTGCCAAGG --trim_umi true --umi_p5 4 --umi_p3 4 --rpkm true



