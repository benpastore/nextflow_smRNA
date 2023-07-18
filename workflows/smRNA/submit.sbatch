#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf -profile cluster --design /fs/ess/PCON0160/ben/projects/dev/disl2_design.txt --results /fs/ess/PCON0160/ben/projects/dev/results_07_12_2023_no_junction --outprefix results --features /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/feature/features.txt --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.slim.bed --contaminant /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.xk.fa --tailor true --transcripts /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/transcripts/transcripts.txt --dge false --multimap 1000 --min_length 12 --format fastq --align_junction false



