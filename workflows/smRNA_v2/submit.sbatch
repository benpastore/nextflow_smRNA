#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA_v2/main.nf -profile cluster --design /fs/ess/PCON0160/ben/projects/2024_MID_domain_mutants/MID_DOMAIN_DESIGN.txt --outdir /fs/ess/PCON0160/ben/projects/2024_MID_domain_mutants/results --outprefix results --features /fs/ess/PCON0160/ben/projects/2024_MID_domain_mutants/features.txt --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.slim.bed --contaminant /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.xk.fa --tailor false --transcripts /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/transcripts/transcripts.txt --dge false --multimap 1000 --min_length 12 --format fastq --primary_aligner bowtie --mismatch 0 --format fastq



