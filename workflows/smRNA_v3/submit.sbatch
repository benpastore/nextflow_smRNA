#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

module load java/12.0.2

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA_v3/main.nf -profile cluster --design /fs/ess/PCON0160/ben/projects/2025_csr1_cleavage/csr_wago_IP.txt --outdir /fs/ess/PCON0160/ben/projects/2025_csr1_cleavage/00AA_AGO_cleavage --outprefix results --features /fs/ess/PCON0160/ben/projects/2025_csr1_cleavage/features.txt --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.slim.bed --contaminant false --tailor true --transcripts /fs/ess/PCON0160/ben/projects/2025_csr1_cleavage/transcripts.txt --dge false --multimap 10 --min_length 10 --mismatch 2 --format fastq



