#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

module load java/12.0.2

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA_v2/main.nf -profile cluster --design /fs/ess/PCON0160/ben/projects/2024_human_dis3l2/design_ago2_bound_miRNA.txt --outdir /fs/ess/PCON0160/ben/projects/2024_human_dis3l2/ago2_bound_miRNA_results --outprefix results --features /fs/ess/PCON0160/Chin/smRNA-seq/DIS3L2_KO/features.txt --align_junction false --genome /fs/ess/PCON0160/Chin/genomes/GRCh38/GRCh38.primary_assembly.genome.fa --tRNA_reference /fs/ess/PCON0160/Chin/genomes/GRCh38/hg38-mature-tRNAs.modified.fa --bed /fs/ess/PCON0160/Chin/genomes/GRCh38/gencode.v45.basic.annotation.miRBASE.fixed.bed --transcripts /fs/ess/PCON0160/Chin/smRNA-seq/DIS3L2_KO/transcripts.txt --tailor true --multimap 1000 --mismatch 0 --min_length 15 --primary_aligner bowtie --format fastq



