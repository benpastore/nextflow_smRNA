#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf -profile cluster --design /fs/ess/PAS1473/deep_sequencing/DISL2_published_seq/CLIP_smRNA_pipeline.txt --results /fs/ess/PAS1473/deep_sequencing/DISL2_published_seq/disl2_clip_smRNA --outprefix results --features /fs/ess/PAS1473/deep_sequencing/DISL2_published_seq/features_GRCh38.txt --bed /fs/ess/PCON0160/ben/genomes/GRCh38/gencode.v42.annotation.GRCh38tRNA.bed --genome /fs/ess/PCON0160/ben/genomes/GRCh38/GRCh38.primary_assembly.genome.fa --align_junction false --contaminant false --tailor true --transcripts false --dge false --multimap 1000 --mismatch 0 --min_length 10 --max_length 100 --format fastq -resume



