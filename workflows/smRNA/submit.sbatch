#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf -profile cluster --design /fs/ess/PCON0160/ben/projects/disl2cde1/design/smRNA_design.txt --results /fs/ess/PCON0160/ben/projects/disl2cde1/design/../smRNA_results_rerun --outprefix piRNA_tailing_trimming --features /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/feature/features.txt --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.bed --contaminant /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.xk.fa --tailor true --dge /fs/ess/PCON0160/ben/projects/disl2cde1/design/smRNA_comparisons.txt --transcripts /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/transcripts/transcripts.txt



