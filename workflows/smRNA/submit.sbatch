#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf -profile cluster --design /fs/ess/PAS1473/deep_sequencing/rde12_rsd6/rde12_rsd_design.txt --results /fs/ess/PAS1473/deep_sequencing/rde12_rsd6/results_v3 --outprefix results --features /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/feature/features.txt --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.slim.bed --contaminant false --artifacts_filter false --tailor true --transcripts /fs/ess/PAS1473/deep_sequencing/rde12_rsd6/transcripts.txt --dge /fs/ess/PAS1473/deep_sequencing/rde12_rsd6/dge.txt --multimap 1000 --min_length 15 --format fasta -resume



