#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf -profile cluster --design /fs/ess/PCON0160/ben/projects/antisense_piRNA_v2/anti_piRNA_design.txt --results /fs/ess/PCON0160/ben/projects/antisense_piRNA_v2/results_07_14_2023 --outprefix analysis --features /fs/ess/PCON0160/ben/projects/antisense_piRNA_v2/features.txt --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.slim.bed --contaminant /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.xk.fa --tailor true --transcripts /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/transcripts/transcripts.txt --format fastq --hardtrim true --hardtrim_5p 4 --hardtrim_3p 0 --adapter CTGTAGGCA --multimap 1 --min_length 12



