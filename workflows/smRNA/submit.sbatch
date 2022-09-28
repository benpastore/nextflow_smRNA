#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf -profile cluster --tailor false --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.exon.protein_coding.window.300.bed --design /fs/ess/PCON0160/ben/projects/disl2cde1/design_smRNA.txt --results /fs/ess/PCON0160/ben/projects/disl2cde1/map_22G_to_window --outprefix DISL2 --run_transcripts false

