#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /users/PAS1473/benpasto1/ben/pipelines/nextflow_smRNA/workflows/smRNA/main.nf \
  --features /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/feature/features.txt \
  --tailor false \
  -with-dag dag.png \
  -profile cluster \
  --design /users/PAS1473/benpasto1/ben/projects/disl2cde1/design_smRNA.txt --results /users/PAS1473/benpasto1/ben/projects/disl2cde1/map22G_window --bed /fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.exon.protein_coding.window.300.bed --outprefix disl2_analysis --dge /fs/ess/PCON0160/ben/projects/disl2cde1/comparisons.txt



