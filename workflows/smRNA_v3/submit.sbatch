#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

module load java/12.0.2

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_smRNA/workflows/smRNA_v3/main.nf -profile cluster --design /fs/ess/PAS1531/ben/2024_oct_cleavage/design_sovic.txt --outdir /fs/ess/PAS1531/ben/2024_oct_cleavage/results_sovic --format fastq --align_genome false --max_length 150 --min_length 15



