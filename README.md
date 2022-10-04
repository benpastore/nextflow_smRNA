# Nextflow smRNA Analysis pipeline

Repo for smRNA pipelines designed using Nextflow

## Installation

1. Clone the repository

```
git clone https://github.com/benpastore/nextflow_smRNA.git
```

2. Install nextflow if necessary.

```
cd nextflow_smRNA

sh setup.sh
```

## Running the pipeline

A basic pipeline run could look like this:
```
Arguments in [ ] are optional

nextflow workflows/smRNA/main.nf
  --design design_file.txt (see example)
  --results /path/to/results/directory
  -profile cluster (To run on HPC otherwise set to standard) \
  [ --mismatch 0 (allowed mismatches for alignment, default = 0 ]
  [ --multimapper 1000 (allowed multimappers for alignment, default = 1000 ]
  [ --features features_sheet.tsv (see ./features/features.txt for example, used in feature counting, default = false) ]
  [ --bed features.bed (IF --features specified must supply bed file with genomic locations of genes..., see example it is a modified bed file), default = false
  [ --features_norm (same format as features.txt, except this file is used to compute normalization constants)
  [ --transcripts transcripts_sheet.tsv (see ./transcripts/transcripts.txt for example, used to align reads to transcriptome, default = false) ]
  [ --dge true (differential gene expression analysis comparisons, format is condition1 vs. condition2 as tab delimmited, conditions must match the design file conditions, default = false) ]
  [ --tailor false (decide to run tailor pipeline to detect nontemplated nucleotide additions at 3' of RNA, default = false, if true will use features.txt to count reads, mutually inclusive with --features --bed) ]

```
