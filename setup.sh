#!/bin/bash

#0. Install singularity
which singularity || echo "Please install singularity"

#1. Install nextflow
which nextflow || curl -s https://get.nextflow.io | bash

# add nextflow to $PATH or add nextflow to directory already in $PATH or make symbolic link (ln -s from to) 