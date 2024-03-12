#!/usr/env/bin python3

import argparse
import os

def index(ref) :

    index_prefix = os.path.basename(ref).replace(".fasta", "").replace(".fa", "")
    cmd = f"bowtie-build --quiet {ref} {index_prefix}" 
    os.system(cmd)
    
    chrom_sizes = os.path.basename(ref).replace(".fasta", "").replace(".fa", "").replace(".ref", "")
    os.system(f"samtools faidx {ref}")
    os.system(f"cut -f1,2 {ref}.fai > {chrom_sizes}.chrom_sizes")

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument(
        "--reference", 
        type = str, 
        required = True)

    return parser.parse_args()

def main() :

    args = get_args()
    index(args.reference)

if __name__ == "__main__" : 

    main()