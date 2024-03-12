#!/usr/env/bin python3

import argparse
import os

def index(ref, path_to_tailor) :

    index_prefix = os.path.basename(ref).replace(".fasta", "").replace(".fa","").replace(".ref", "")
    cmd = f"{path_to_tailor} build -i {ref} -p {index_prefix}" 
    os.system(cmd)
    
    chrom_sizes = os.path.basename(ref).replace(".fasta", "").replace(".fa","").replace(".ref", "")
    os.system(f"samtools faidx {ref}")
    os.system(f"cut -f1,2 {ref}.fai > {chrom_sizes}.chrom_sizes")

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument(
        "-t", 
        "--transcripts", 
        type = str, 
        required = True,
        help="tab delimited transcripts table")

    parser.add_argument(
        "--tailor", 
        type = str, 
        required = True,
        help="path to tailor software")

    return parser.parse_args()

def main() :

    args = get_args()

    with open(args.transcripts, 'r') as f : 
        for line in f :
            if not line.startswith("Name") : 
                info = line.strip().split("\t")
                reference = info[8]
                print(f"Indexing {os.path.basename(reference)}")

                index(reference, args.tailor)
    f.close()

if __name__ == "__main__" : 

    main()