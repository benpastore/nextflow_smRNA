#!/usr/bin/env python3

import sys
import os

def process(fasta) :

    val = 0
    with open(fasta, 'r') as f : 
        for line in f : 
            if line.startswith(">") :
                val += float( line.strip().split(":")[1] )
    f.close()

    print(f"total\t{val}")

def main() : 
    
    process(sys.argv[1])

if __name__ == "__main__" : 

    main()
