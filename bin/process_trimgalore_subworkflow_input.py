#!/usr/bin/env python3

import argparse
import os

def process(fasta) : 

    info = [ i for i in my_string.split(",") ]

    for i in info : 
        print(f"{os.path.basename(i).split('.')[0]},{i}")

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-input", type = str, required = False)

    return parser.parse_args()

def main() : 

    args = get_args() 
    process(args.input)

if __name__ == "__main__" : 

    main()