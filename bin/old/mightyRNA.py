#!/usr/bin/python3

import os
import argparse

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-p", "--paramsfile", type = str, required = True)
    required.add_argument("-r", "--resume", action = 'store_true', required = False )
    
    return parser.parse_args()

def main() : 

    # get args
    args = get_args()

    # get directories
    dir = os.path.dirname(os.path.realpath(__file__))

    # define main
    main = os.path.join(dir, "main.nf")
    params = args.paramsfile

    # run the pipeline 
    if args.resume : 
        cmd = f"nextflow {main} -params-file {params} -resume"
        os.system(cmd)
    else : 
        cmd = f"nextflow {main} -params-file {params}"
        os.system(cmd)
    

if __name__ == "__main__" : 

    main()
