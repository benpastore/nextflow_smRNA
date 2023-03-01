#!/usr/bin/env python3 

import pandas as pd 
import numpy as np 
import argparse 
import time 
import logging
from utilities import timing

@timing
def parse_normalization_constants(spike_ins) : 

    norm_factors = {}
    for spk in spike_ins :
        with open(spk, 'r') as f : 
            for line in f : 
                info = line.strip().split("\t")
                sample = info[0]
                count = float(info[1])
                norm_factors[sample] = count
        f.close()
    
    return norm_factors

@timing
def normalize(counts, spike_ins) : 

    # read in counts df 
    df = pd.read_csv(counts, sep = '\t')

    # parse spike in normalization factors 
    norm_factors = parse_normalization_constants(spike_ins)

    df_non = df.select_dtypes(exclude=[np.float]).reset_index(drop = True)
    df_num = df.select_dtypes(include=[np.float]).reset_index(drop = True)

    cols = df_num.columns

    for col in cols : 
        df_non[col] = df_num[col]/norm_factors[col]
    
    return df_non

@timing
def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument(
        "-c", 
        "--counts", 
        type = str, 
        required = True,
        help="tab delimited table of unnormalized counts")
    
    parser.add_argument(
        "-s", 
        "--spike_ins", 
        type = str, 
        required = True,
        help="list of files with spike in normalization factors OR directory with spike in normalization files, number of files should be equal to the number of samples in data table -c")
    
    parser.add_argument(
        "-d", 
        required = False,
        action = 'store_true',
        help="indicator on whether or not -s is representative of a directory. Default : False/None")

    parser.add_argument(
        "-o", 
        "--output", 
        type = str, 
        required = True,
        help="output")

    return parser.parse_args()
    
@timing
def main() : 
    
    logging.basicConfig(level=logging.DEBUG)
    
    args = get_args()

    if args.d : 
        files = glob.glob(f"{args.spike_ins}")
    else : 
        if "[" in args.spike_ins : 
            files = args.spike_ins.replace("[","").replace("]","").split(", ")
        else : 
            files = args.spike_ins.replace("[","").replace("]","").split(",")
            
    res = normalize(args.counts, files)
    res.to_csv(args.output, sep = '\t', header = True, index = False)


if __name__ == "__main__" : 

    main()