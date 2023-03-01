#!/usr/bin/env python3

import os 
import sys
import argparse
import time
import logging
import pandas as pd
import numpy as np
from scipy import stats
from compare import Compare

def timing(f):
    """
    Helper function for timing other functions
    Parameters
    ----------
    f : function
    Returns
    -------
    function
        new function wrap with timer and logging 
    """
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        logging.debug('{:s} function took {:.10f} s'.format(f.__name__, (time2-time1)))
        return ret
    
    return wrap
    
@timing
def comp(reps, comparisons, counts_file, pseudocount, outdir) :
    
    replicates = {}
    with open(reps, 'r') as f : 
        for line in f :
            if not line.startswith("simple_name") : 
                info = line.strip().split(",")
                sample, condition = info[0], info[1]

                if not condition in replicates.keys() : 
                    replicates[condition] = [sample]
                else :
                    replicates[condition].append(sample)
                print(replicates)
    f.close() 
    
    print(replicates)
    
    with open(comparisons, 'r') as f : 
        for line in f : 
            if not line.startswith("X") : 
                info = line.strip().split("\t")
                x_val = info[0]
                y_val = info[1]
                
                x_samples = replicates[x_val]
                y_samples = replicates[y_val]
                
                print(f"Comparing {x_val} vs {y_val}...")
                
                @timing
                def run_compare(counts_file, x_samples, y_samples, x_val, y_val, outdir, pseudocount) : 
                    curr = Compare(
                        counts = counts_file, 
                        x = x_samples,
                        y = y_samples,
                        x_name = x_val,
                        y_name = y_val,
                        outdir = outdir, 
                        identifier = None,
                        pseudo = pseudocount)
                    curr.compare()
                    
                run_compare(counts_file, x_samples, y_samples, x_val, y_val, outdir, pseudocount)
                
    f.close() 

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-c", "--comparisons", type = str, required = True, help = f"Comparisons table")
    required.add_argument("-r", "--replicates", type = str, required = True, help = f"Replicates table")
    required.add_argument("-f", "--counts", type = str, required=True, help = "Counts table")
    required.add_argument("-o", "--outdir", type = str, required=False)
    required.add_argument("-pseudo", type = float, required = False)

    return parser.parse_args()
 
@timing
def main() :
    
    logging.basicConfig(level=logging.DEBUG)

    args = get_args()

    comp(
        args.replicates, 
        args.comparisons, 
        args.counts, 
        args.pseudo, 
        args.outdir
    )

if __name__ == "__main__" : 

    main()
    
    

'''
    #x = Compare_multi(args.comparisons, args.file, args.output, args.pseudo)
    #x.compare_multi()
class Compare_multi() :

    def __init__(self, input, file, output, pseudo) : 

        self._input = input
        self._file = file
        self._output = output
        self._pseudo = pseudo
        
    def compare_multi(self) :
        
        with open(self._input, 'r') as f : 
            for line in f : 
                params = parse(line)
                print(f"Comparing {params[0]} vs {params[1]}...")
                curr = Compare(
                    counts=self._file, 
                    x=params[0],
                    y=params[1],
                    x_name=params[2],
                    y_name=params[3],
                    outdir = self._output, 
                    identifier = None,
                    pseudo = self._pseudo)

                curr.compare()
'''