#!/usr/bin/env python3

import os 
import sys
import argparse
import pandas as pd
import numpy as np
from time import time 
from scipy import stats

from compare import Compare


def parse(line) : 
        
    info = line.strip().split("\t")
    x_samples = info[0].replace("\"", "").split(",")
    y_samples = info[1].replace("\"", "").split(",")
    x_name = str(info[2])
    y_name = str(info[3])

    return x_samples, y_samples, x_name, y_name

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

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-c", "--comparisons", type = str, required = True, help = f"Comparisons table")
    required.add_argument("-f", "--file", type = str, required=True, help = "Counts table")
    required.add_argument("-o", "--output", type = str, required=False)
    required.add_argument("-pseudo", type = float, required = False)

    return parser.parse_args()
 
def main() :

    args = get_args()

    x = Compare_multi(args.comparisons, args.file, args.output, args.pseudo)
    x.compare_multi()

if __name__ == "__main__" : 

    main()