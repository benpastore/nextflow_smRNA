#!/usr/bin/env python3

import argparse
import glob
import pandas as pd
import numpy as np 
import os

class Filter() : 

    def __init__(self, directory, condition, output, feature) : 

        self._directory = directory
        self._condition = condition
        self._output = output
        self._feature = feature
    
    def execute(self) : 

        self.filter()
    
    def filter(self) : 

        """
        Given N number of files with tailed alignments, filter alignments that
        appear a minimum of two times. For instance the sequence: 
            AATTGGAGAGACCT with the tail "TT" must appear in two separate replicates.        
        """

        conditions = {}
        with open(self._condition, 'r') as f : 
            for line in f : 
                info = line.strip().split("\t")
                sample = info[0]
                condition = info[1]
                conditions[sample] = condition
        f.close()
        
        files = glob.glob(f"{self._directory}/*.bed.tsv")
        for i,f in enumerate(files) : 
            k = os.path.basename(f).split(".")[0]
            data = pd.read_csv(f, sep = "\t")
            
            v = conditions[k]
            
            data['sample'] = k
            data['condition'] = v
            
            if self._feature is None : 
                curr = data
            else : 
                curr = data[ data['feature'] == self._feature ].reset_index()

            if i == 0 : 
                master = curr
            else : 
                master = master.append(curr)
            
            i += 1

        # filter by replicate
        master['size'] = master.groupby(['gene_name', 'seq_id', 'locus_id', 'biotype', 'feature', 'class', 'seq', 'tail', 'condition'])['seq'].transform('size')
        passing = master[ master['size'] >= 2]
        passing = passing.reset_index(drop = True)

        # filter by conditional 
        conditional = master[ master['condition'] == master['sample'] ]
        conditional = conditional.reset_index(drop = True)
        conditional['condition'] = conditional['sample']
        
        # combine passing and conditional
        final = passing.append(conditional)
        final = final.reset_index(drop = True)
        
        final = final.drop(columns=['index'])
        final.to_csv(self._output, sep = "\t", header = True, index = False)
        
def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument("-c", "--condition", type = str, required = True)
    parser.add_argument("-d", "--directory", type = str, required = True)
    parser.add_argument("-o", "--output", type = str, required = True)
    parser.add_argument("-f", "--feature", type = str, required = False)

    return parser.parse_args()

##################################################
# main function
def main() : 

    args = get_args()

    directory = args.directory
    condition = args.condition
    output = args.output
    feature = args.feature

    run = Filter(directory, condition, output, feature)
    run.execute()

if __name__ == '__main__' : 

    main()


