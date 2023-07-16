#!/usr/bin/env python3

import os 
import sys
import argparse
import time
import logging
import pandas as pd
import numpy as np
from scipy import stats
import multiprocessing
from multiprocessing import  Pool

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
    
def fold_change(x, y, pseudocount) : 

    if x > 0 and y > 0 : 
        return np.log2( y/x )
    else : 
        return np.log2( (y+pseudocount)/(x+pseudocount) )
        
def correct_pvalue(rank, number_significant) :

    return min(rank/number_significant, 1)

def parallelize_dataframe(df, func, n_cores=20):
    
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()

    return df

def ttest(df) :
    
    df['pvalue'] = df.apply(lambda row: stats.ttest_ind(row['IP'], row['input'], nan_policy = 'propagate', equal_var = True)[1], axis = 1)
    
    return df

class Compare() : 

    def __init__(self, counts, x, y, x_name, y_name, outdir, identifier, pseudo) :

        self._counts = pd.read_table(counts, sep = "\t")
        
        self._x = x
        self._y = y

        self._x_name = x_name
        self._y_name = y_name
        
        self._identifier = identifier

        if identifier is None : 
            self._ID = self._counts.loc[:, (self._counts.dtypes != int) & (self._counts.dtypes != float) ]
        else : 
            self._ID = self._identifier

        self._outprefix = os.path.basename(counts).replace(".tsv", "")
        
        if outdir is None : 
            self._outdir = os.getcwd()
        else : 
            if not os.path.exists(outdir) : 
                os.mkdir(outdir)
            self._outdir = outdir

        if "count" in os.path.basename(counts) : 
            norm_method = os.path.basename(counts).split("count")[1].replace(".tsv", "").replace("-", "_")
            self._norm_method = f"{norm_method}"
        else : 
            self._norm_method = ""

        self._pseudo = pseudo

    def compare(self) : 

        # subset the dataframe to have the gene, count
        # comparisons should be semi-colon separated list of genotypes in x and y
        # get x and y dataframes
        print("\tCalculating log2FoldChange & pvalue...")
        df = self._counts[[c for c in self._counts.columns if c in self._x or c in self._y or c in self._ID]]
        df = df.reset_index(drop=True)

        df[f"{self._x_name}"] = df[self._x].mean(axis = 1)
        df[f"{self._y_name}"] = df[self._y].mean(axis = 1)

        df["basemean"] = df[[f"{self._x_name}", f"{self._y_name}"]].mean(axis = 1)
        df = df[ df['basemean'] > 0 ].reset_index(drop=True)
        
        median_basemean = df['basemean'].median()
        
        pseudocount = 0.1 if self._pseudo is None else self._pseudo
        
        pseudocount = round(pseudocount, 2)
        
        print(f"\tPseudocount: {pseudocount}")

        # calculate log2FoldChange & pvalue
        df['log2FoldChange'] = df.apply(lambda row : fold_change(row[f'{self._x_name}'], row[f'{self._y_name}'], pseudocount), axis = 1 )
        
        fin = parallelize_dataframe(df, ttest, n_cores = 30)
        #try : 
        #    df['pvalue'] = df.apply(lambda row : stats.ttest_ind( row[self._x], row[self._y], alternative = "two-sided", nan_policy = 'propagate', equal_var=True)[1], axis = 1)
        #except : 
        #    df['pvalue'] = df.apply(lambda row : stats.ttest_ind( row[self._x], row[self._y], nan_policy = 'propagate', equal_var=True)[1], axis = 1)
            
        fin = fin.sort_values(by = 'gene').reset_index(drop=True)

        # export
        print(f"\tWriting files... {self._outdir}/{self._x_name}_vs_{self._y_name}{self._norm_method}.tsv")
        
        fin.to_csv(f"{self._outdir}/{self._x_name}_vs_{self._y_name}.{self._norm_method}.tsv", sep = "\t", header = True, index = False)
        
        print("\tAll done!")
        
def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-c", "--counts", type = str, required = True)
    required.add_argument("-x", "--xlist", type=lambda s: [item for item in s.replace("[","").replace("]","").replace(", ",",").split(",") ], action='store',   required = True)
    required.add_argument("-y", "--ylist", type=lambda s: [item for item in s.replace("[","").replace("]","").replace(", ",",").split(",") ], action='store',   required = True)
    required.add_argument("-nx", "--name_x", type = str, required = True)
    required.add_argument("-ny", "--name_y", type = str, required = True)
    required.add_argument("-o", "--outdir", type = str, required = False)
    required.add_argument("-i", "--identifier", type=lambda s: [item for item in s.split(',')], action='store',   required = False)
    required.add_argument("-pseudo", type = float, required = False)

    return parser.parse_args()

@timing
def main() : 
    
    logging.basicConfig(level=logging.DEBUG)

    args = get_args()   
    comp = Compare(args.counts, args.xlist, args.ylist, args.name_x, args.name_y, args.outdir, args.identifier, args.pseudo)
    comp.compare()
    
if __name__ == "__main__" : 

    main()