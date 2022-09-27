#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

"""
Given a dataframe of UNnormalized counts: 

1. calculate the reference for each gene (n1*n2...ni)^1/n
2. Calculate the ratio SampleN / reference
3. Find the median ratio per sample
4. Divide the median ratio per sample.

"""

def pseudo_reference(x) : 
    
    l = [ i for i in x if i > 0 ]
    
    n_samples = len(l)
    
    prod = np.prod(l)
    
    pseudo = np.power(prod, 1/n_samples)
    
    return(pseudo)

    np.power(np.prod(list(row)), 1/n_samples)

def geometric_normalization(file, outname) : 
    
    dataframe = pd.read_csv(file, sep = "\t")
    
    info = dataframe.select_dtypes(include = ['object'])
    
    df = dataframe.select_dtypes(exclude=['object']).reset_index(drop = True)
    
    n_samp = df.shape[1]
    
    df['pseudo_reference'] = df.apply(lambda row : pseudo_reference(list(row)), axis = 1)
    
    df_ratios = df.iloc[:,0:n_samp+1].divide(df.iloc[:,-1], axis = 'rows')
    
    df_medians = df_ratios.median().to_frame().T
    
    df_norm = df.div(df_medians.iloc[0], axis = 'columns')
    
    res = info.merge(df_norm, left_index = True, right_index = True)
    
    df_medians.to_csv(f"{outname}_ratio_of_medians.tsv", sep = "\t", header = True, index = False)
    res.to_csv(f"{outname}.count_geometric_norm.tsv", sep = "\t", header = True, index = False)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", type = str, required = True)
    required.add_argument("-o", "--outname", type = str, required = True)

    return parser.parse_args()

def main() : 
    
    args = get_args() 
    geometric_normalization(args.input, args.outname)

if __name__ == "__main__" : 
    
    main()