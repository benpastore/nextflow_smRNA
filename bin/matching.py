#!/usr/bin/env python3

import pandas as pd

def selector_match(x,y) :

    """
    General function. Searches if x in y for tuples. 
    y is a tuple of typles/lists. x is an int/str.

    All x's must be in the respective y tuple/list to pass filters.
    """
    
    return all( True if a in b else True if "*" in b else True if (type(a) == tuple and any(q in b for q in a)) else False for a,b in zip(x,y) )

def read_features(features) : 

    """
    Read features into pandas dataframe
    """

    # 0. read in the features dataframe with pandas
    df = pd.read_csv(features, sep = "\t",
        names = [ 'rule_name','source','feature','biotype','gene','Class','firstNT','length','match','strand','distance', 'multimap', 'mismatch', 'rank'],
        dtype=str, header = 0)

    # 1. split feature params into tuples for searching
    for index,row in df.iteritems() :
        if not index == "rule_name" and not index == "rank" : 
            # if length
            if index == "length" :
                for i,k in enumerate(df[f'{index}']) :
                    if "-" in k : 
                        Min = int(df[f'{index}'][i].split("-")[0])
                        Max = int(df[f'{index}'][i].split("-")[1]) + 1
                        df[f'{index}'][i] = tuple(range(Min,Max))
                    elif "*" in k : 
                        df[f'{index}'][i] = tuple("*")
                    else : 
                        df[f'{index}'][i] = (int(k),)

            # if match
            elif index == "match" : 
                for i,k in enumerate(df[f'{index}']) :
                    if k in ['*', 'partial'] : 
                        df[f'{index}'][i] = tuple("*")
                    else : 
                        df[f'{index}'][i] = tuple(range(int(k), 100))

            # if distance
            elif index == "distance": 
                for i,k in enumerate(df[f'{index}']) : 
                    if "*" in k : 
                        df[f'{index}'][i] = tuple("*")
                    elif "_" in k : 
                        x1 = int(k.split("_")[0])
                        x2 = int(k.split("_")[1])
                        df[f'{index}'][i] = tuple(range(x1,x2+1))
                    elif ";" in k : 
                        df[f'{index}'][i] = tuple(i for i in k.split(";"))
                    else :
                        df[f'{index}'][i] = (int(k),)
            
            elif index == "mismatch" or index == "multimap" : 
                for i,k in enumerate(df[f'{index}']) :
                    if k in ['*'] : 
                        df[f'{index}'][i] = tuple("*")
                    else : 
                        df[f'{index}'][i] = tuple(range(int(k)+1))
            
            #lif index == "Class" : 
            #    for i,k in enumerate(df[f'{index}']) :
            #        if k in ['*'] : 
            #            df[f'{index}'][i] = tuple("*")
            #        else : 
            #            df[f'{index}'][i] = tuple(df[f"{index}"].str.split(";"))

            # if anything else
            else : 
                df[f'{index}'] = tuple(df[f"{index}"].str.split(";"))
        
    feature_selectors = []
    feature_dict = {}

    for index,row in df.iterrows() :
        name = row['rule_name']
        rank = row['rank']
        feature_dict[name] = row

        selector = ( 
            (
            feature_dict[name]['source'], 
            feature_dict[name]['feature'], 
            feature_dict[name]['biotype'], 
            feature_dict[name]['Class'],
            feature_dict[name]['gene'],
            feature_dict[name]['strand'], 
            feature_dict[name]['firstNT'],
            feature_dict[name]['length'],
            feature_dict[name]['match'],
            feature_dict[name]['distance'],
            feature_dict[name]['multimap'],
            feature_dict[name]['mismatch']
            ), name, rank
        )

        feature_selectors.append(selector)

    # Note -> selector has the form [ (((source), (feature), (biotype), (strand), (firstNT), (length), (match), (distance), (multimap) ), rule_name) ]
    return feature_selectors

def read_normalization(norm) : 

    norm_const = {}
    with open(norm, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            norm_const[info[0]] = float(info[1])
    f.close()

    return norm_const
