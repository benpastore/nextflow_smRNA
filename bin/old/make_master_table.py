
import pandas as pd 
import os 
import sys
import glob 

def read_tables(list_of_files, sample) :

    """
    A function to take a list of files and read them in,
    as they are read in combine the dataframes based on colname.

    Input : a list of files

    Output : a dataframe of "rbinded" files

    """

    for i,f in enumerate(list_of_files) : 

            d = pd.read_csv( f, sep = "\t", names = ['stable_ID', '{}'.format(sample)] )

            if i == 0 : 
                df_out = d
            else : 
                df_out = df_out.append(d, ignore_index = False)
    
    return df_out



def make_tables(directory, outfile) :

    """
    A function to take a directory containing gene counts for various subgroups of 
    RNAs and genrate a comprehensive counts table.

    1) First, need to rbind all tables that have the same prefix
    2) Next, outer join all rbinded dables

    Input : directory containing the genes/counts lists. Must be in the form gene (tab) count 

    Output : single master table. RowNames will be gene, column names will be sample 
    
    """

    #1. rbind tables with same prefix 
    os.chdir(directory)
    files = glob.glob("*.counts")
    prefixes = [ i.split(".")[0] for i in files ]
    prefixes_uniq = list(set(prefixes))

    for i,p in enumerate(prefixes_uniq) :

        files = glob.glob("{}*".format(p))

        df = read_tables(files, p)

        if i == 0 :
            df_out = df
        else : 
            df_out = df_out.merge(df, on = ['stable_ID'], how = 'outer')
    
    df_out = df_out[~df_out['stable_ID'].str.contains("\*") ]
    return df_out


def main() : 

    directory = sys.argv[1]
    outfile = sys.argv[2]

    res = make_tables(directory, outfile)
    res.to_csv(outfile, sep = '\t', index = False, header = True)

if __name__ == "__main__" : 

    main()