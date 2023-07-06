
import os
import sys
import glob
import pandas as pd
import argparse

class create_table() : 

    def __init__(self, files, outprefix) :
        self._files = files
        self._outprefix = outprefix
    
    def cbind_tables(self) :
        
        sample = pd.read_table(self._files[0], sep = "\t")
        cols = []
        non_count_cols = []
        
        for col in sample.columns :
            if "count" in col :
                cols.append(col)
            else : 
                non_count_cols.append(col)
                
        print(cols)
                
        for col in cols : 
            names = []
            for i,f in enumerate(self._files) : 
                sample_name = os.path.basename(f).split(".")[0]
                print(sample_name)
                names.append(sample_name)
                df = pd.read_table(f, sep = "\t")
                print(df.head())
                save = non_count_cols + [ col ] 
                sub = df[save]
                my_names = non_count_cols + [ f'{sample_name}' ] 
                sub.columns = my_names
            
                if i == 0 : 
                    df_out = sub
                else : 
                    df_out = df_out.merge(sub, on = non_count_cols, how = 'outer')
                    df_out = df_out.fillna(0)
            
               
            #names = names.sort()
            col_order = non_count_cols + names
            df_out = df_out[col_order]
            df_dedup = df_out.drop_duplicates().reset_index(drop = True)
            df_dedup.to_csv(f"{self._outprefix}.{col}.tsv", header = True, index = False, sep = "\t")

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)
    
    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--files", type = str, required = True)
    required.add_argument("-o", "--outname", type = str, required = True)
    required.add_argument("-d", action = 'store_true', required = False)

    return parser.parse_args()

def main() : 
    
    args = get_args() 
    
    if args.d : 
        files = glob.glob(f"{args.files}/*.tsv")
        print(files)
    else : 
        if "[" in args.files : 
            files = args.files.replace("[","").replace("]","").split(", ")
        else : 
            files = args.files.replace("[","").replace("]","").split(",")

    outprefix = args.outname
    
    x = create_table(files, outprefix)
    x.cbind_tables()

if __name__ == "__main__" : 

    main()

        

