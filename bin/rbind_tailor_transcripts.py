import os
import sys
import argparse
import pandas as pd
from time import time
import glob

class Rbind_tables() : 

    def __init__(self, files, outname, feature) : 
        self._files = files
        self._outname = outname
        self._feature = feature
    
    def rbind(self) : 

        total_files = len(self._files)
        for i,f in enumerate(self._files) : 
            df = pd.read_csv(f, sep = "\t")
            if self._feature : 
                df = df.query('feature == @self._feature').reset_index(drop = True)
                
            name = os.path.basename(f).split(".")[0]
            df['sample'] = name
            if i == 0 :
                df_out = df
            else : 
                df_out = pd.concat([df_out, df], ignore_index = True)
            
            print(f"processed {i+1} out of {total_files}")
                
        df_dedup = df_out.drop_duplicates().reset_index(drop = True)
        df_dedup.to_csv(f"{self._outname}", sep = "\t", index = False, header = True)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--files", type = str, required = True)
    required.add_argument("-o", "--outname", type = str, required = True)
    required.add_argument("-d", action = 'store_true', required = False)
    required.add_argument("-feature", required = False, type = str)

    return parser.parse_args()

def main() : 

    args = get_args()
    if args.d : 
        files = glob.glob(f"{args.files}/*.bed.tsv")
    else : 
        files = args.files.replace("[","").replace("]","").split(", ")

    rbind_obj = Rbind_tables(files, args.outname, args.feature)
    rbind_obj.rbind()

if __name__ == "__main__" : 

    main()



