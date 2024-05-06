import os
import sys
import argparse
import pandas as pd
from time import time


class Rbind_tables() : 

    def __init__(self, files, outname) : 
        self._files = files
        self._outname = outname
    
    def rbind(self) : 

        for i,f in enumerate(self._files) : 
            df = pd.read_csv(f, sep = "\t")
            name = os.path.basename(f).split(".")[0]
            df['sample'] = name
            if i == 0 :
                df_out = df
            else : 
                df_out = df_out.append(df)
        
        df_dedup = df_out.drop_duplicates().reset_index(drop = True)
        df_dedup.to_csv(f"{self._outname}", sep = "\t", index = False, header = True)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--files", type = str, required = True)
    required.add_argument("-o", "--outname", type = str, required = True)

    return parser.parse_args()

def main() : 

    args = get_args()
    files = args.files.replace("[","").replace("]","").split(", ")

    rbind_obj = Rbind_tables(files, args.outname)
    rbind_obj.rbind()

if __name__ == "__main__" : 

    main()



