
import os
import sys
import pandas as pd

class create_table() : 

    def __init__(self, files, outprefix) :
        self._files = files
        self._outprefix = outprefix
    
    def cbind_tables(self) :

        sample = pd.read_table(self._files[0], sep = "\t")
        print(sample.head())

        for i,f in enumerate(self._files) : 
            sample_name = os.path.basename(f).split("-normalization-")[0]
            df = pd.read_table(f, sep = "\t", names = ['method', f"{sample_name}"])

            if i == 0 : 
                df_out = df
            else : 
                df_out = df_out.merge(df, on = ['method'], how = 'outer')
                df_out = df_out.fillna(0)

            df_out.to_csv(f"{self._outprefix}-normalization-methods.tsv", header = True, index = False, sep = "\t")


def main() : 

    files = sys.argv[1].replace("[","").replace("]","").split(", ")
    outprefix = sys.argv[2]
    
    x = create_table(files, outprefix)
    x.cbind_tables()

if __name__ == "__main__" : 

    main()

        

