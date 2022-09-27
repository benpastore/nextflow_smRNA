import os
import sys
import pandas as pd
import argparse

class process() :

    def __init__(self, alignment, normalization, outprefix, biotype, rule_name) :
        self._alignment = alignment
        self._normalization = normalization
        self._outprefix = outprefix
        self._biotype = biotype
        self._rule_name = rule_name
        
    def count(self) : 

        genes = {}
        with open(self._alignment, 'r') as f : 
            for line in f : 
                info = line.strip().split("\t") 

                gene = info[0]
                count = float(info[4])

                if gene in genes.keys() :     
                    genes[gene] += count 
                else : 
                    genes[gene] = count
        f.close()

        self._counts = pd.DataFrame.from_dict(genes, orient = 'index').reset_index()
        self._counts.columns = ['gene_name', 'count']

    def parse_normalizaiton(self) : 
        
        norm_consts = {}
        with open(self._normalization, 'r') as f : 
            for line in f : 
                info = line.strip().split("\t")
                norm_method = info[0]
                constant = float(info[1])
                norm_consts[norm_method] = constant
        f.close() 

        self._norm_consts = norm_consts

    def normalize(self) : 

        for k,v in self._norm_consts.items() : 
            self._counts[f"count-{k}-norm"] = self._counts.apply(lambda row : 1000000*(row['count']/v), axis = 1)

    def add_attributes(self) : 

        self._counts['seq_id'] = self._counts['gene_name']
        self._counts['locus_id'] = self._counts['gene_name']
        self._counts['class'] = 'unknown'
        self._counts['biotype'] = f"{self._biotype}"
        self._counts['feature'] = f"{self._rule_name}"
    
    def write(self) : 

        self._counts.to_csv(f"{self._outprefix}.counts.tsv", sep = "\t", header = True, index = False)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", type = str, required = True)
    required.add_argument("-n", "--normalization", type = str, required = True)
    required.add_argument("-o", "--output", type = str, required = True)
    required.add_argument("-b", "--biotype", type = str, required = True)
    required.add_argument("-r", "--rule_name", type = str, required=False, default=os.getcwd())
    
    return parser.parse_args()

def main() :

    args = get_args()

    run = process(args.input, args.normalization, args.output, args.biotype, args.rule_name)
    run.count() 
    run.parse_normalizaiton()
    run.normalize()
    run.add_attributes()
    run.write()

if __name__ == "__main__" : 

    main()
        




                

