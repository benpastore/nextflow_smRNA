import sys
import os
import pandas as pd
import itertools
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from matching import selector_match, read_features, read_normalization
from tailor_match import filter_tails, parse_normalization_constants

"""
    usage: tailor_count.py [-h] -i INPUT -a ANNOTATION -f FEATURES [-n NORMALIZETO] -o OUTPUT

    optional arguments:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
    -a ANNOTATION, --annotation ANNOTATION
    -f FEATURES, --features FEATURES
    -n NORMALIZETO, --normalizeTo NORMALIZETO
    -o OUTPUT, --output OUTPUT

    This script is the same as count.py, however it is adapted to filter tailed reads.
"""

#### Functions
'''
def intersect_to_best(input, output) : 
    
    """
    Iterate through intersect for each sequence find best overlap, only keep line if overlap == best.
    """
    
    output_dict = {}
    best = {}
    with open(input, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")

            seq = info[7]
            overlap = int(info[19])
            
            if seq in best.keys() : 
                if overlap >= best[seq] :
                    best[seq] = overlap
                    output_dict[seq] += line
            else : 
                output_dict[seq] = line
                best[seq] = overlap
    f.close() 
    
    out = open(output, 'w')
    with open(input, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            seq = info[7]
            overlap = int(info[19])
            
            if overlap == best[seq] : 
                out.write(line)
            else : 
                pass
    f.close() 
    out.close()
'''

##################################################
# Class of features

class Features() :

    def __init__(self, alignment, annotation, features, norm_constants, outprefix, fasta) : 

        self._alignment = alignment
        self._annotation = annotation
        self._features = features
        self._norm_constants = norm_constants
        self._normalize = True if norm_constants is not None else False
        self._outprefix = outprefix
        self._records = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta')) 

        if os.path.exists(f"{self._outprefix}.log") : 
            os.remove(f"{self._outprefix}.log")

        self._log = open(f"{self._outprefix}.log", 'a')
    
    def intersect(self) :

        self._intersect = f"{self._outprefix}.intersect"
        cmd = f"bedtools intersect -wo -a {self._annotation} -b {self._alignment} > {self._intersect}"
        
        if not os.path.exists(self._intersect) : 
            os.system(cmd)
            #intersect_to_best("tmp", self._intersect)
        else : 
            pass
    
    def features(self) : 

        """
        Parse feature selectors and normalization selectors into quickly searchable 
        tuples.
        """
        self._feature_selectors = read_features(self._features)      

    def select_matching(self) : 

        """
        First layer of selection. For each alignment check to see if ANY rule is satasfied.
        If a rule is satisfied store the bed information (modified) and rule in a list.

        To search for matching rules for each alignment we use itertools to search the cartesian product 
        of the alignment selector (source, feature, biotype, strand orientation, firstNT) against a list 
        of feature selectors [] ( [source], [feature], [biotype], [strand orientation], ['firstNT'] ), ... ].
        If a hit is found, the for loop is broken and the [alignment, feature] information is appended to a 
        hits` list.

        Once the reads of interest have been selected, go on to do specific feature based selection. 
        This will evaluate the length, first NT, overlap, 5p to 5p distance. Reads passing this layer 
        of selection will be added to `winners` list. This list will be turned into a pandas dataframe 
        and then will be used for normalization.
        """

        hits = []
        counter = {}
        with open(self._intersect, 'r') as f : 
            for line in f : 

                info = line.strip().split("\t")

                source = 3
                feature = 4
                ref_strand = 5
                gene = 6
                seq = 7
                locus = 8
                biotype = 9
                Class = 10
                sequence = 14
                strand = 16
                tail = 17
                multimap = 18
                match = 19

                mismatch = 0

                orientation = "sense" if info[ref_strand] == info[strand] else "anti"

                if info[strand] == "+" and orientation == "sense" : 
                    distance = int(info[12]) - int(info[1])
                elif info[strand] == "-" and orientation == "sense" : 
                    distance = int(info[2]) - int(info[13])
                else : 
                    distance = 0

                selector = ( 
                    info[source], 
                    info[feature], 
                    info[biotype], 
                    info[Class],
                    (info[gene], info[seq], info[locus]),
                    orientation, 
                    str(info[sequence][0]), 
                    len(info[sequence])-len(info[tail]),
                    (int(info[match])),
                    distance,
                    int(info[multimap]),
                    mismatch
                    )
                
                recorder = (
                    info[source], 
                    info[feature], 
                    info[biotype], 
                    "sense" if info[ref_strand] == info[strand] else "anti", 
                    str(info[sequence][0])
                    )

                for x,y in itertools.product( [selector], self._feature_selectors ) :
                    if selector_match(x,y[0]) :
                        if filter_tails(info[sequence], info[tail], len(info[tail]), info[11], int(info[12]), int(info[13]), info[strand], self._records, 4) :
                            hits.append(
                                    [ 
                                    info[11], 
                                    int(info[12]),
                                    int(info[13]), #+len(info[tail]) 
                                    info[14], 
                                    float(info[15]), 
                                    info[16], 
                                    "sense" if info[ref_strand] == info[strand] else "anti", 
                                    (info[6], info[7], info[8]), 
                                    info[biotype], 
                                    info[10], 
                                    info[tail],
                                    int(info[match])
                                    ] + [ y[1] ] + [ y[2] ]
                            )
                            
                            if recorder in counter.keys() : 
                                counter[recorder] += float(info[15])
                            else :
                                counter[recorder] = float(info[15])
        f.close()

        # make counter an instance of the class
        self._counter = counter

        # turn hits into pandas dataframe in BED format
        tmp = pd.DataFrame(hits, columns = ['chrom', 'start', 'end', 'seq', 'count', 'strand', 'orientation', 'gene', 'biotype', 'class', 'tail' , 'overlap', 'feature', 'rank']) 

        idx = tmp.groupby( ['seq'] )['rank'].transform(min) == tmp['rank']
        bed = tmp[idx]
        bed = bed.reset_index(drop = True)

        # normalize to number features mapped
        bed['number_feature_mapped'] = bed.groupby( ['chrom', 'start', 'end', 'seq'] )['seq'].transform('size')
        bed['count_nfm'] = bed.apply(lambda row : row['count'] / row['number_feature_mapped'], axis = 1)
        bed['sequence'] = bed.apply(lambda row: (f'{row["seq"]}:{row["tail"]}'), axis = 1)
        
        # select desired columns
        bed_final = bed[['chrom', 'start', 'end', 'sequence', 'count_nfm', 'strand', 'orientation', 'gene', 'biotype', 'class', 'feature', 'rank', 'overlap']]

        # change name of count_nfm to count
        bed_final.columns = ['chrom', 'start', 'end', 'seq', 'count', 'strand', 'orientation', 'gene', 'biotype', 'class', 'feature', 'rank', 'overlap']

        # split gene into gene_name, seq_id, locus_id
        bed_final[['gene_name', 'seq_id', 'locus_id']] = pd.DataFrame(bed_final['gene'].tolist(), index=bed_final.index)
        
        # remove `gene` column
        self._bed_final = bed_final.drop(['gene'], axis = 1)

    def combine_counts(self) : 

        """
        Combine counts by gene, biotype, class, feature, tail
        """

        self._counts = self._bed_final.groupby( ['gene_name', 'seq_id', 'locus_id', 'biotype', 'class', 'tail', 'feature'] )['count'].sum().reset_index()

    def normalize_counts(self) :

        """
        Normalize counts per gene as set forth by the normalization constants features
        """ 
        
        if self._normalize : 
            self._normalization_factors = parse_normalization_constants(self._norm_constants)

            for k,v in self._normalization_factors.items() : 
                if not v == 0 : 
                    #self._counts[f'count_{k}_norm'] = self._counts.apply(lambda row : 1e6*(row['count']/v), axis = 1)
                    self._bed_final[f'count_{k}_norm'] = self._bed_final.apply(lambda row : 1e6*(row['count']/v), axis = 1)
    
        #self._counts.to_csv(f"{self._outprefix}.counts.tsv", sep = "\t", index = False, header = True)
        self._bed_final.to_csv(f"{self._outprefix}.bed.tsv", sep = "\t", index = False, header = True)

##################################################
# parse args

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument(
        "-i", 
        "--input", 
        type = str, 
        required = True,
        help="BED format file (chrom, start, end, seq, count, strand). Tab delimmited. In pipeline count it normalized to number of times mapped.")

    parser.add_argument(
        "-a", 
        "--annotation", 
        type = str, 
        required = True)

    parser.add_argument(
        "-f", 
        "--features", 
        type = str, 
        required = True)
    
    parser.add_argument(
        "-n",
        "--normalize",
        type = str,
        required=False
    )

    parser.add_argument(
        "-o",
        "--outprefix",
        type = str,
        required=True
    )

    parser.add_argument(
        "-g",
        "--fasta",
        type = str,
        required=True
    )

    return parser.parse_args()

##################################################
# main function
def main() : 

    args = get_args()

    counter = Features(args.input, args.annotation, args.features, args.normalize, args.outprefix, args.fasta)
    counter.intersect()
    counter.features()
    counter.select_matching()
    #counter.combine_counts()
    counter.normalize_counts()

##################################################
if __name__ == "__main__" : 

    main()