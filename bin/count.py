import sys
import os
import pandas as pd
import itertools
import argparse
from matching import selector_match, read_features
from utilities import timing

"""
    usage: count.py [-h] -i INPUT -a ANNOTATION -f FEATURES [-n NORMALIZETO] -o OUTPUT

    optional arguments:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
    -a ANNOTATION, --annotation ANNOTATION
    -f FEATURES, --features FEATURES
    -n NORMALIZETO, --normalizeTo NORMALIZETO
    -o OUTPUT, --output OUTPUT
"""

#### Functions
def intersect_to_best(input, output) : 
    
    """
    Iterate through intersect for each sequence find best overlap, only keep line if overlap == best.
    """
    
    best = {}
    with open(input, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            seq = info[7]
            overlap = int(info[19])
            
            if seq in best.keys() : 
                if overlap < best[seq] : 
                    best[seq] = overlap
            else : 
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

##################################################
# Class of features

class Features() :

    def __init__(self, alignment, annotation, features, norm_features, outprefix, rpkm) : 

        self._alignment = alignment
        self._annotation = annotation
        self._features = features
        self._norm_feats = norm_features
        self._normalize = True if norm_features is not None else False
        self._outprefix = outprefix
        self._normalize_rpkm = True if rpkm else False

        if os.path.exists(f"{self._outprefix}.log") : 
            os.remove(f"{self._outprefix}.log")

        self._log = open(f"{self._outprefix}.log", 'a')
    
    def intersect(self) :

        """
        Intersect reads to genomic features

        3/31/2022
            For any given sequence, group by sequence and pick alignments with the least number of mismatches.
        """
        self._intersect = f"{self._outprefix}.intersect"
        cmd = f"bedtools intersect -wo -a {self._annotation} -b {self._alignment} > {self._intersect}"
        if not os.path.exists(self._intersect) : 
            os.system(cmd)
        else : 
            pass
        
    def total_reads(self) :

        """
        Count the total number of reads. If a read maps to multiple features per alignment
        only the first instance will be counted.
        """ 

        seqs = {}
        self._total = 0
        self._normalization_factors = {}
        with open(self._alignment, 'r') as f :
            for line in f :  
                info = line.strip().split("\t")
                self._total += float(info[4])

        self._log.write(f"total number of reads = {round(self._total)}\n")
        self._normalization_factors['total'] = self._total
    
    def features(self) : 

        """
        Parse feature selectors and normalization selectors into quickly searchable 
        tuples.
        """

        self._feature_selectors = read_features(self._features)
        if self._normalize : 
            self._norm_selectors = read_features(self._norm_feats)  

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
                #seq = 7
                #locus = 8
                biotype = 7 #9
                Class = 8 #10
                sequence = 12 #14
                strand = 14 #16 
                multimap = 15 #17
                mismatch = 16 #18
                match = 17 #19

                orientation = "sense" if info[ref_strand] == info[strand] else "anti"

                if info[strand] == "+" and orientation == "sense" : 
                    distance = int(info[10]) - int(info[1])
                elif info[strand] == "-" and orientation == "sense" : 
                    distance = int(info[2]) - int(info[11])
                else : 
                    distance = 0

                classes = tuple(info[Class].split(","))
                
                selector = ( 
                    info[source], 
                    info[feature], 
                    info[biotype], 
                    classes,
                    info[gene].split(","),
                    #(info[gene], info[seq], info[locus]),
                    orientation, 
                    str(info[sequence][0]), 
                    len(info[sequence]),
                    int(info[match]),
                    distance, 
                    int(info[multimap]),
                    int(info[mismatch])
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

                        # in the form [ (chrom, start, end, seq, count [ntm], strand, strand orientation, (gene, seq, locus), biotype, )
                        hits.append(
                                [ 
                                info[9], 
                                info[10], 
                                info[11], 
                                info[12], 
                                float(info[13]), 
                                info[14], 
                                "sense" if info[ref_strand] == info[strand] else "anti", 
                                info[gene], 
                                info[biotype], 
                                info[8],
                                int(info[match])
                                ] + [ y[1] ] + [ y[2] ]
                        )
                        
                        if recorder in counter.keys() : 
                            counter[recorder] += float(info[13])
                        else :
                            counter[recorder] = float(info[13])
                
                if self._normalize : 
                    for x,y in itertools.product( [selector], self._norm_selectors ) : 
                        if selector_match(x, y[0]) :
                            name = y[1]
                            if name in self._normalization_factors.keys() : 
                                self._normalization_factors[name] += float(info[13])
                            else : 
                                self._normalization_factors[name] = float(info[13])
    
        f.close()

        # make counter an instance of the class
        self._counter = counter

        # turn hits into pandas dataframe in BED format
        tmp = pd.DataFrame(hits, columns = ['chrom', 'start', 'end', 'seq', 'count', 'strand', 'orientation', 'gene', 'biotype', 'class','overlap', 'feature', 'rank'])
        
        # Filter sequences by rank
        idx = tmp.groupby( ['seq'] )['rank'].transform(min) == tmp['rank']
        bed = tmp[idx]
        bed = bed.reset_index(drop = True)

        bed['number_feature_mapped'] = bed.groupby( ['chrom', 'start', 'end', 'seq'] )['seq'].transform('size')
        bed['count_nfm'] = bed.apply(lambda row : row['count'] / row['number_feature_mapped'], axis = 1)

        # select desired columns
        bed_final = bed[['chrom', 'start', 'end', 'seq', 'count_nfm', 'strand', 'orientation', 'gene', 'biotype', 'class', 'feature', 'rank', 'overlap']]

        # change name of count_nfm to count
        bed_final.columns = ['chrom', 'start', 'end', 'seq', 'count', 'strand', 'orientation', 'gene', 'biotype', 'class', 'feature', 'rank',  'overlap']

        # split gene into gene_name, seq_id, locus_id
        #bed_final[['gene_name', 'seq_id', 'locus_id']] = pd.DataFrame(bed_final['gene'].tolist(), index=bed_final.index)

        bed_final = bed_final[['chrom', 'start', 'end', 'seq', 'count', 'strand', 'orientation', 'biotype', 'class', 'feature', 'gene', 'rank', 'overlap']]
        self._bed_final = bed_final

        # output normalization constants as tab delimited dataframe
        ndf = pd.DataFrame.from_dict(self._normalization_factors, orient='index')
        ndf = ndf.reset_index()
        ndf.to_csv(f"{self._outprefix}.normalization.constants.tsv", sep = "\t", index = False, header = False)

    def report_statistics(self) : 

        # Report number of reads passing filters.
        self._log.write("\nAssignment statistsics...(Note: This is pre-rank filtering):\n") 
        total = 0
        self._log.write("")
        for k,v in self._counter.items() : 
            self._log.write(f"{k}\t{ round(100*(v/self._total),4) }% \t {round(v,4)}\n")
            total += 100*(v/self._total)
            
        self._log.write(f"Unannotated/unassigned/failing filters = {round(100 - total, 4)}%\n")
        self._log.write("-------------------------------------------------------------------------\n")

        self._log.write("\nNormalization constant based statistics (No rank filtering):\n") 
        for k,v in self._normalization_factors.items() : 
            self._log.write(f"{k}\t {round(v,4)}\n")
        self._log.write("-------------------------------------------------------------------------\n")

    def combine_counts(self) : 

        """
        Combine counts by gene, biotype, class, feature
        """
        if self._normalize_rpkm : 
            self._bed_final['count_kmer'] = self._bed_final.apply(lambda row: row['count']/len(row['seq']), axis = 1)
            self._counts = self._bed_final.groupby( ['gene', 'biotype', 'class', 'feature'] )['count_kmer'].sum().reset_index()
        else :
            self._counts = self._bed_final.groupby( ['gene', 'biotype', 'class', 'feature'] )['count'].sum().reset_index()

    def normalize_counts(self) :

        """
        Normalize counts per gene as set forth by the normalization constants features
        """
        for k,v in self._normalization_factors.items() : 
            if not v == 0 :
                if self._normalize_rpkm : 
                    self._counts[f'count_{k}_kmer_norm'] = self._counts.apply(lambda row : 1e6*(row['count_kmer']/v), axis = 1)
                    self._bed_final[f'count_{k}_kmer_norm'] = self._bed_final.apply(lambda row : 1e6*(row['count_kmer']/v), axis = 1)
                else : 
                    self._counts[f'count_{k}_norm'] = self._counts.apply(lambda row : 1e6*(row['count']/v), axis = 1)
                    self._bed_final[f'count_{k}_norm'] = self._bed_final.apply(lambda row : 1e6*(row['count']/v), axis = 1)

        self._counts.to_csv(f"{self._outprefix}.counts.tsv", sep = "\t", index = False, header = True)
        self._bed_final.to_csv(f"{self._outprefix}.bed.tsv", sep = "\t", index = False, header = True)


##################################################
# parse args

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-i", "--input", type = str, required = True, help="BED format file (chrom, start, end, seq, count, strand). Tab delimmited. In pipeline count it normalized to number of times mapped.")
    parser.add_argument("-a", "--annotation", type = str, required = True)
    parser.add_argument("-f", "--features", type = str, required = True)
    parser.add_argument("-n", "--normalizeTo", type = str, required=False)
    parser.add_argument("-o", "--output", type = str, required=True)
    parser.add_argument("-rpkm", action = 'store_true', required=False)
    return parser.parse_args()

##################################################
# main function
@timing
def main() : 

    args = get_args()
    
    if args.normalizeTo == "" or args.normalizeTo == '' : 
        args.normalizeTo = None

    counter = Features(args.input, args.annotation, args.features, args.normalizeTo, args.output, args.rpkm)
    counter.intersect()
    counter.features()
    counter.total_reads()
    counter.select_matching()
    counter.combine_counts()
    counter.normalize_counts()
    #counter.report_statistics()

##################################################
if __name__ == "__main__" :

    main()
