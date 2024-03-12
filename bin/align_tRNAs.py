#!/usr/bin/env python3

import os 
import re 
import glob
import argparse
import pandas as pd

def align(idx, fasta, basename, normalization) : 

    mismatches = [0,1,2,3]
    multimap = 1000
    input_fa = fasta

    for mismatch in mismatches : 
        cmd = f"bowtie -x {idx} -f {input_fa} -v {mismatch} -m {multimap} -p 8 --un {basename}.unmapped.v{mismatch}.m{multimap}.fa -a --best --strata -S {basename}.aligned.v{mismatch}.m{multimap}.sam"
        os.system(cmd)
        input_fa = f"{basename}.unmapped.v{mismatch}.m{multimap}.fa"
    
    # get samfiles and convert to bed files, for top 1% of reads aligning with 2/3 mismatches, remove reads that align to multiple types of tRNA
    lines = ''
    sam = glob.glob("./*.sam")
    pattern = r'\.v(\d+)\.m'
    chroms = []
    starts = []
    ends = []
    seqs = []
    counts = []
    strands = []
    mismatches = []
    for s in sam : 
        match = re.search(pattern, s)
        n_mismatch = match.group(1)
        cmd = f"sam2bed < {s} | cut -f 1-6 > tmp.bed"
        os.system(cmd)
        with open('tmp.bed', 'r') as f : 
            for line in f : 
                info = line.strip().split("\t")
                chrom = info[0]
                start = int(info[1])
                end = int(info[2])
                seq = str(info[3]).split(":")[0]
                count = float(info[3].split(":")[1])
                strand = str(info[5])
                mismatch = int(n_mismatch)

                chroms.append(chrom)
                starts.append(start)
                ends.append(end)
                seqs.append(seq)
                counts.append(count)
                strands.append(strand)
                mismatches.append(mismatch)
        f.close()
    
    df = pd.DataFrame({
        'gene' : chroms, 
        'start' : starts, 
        'end' : ends, 
        'seq' : seqs, 
        'count' : counts, 
        'strand' : strands,
        'mismatch' : mismatches
    })

    # normalize to the number of mapped reads per sequence
    df_agg = df.groupby(['seq']).agg(ntm = ('seq', 'size')).reset_index()

    df = df.merge(df_agg, on = ['seq'], how = 'left')

    df['ntm'] = df.apply(lambda x : x['count'] / x['ntm'], axis = 1)

    res = df[['gene', 'start', 'end', 'seq', 'ntm', 'strand']].reset_index(drop = True)
    
    # normalize to normalization constants
    normalization_features = {}
    with open(normalization, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t") 
            normalization_features[info[0]] = float(info[1])

    res['biotype'] = 'tRNA'
    res['class'] = 'none'
    res['feature'] = 'tRNA_alignment_pipeline'

    counts = res.groupby(['gene', 'biotype', 'class', 'feature']).agg(count = ('ntm', 'sum')).reset_index()

    for k,v in normalization_features.items() :
        res[f'count_{k}_norm'] = res.apply(lambda row : 1e6*(row['ntm']/v), axis = 1)
        counts[f'count_{k}_norm'] = counts.apply(lambda row : 1e6*(row['count']/v), axis = 1)

    res.to_csv(f'{basename}.tRNA.aligned.bed.tsv', sep = "\t", header = True, index = False)
    counts.to_csv(f'{basename}.tRNA.aligned.counts.tsv', sep = "\t", header = True, index = False)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument("--index", type = str, required = True)
    parser.add_argument("--fasta", type = str, required = True)
    parser.add_argument("--basename", type = str, required = True)
    parser.add_argument("--norm_consts", type = str, required = True)

    return parser.parse_args()

def main() :

    args = get_args()

    align(args.index, args.fasta, args.basename, args.norm_consts)

if __name__ == "__main__" : 

    main()