#!/usr/env/bin python3

import argparse
import os
import pandas as pd 
from bed_to_ntm import bed_to_ntm

def align(ref, fasta, mism, multim, tailor_software) : 

    tmp = os.path.basename(fasta).replace(".fa", ".tmp")
    cmd = f"{tailor_software} map -p {ref} -i {fasta} -n 25  > {tmp}.sam"
    os.system(cmd)

    lines = ''
    with open(f"{tmp}.sam") as f : 
        for line in f : 
            info = line.strip().split("\t")
            chrom = info[0]
            start = int(info[1])
            end = int(info[2])
            seq = info[3].split(":")[0]
            count = float(info[3].split(":")[1])
            strand = info[5]
            tail = info[6]
            nmap = float(info[7])
            ntm = count / nmap

            lines += f"{chrom}\t{start}\t{end}\t{seq}\t{ntm}\t{strand}\t{tail}\t{nmap}\n"
    f.close() 
    op = open(f"{tmp}.ntm", 'w')
    op.write(lines)
    op.close()

    return f"{tmp}.ntm"

#########################################
class Transcripts() :

    """
    A class to align reads (fasta format >seq:count) to various transcript references, 
    useful for implementation in pipelines to align smRNA reads to miRBase & repBase annotation
    """ 

    def __init__(self, fasta, transcripts, normalization, outname, mismatch, multimap, index_path, build_index, normalize_rpkm, tailor_path) : 

        self._fasta = fasta
        self._transcripts = transcripts 
        self._outname = outname
        self._normalization = normalization
        self._mismatch = mismatch
        self._multimap = multimap
        self._build_index = build_index
        self._normalize_rpkm = True if normalize_rpkm else False
        self._tailor_path = tailor_path

        if self._normalization is not None : 
            self._normalization_features = {}
            with open(self._normalization, 'r') as f : 
                for line in f : 
                    info = line.strip().split("\t") 
                    self._normalization_features[info[0]] = float(info[1])
        
        self._index_path = index_path if index_path is not None else os.path.join( os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "index/transcripts")
    
    def process(self) : 
        
        gene_ids,  seq_ids, locus_ids, biotypes, Classes, features, counts, starts, ends, strands, seqs , tails, rawtails = [], [], [], [], [], [], [], [], [], [], [], [], []
        with open(self._transcripts, 'r') as q : 
            for line in q :
                if not line.startswith("Name") : 
                    info = line.strip().split("\t")
                    
                    # feature name
                    feature = info[0]
                    biotype = info[1]
                    Class = info[2]

                    # five prime nucleotide
                    if info[3] == "*" : 
                        nt = [ "*" ]
                    else : 
                        nt = info[3].split(";") 

                    # length
                    if "-" in info[4] : 
                        min = int(info[4].split("-")[0])
                        max = int(info[4].split("-")[1]) + 1
                        length = [ i for i in range(min,max) ]
                    elif info[4] == "*" : 
                        length = ["*"]
                    else : 
                        length = [ int(info[4]) ]
                    
                    # strand oritentation
                    rule_orientation = info[5]
                    reference = info[8]

                    # alignment settings 
                    mismatch = [ i for i in range(0, int(info[6])+1 ) ]
                    multimap = [ i for i in range(0, int(info[7])+1 ) ]

                    # rule selector
                    rule_selector = (nt, length, rule_orientation, multimap, mismatch)
                    
                    # index genome
                    if self._build_index : 
                        idx = index(reference)
                    else : 
                        idx = os.path.join(self._index_path, os.path.basename(reference).replace(".fa",""))

                    # align

                    ntm = align(idx, self._fasta, self._mismatch, self._multimap, self._tailor_path)

                    # process alignment
                    with open(ntm, 'r') as f : 
                        for line in f : 
                            info = line.strip().split("\t")
                            gene = info[0]
                            start = int(info[1])
                            end = int(info[2])
                            seq = info[3]
                            count = float(info[4])
                            strand = info[5]
                            tailtmp = int(info[6])
                            tail = list(set(tailtmp))[0] if len(list(set(tailtmp))) == 1 else "Other"
                            read_multimap = int(info[7])
                            read_mismatch = 0

                            orientation = "sense" if strand == "+" else "anti" 

                            seq_len = len(seq)
                            seq_nt = seq[0]

                            aln_selector = (seq_nt, seq_len, orientation, read_multimap)

                            if all( True if a in b else True if "*" in b else False for a,b in zip(aln_selector,rule_selector) ) :
                                gene_ids.append(gene)
                                biotypes.append(biotype)
                                Classes.append(Class)
                                features.append(feature)
                                counts.append(count)
                                
                                starts.append(start)
                                ends.append(end)
                                strands.append(strand)
                                seqs.append(seq)
                                tails.append(tail)
                                rawtails.append(tailtmp)
                                
                            else : 
                                pass
                    f.close()
        q.close()
                            
        results = pd.DataFrame({
            'gene':gene_ids,
            'biotype':biotypes,
            'class':Classes,
            'feature':features,
            'count':counts,
            'tail':tails
        })
        
        bed = pd.DataFrame({
            'gene':gene_ids,
            'start':starts,
            'end':ends,
            'seq':seqs,
            'count':counts,
            'strand':strands,
            'feature':features,
            'rawtail' : rawtails,
            'tail':tails
        })
        
        if self._normalize_rpkm :
            self._bed_final['count_kmer'] = bed.apply(lambda row: row['count']/len(row['seq']), axis = 1)
            self._counts = results.groupby( ['gene', 'biotype', 'class', 'feature', 'tail'] )['count_kmer'].sum().reset_index()
        else : 
            self._bed_final = bed
            self._counts = results.groupby( ['gene', 'biotype', 'class', 'feature', 'tail'] )['count'].sum().reset_index()

        for k,v in self._normalization_features.items() : 
            if not v == 0 :
                if self._normalize_rpkm : 
                    self._counts[f'count_{k}_kmer_norm'] = self._counts.apply(lambda row : 1e6*(row['count_kmer']/v), axis = 1)
                    self._bed_final[f'count_{k}_kmer_norm'] = self._bed_final.apply(lambda row : 1e6*(row['count_kmer']/v), axis = 1)
                else : 
                    self._counts[f'count_{k}_norm'] = self._counts.apply(lambda row : 1e6*(row['count']/v), axis = 1)
                    self._bed_final[f'count_{k}_norm'] = self._bed_final.apply(lambda row : 1e6*(row['count']/v), axis = 1)

        #if self._normalization is not None : 
        #    for k,v in self._normalization_features.items() :
        #        if not v == 0 :
        #            if self._normalize_rpkm :
        #                results_grouped[f'count_{k}_kmer_norm'] = results_grouped.apply(lambda row : 1e6*(row['count_kmer']/v), axis = 1)
        #                bed[f'count_{k}_kmer_norm'] = bed.apply(lambda row : 1e6*(row['count_kmer']/v), axis = 1)
        #            else : 
        #                results_grouped[f'count_{k}_norm'] = results_grouped.apply(lambda row : 1e6*(row['count']/v), axis = 1)
        #                bed[f'count_{k}_norm'] = bed.apply(lambda row : 1e6*(row['count']/v), axis = 1)
        
        self._counts.to_csv(f"{self._outname}.counts.tsv", sep = "\t", index = False, header = True)
        self._bed_final.to_csv(f"{self._outname}.bed.tsv", sep = "\t", index = False, header = True)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument("-f", "--fasta", type = str, required = True, help="fasta file of reads")
    parser.add_argument( "-t", "--transcripts", type = str, required = True, help="tab delimited transcripts table")
    parser.add_argument("-n", "--normalization", type = str, required=False, help='normalization constants')
    parser.add_argument("-o", "--outname", type = str, required=True, help='outname prefix')
    parser.add_argument("-v", "--mismatch", type = str, required=True, help='alignment mismatch')
    parser.add_argument("-m", "--multimap", type = str, required=True, help='alignment multimap')
    parser.add_argument("-idx", "--index_path", type = str, required=True, help='path to where transcript fastas are indexed')
    parser.add_argument("-index", required=False, action='store_true', help='indicates to build index')
    parser.add_argument("-tailor", required=False, action='store_true', help='indicates to build index')
    parser.add_argument("-rpkm", action = 'store_true', required=False)

    return parser.parse_args()

def main() : 

    args = get_args()

    run = Transcripts(args.fasta, args.transcripts, args.normalization, args.outname, args.mismatch, args.multimap, args.index_path, args.index, args.rpkm, args.tailor)
    run.process()

if __name__ == "__main__" : 

    main()
                            










