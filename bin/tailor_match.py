import sys
import os
import pandas as pd
import itertools
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from alternating_string import is_alternating

def filter_tails(sequence, tail, chrom, start, end, strand, records, N = False) :

    """
    A function to filter reads where the length of the tail is not equal to the edit distance to the reference seq.

    Also I will filter sequences based on 4 preceeding nucleotides before tail to remove possible sequencing errors. 
    (must be the same NT repeated N times)

    Exceptions : 
    1. If the tail length is <= 2 keep it 
    2. If the tail is homogenously one nucleotide aka TTTTTTT or AAAAAAAA keep it 
    """

    tail_length = len(tail)

    if N : 
        N = N
    else : 
        N = 4

    #print(f"Sequence: {sequence}")
    #print(f"Tail: {tail}")

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    seq_minus_tail = sequence[:-int(tail_length)]
    
    lastN = seq_minus_tail[-int(N):]
    #print(f"Last {N} nucleotides: {lastN}")

    if len(set(lastN)) == 1 :
        #print(f"Last 4 nucleodies are the same ({lastN})....fail")
        return False

    # len(set(tail)) == 1 or 
    if tail_length == 1 or is_alternating(tail) : 
        return True
    
    else :
        if strand == "+" :
            ref = str(records[chrom][start : end + tail_length].seq)
        else : 
            fwd = records[chrom][start - tail_length : end ].seq
            seq_obj = Seq(fwd)
            ref = str(seq_obj.reverse_complement())
        
        edit_distance = sum(1 for a, b in zip(sequence, ref) if a != b)

        if len(tail) >= 4 and len(set(tail)) == 1 : 
            return True 

        if edit_distance == tail_length :
            #print(f"Edit distance == tail length {sequence} (seq) vs. {ref} (ref) with tail {tail}..pass")
            return True
            
        else : 
            #print(f"Edit distance does not equal tail length: {sequence} (seq) vs. {ref} (ref) with tail {tail}\ncheck for alternating tail...")
            return False

    #print("Other tail not meeting some specification...fail")
    return False

def parse_normalization_constants(file) : 

    constants = {}
    with open(file, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            key = info[0]
            value = float(info[1])
            constants[key] = value
    f.close() 
    return constants
    
#######################
# if longer than 4
#elif tail_length >= 5 : 

#    error_rate = 0.25
    
#    allowance = round(tail_length*error_rate)

#    ref = str(records[chrom][start : end + tail_length].seq)
#    
#    edit_distance = sum(1 for a, b in zip(sequence, ref) if a != b)
#    
#    difference = tail_length - edit_distance
#    
#    if not difference <= allowance :
#        
#        return False
#        
#else : 
#    return False

    
    