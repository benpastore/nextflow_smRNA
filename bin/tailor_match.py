import sys
import os
import pandas as pd
import itertools
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

def filter_tails(sequence, tail, tail_length, chrom, start, end, strand, records, N) :

    """
    A function to filter reads where the length of the tail is not equal to the edit distance to the reference seq.

    Also I will filter sequences based on 4 preceeding nucleotides before tail to remove possible sequencing errors. 
    (must be the same NT repeated N times)

    Exceptions : 
    1. If the tail length is <= 2 keep it 
    2. If the tail is homogenously one nucleotide aka TTTTTTT or AAAAAAAA keep it 
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    seq_minus_tail = sequence[:-int(tail_length)]
    lastN = seq_minus_tail[-int(N):]

    if len(set(lastN)) == 1 : 
        return False

    if tail_length <= 2 :
        return True
    elif len(set(tail)) == 1 : 
        return True
    else : 
        # find genomic sequence, take rev comp if on minus strand
        if strand == "+" :
            ref = str(records[chrom][start : end + tail_length].seq)
        else : 
            fwd = records[chrom][start : end + tail_length].seq
            seq_obj = Seq(fwd)
            ref = seq_obj.reverse_complement()
            
        # if tail length(x) 4 >= x >= 3
        if tail_length >= 3 and tail_length <= 4 : 

            edit_distance = sum(1 for a, b in zip(sequence, ref) if a != b)

            if not edit_distance == tail_length :
                return False

        # if longer than 4
        elif tail_length > 4 : 

            error_rate = 0.25
            allowance = round(tail_length*error_rate)

            ref = str(records[chrom][start : end + tail_length].seq)
            edit_distance = sum(1 for a, b in zip(sequence, ref) if a != b)
            difference = tail_length - edit_distance
            if not difference <= allowance :
                return False
                
        else : 
            return False

    return True


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