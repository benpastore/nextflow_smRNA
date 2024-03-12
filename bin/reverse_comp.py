#!/usr/bin/env python3

import argparse

def rev_comp(string) : 
    
    rev = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G'}
    
    rev_comp = ''
    for i in string : 
        rev_comp += rev[i]
    
    return rev_comp[::-1]

def run(file, outname) :
    
    fq = open(outname, 'w')
    header = None
    with open(file, 'r') as f : 
        for line in f : 
            if line.startswith("@") : 
                if header : 
                    
                    seq = entries[0]
                    re_seq = rev_comp(seq)
                    fq.write(f'@{re_seq}:{header[1]}\n{re_seq}\n+\n{"I"*len(re_seq)}\n')

                    header = line.strip().replace(">", "").split(":")
                    entries = []
                    
                else : 
                    header = line.strip().replace(">", "").split(":")
                    entries = []
            else : 
                entries.append(line.strip())
        else :
            seq = entries[0]
            re_seq = rev_comp(seq)
            fq.write(f'@{re_seq}:{header[1]}\n{re_seq}\n+\n{"I"*len(re_seq)}\n')
    f.close()

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    ar = parser.add_argument_group('required arguments')
    ar.add_argument("-i", "--input", type = str, required = True)
    ar.add_argument("-o", "--output", type = str, required = True)

    return parser.parse_args()

if __name__ == "__main__" : 

    args = get_args()

    run(args.input, args.output)