#!/usr/bin/env python3

import argparse

def trim_a(seq) : 
    
    rev_seq = seq[::-1]
    remove = 0
    for i,nt in enumerate(rev_seq) : 
        if nt == "A" : 
            remove += 1
        else :
            out = rev_seq[remove:len(rev_seq)][::-1]
            return out

def run(input, output) : 

    """
    Read through input fasta file, trim any 3' As

    After 3' trimming add together the scores of the same sequence
    """

    count_dict = {}
    header = None
    with open(input, 'r') as f : 
        for line in f : 
            if line.startswith(">") : 
                if header : 
                    
                    trimmed = trim_a(seq)
                    count = int(header[1])
                    
                    if trimmed in count_dict.keys() : 
                        count_dict[trimmed] += count
                    else : 
                        count_dict[trimmed] = count

                    header = line.strip().replace(">", "").split(":")
                    seq = ''
                    
                else : 
                    header = line.strip().replace(">", "").split(":")
                    seq = ''
            else : 
                seq += line.strip()
        else : 
            trimmed = trim_a(seq)
            count = int(header[1])
            
            if trimmed in count_dict.keys() : 
                count_dict[trimmed] += count
            else : 
                count_dict[trimmed] = count
    f.close()

    # output fasta of trimmed A's....
    op = open(output, 'w')
    for k,v in count_dict.items() : 
        op.write(f">{k}:{v}\n{k}\n")
    op.close()

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