#!/usr/bin/env python3

import argparse

def trim(seq, p5 = 0, p3 = 0) :
    
    trim_seq = seq[ 0+p5 : len(seq)-p3 ]

    return trim_seq

def run(inpt, output, p5, p3) : 

    """
    Read through input fasta file, hard trim nucleotides from eith 5' or 3' end

    After timming add together sequence counts for identical sequences
    """
    
    count_dict = {}
    header = None
    with open(inpt, 'r') as f : 
        for line in f : 
            if line.startswith(">") : 
                if header : 
                    
                    info = trim(seq, p5, p3)
                    
                    trim_seq = info

                    count = int(header[1])
                    
                    if trim_seq in count_dict.keys() : 
                        count_dict[trim_seq] += count
                    else : 
                        count_dict[trim_seq] = count

                    header = line.strip().replace(">", "").split(":")
                    seq = ''
                    
                else : 
                    header = line.strip().replace(">", "").split(":")
                    seq = ''
            else : 
                seq += line.strip()
        else :
            
            info = trim(seq, p5, p3)
                    
            trim_seq = info

            count = int(header[1])
            
            if trim_seq in count_dict.keys() : 
                count_dict[trim_seq] += count
            else : 
                count_dict[trim_seq] = count

    f.close()

    outname = inpt.replace(".fa", f".x5p{p5}bp.x3p{p3}.fa") if not output else output
    op = open(outname, 'w')
    for k,v in count_dict.items() : 
        op.write(f">{k}:{v}\n{k}\n")
    op.close()

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    ar = parser.add_argument_group('required arguments')
    ar.add_argument("-i", "--input", type = str, required = True)
    ar.add_argument("-o", "--output", type = str, required = False)
    ar.add_argument("-p5", type = int, required = False, default = 0)
    ar.add_argument("-p3", type = int, required = False, default = 0)

    return parser.parse_args()

if __name__ == "__main__" : 

    args = get_args()

    run(args.input, args.output, args.p5, args.p3)