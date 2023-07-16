#!/usr/bin/env python3

# trims N nucleotides from the 5' end to remove 5' barcode sequences
# This assumes file was already demultiplexed, this script is only to trim 
# the 5' barcode sequences from the reads. 

import argparse

def trim_5p_barcode(seq, p5 = 0, barcode = None) : 
    
    p5_barcode = seq[ 0 : p5 ]

    if barcode : 
        if p5_barcode == barcode : 
            seq_minus_5p_bar = seq[ 0+p5 : len(seq)]
            return seq_minus_5p_bar
        else : 
            return False
    else : 
        seq_minus_5p_bar = seq[ 0+p5 : len(seq)]
        return seq_minus_5p_bar

def run(inpt, output, p5, barcode = None) : 

    """
    Read through input fasta file, trim any 3' As

    After 3' trimming add together the scores of the same sequence
    """
    
    xumi_fa = ''
    xumi_fa_collapsed = ''

    count_dict = {}
    header = None
    with open(inpt, 'r') as f : 
        for line in f : 
            if line.startswith(">") : 
                if header : 
                    
                    info = trim_5p_barcode(seq, p5, barcode)
                    seq_x5p_bar = info
                    if seq_x5p_bar :
                        count = int(header[1])
                        
                        if seq_x5p_bar in count_dict.keys() : 
                            count_dict[seq_x5p_bar] += count
                        else : 
                            count_dict[seq_xumi] = count

                    header = line.strip().replace(">", "").split(":")
                    seq = ''
                    
                else : 
                    header = line.strip().replace(">", "").split(":")
                    seq = ''
            else : 
                seq += line.strip()
        else :
            info = trim_5p_barcode(seq, p5, barcode)
            seq_x5p_bar = info
            count = int(header[1])
            
            if seq_x5p_bar : 

                if seq_x5p_bar in count_dict.keys() : 
                    count_dict[seq_x5p_bar] += count
                else : 
                    count_dict[seq_xumi] = count

    f.close()

    barname = '' if barcode is None else f'.{barcode}'
    outname_collapsed = inpt.replace(".fa", f".x5pBar{barname}.uniq.fa") if not output else output
    op = open(outname_collapsed, 'w')
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
    ar.add_argument("-barcode", type = str, required = False, default = None)

    return parser.parse_args()

if __name__ == "__main__" : 

    args = get_args()

    run(args.input, args.output, args.p5, args.barcode)