#!/usr/bin/env python3

import argparse

def trim_umi(seq, p5 = 0, p3 = 0) : 
    
    p5_umi = seq[ 0 : p5 ] if p5 > 0 else ""
   
    p3_umi = seq[ len(seq)-p3 : len(seq) ] if p3 > 0 else ""

    umi_string = f"{p5_umi}~{p3_umi}"

    seq_minus_umi = seq[ 0+p5 : len(seq)-p3 ]

    return [seq_minus_umi, umi_string]

def run(inpt, output, p5, p3) : 

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
                    
                    info = trim_umi(seq, p5, p3)
                    
                    seq_xumi = info[0]
                    umi_string = info[1]

                    count = int(header[1])
                    xumi_fa += f">{seq}:{count}:{umi_string}\n{seq_xumi}\n"
                    
                    if seq_xumi in count_dict.keys() : 
                        count_dict[seq_xumi] += 1
                    else : 
                        count_dict[seq_xumi] = 1

                    header = line.strip().replace(">", "").split(":")
                    seq = ''
                    
                else : 
                    header = line.strip().replace(">", "").split(":")
                    seq = ''
            else : 
                seq += line.strip()
        else :

            info = trim_umi(seq, p5, p3)
                    
            seq_xumi = info[0]
            umi_string = info[1]

            count = int(header[1])

            xumi_fa += f">{seq}:{count}:{umi_string}\n{seq_xumi}\n"
            
            if seq_xumi in count_dict.keys() : 
                count_dict[seq_xumi] += 1
            else : 
                count_dict[seq_xumi] = 1

    f.close()

    outname_xumi_collapsed = inpt.replace(".fa", ".xumi.uniq.fa") if not output else output
    op = open(outname_xumi_collapsed, 'w')
    for k,v in count_dict.items() : 
        op.write(f">{k}:{v}\n{k}\n")
    op.close()

    outname_xumi = inpt.replace(".fa", ".xumi.fa")
    op = open(outname_xumi, 'w') 
    op.write(xumi_fa)
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