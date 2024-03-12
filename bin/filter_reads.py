#!/usr/bin/env python3

import argparse

def read_sequences_bed(bed) : 

    seqs = []
    with open(bed, 'r') as f :
        for line in f : 
            info = line.strip().split("\t")
            seq = info[3]
            seqs.append(seq)
    f.close()

    return seqs

def read_sequences_fastaq(fastq, filetype = 'fastq') : 

    if filetype == 'fastq' : 
        header = False
        entry = []
        seqs = []
        seqdict = {}
        with open(fastq, 'r') as f :
            for line in f : 
                if line.startswith("@") : 
                    if header :
                        seqs.append(str(entry[0]))
                        seqdict[str(entry[0])] = f"{header}\n{entry[0]}\n{entry[1]}\n{entry[2]}\n"
                        header = line.strip()
                        entry = []
                    else : 
                        header = line.strip()
                else : 
                    entry.append(line.strip())
            else : 
                seqdict[str(entry[0])] = f"{header}\n{entry[0]}\n{entry[1]}\n{entry[2]}\n"
                seqs.append(str(entry[0]))
        f.close()
    else : 
        header = False
        seq = ''
        seqs = []
        seqdict = {}
        with open(fastq, 'r') as f :
            for line in f : 
                if line.startswith(">") : 
                    if header :
                        seqs.append(seq)
                        seqdict[seq] = f"{header}\n{seq}\n"
                        header = line.strip()
                        seq = ''
                    else : 
                        header = line.strip()
                else : 
                    seq += line.strip()
            else :
                seqs.append(seq) 
                seqdict[seq] = f"{header}\n{seq}\n"
    f.close()


    return [seqs, seqdict]

def filter_reads(fastq, bed, outfile) : 

    out = open(outfile, 'w')

    original_seqs = read_sequences_fastq(fastq)
    aligned = read_sequences_bed(bed)
    difference = original_seqs - aligned[0]

    for s in difference : 
        out.write(aligned[1][s])
   
    out.close()

def get_args() : 

    """Parse command line parameters from input"""
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-i", "--input", type = str, required = True)
    parser.add_argument("-o", "--output", type = str, required = True)
    parser.add_argument("-b", "--bed", type = str, required = True)

    return parser.parse_args()

def main() : 

    args = get_args()
    filter_reads(args.input, args.bed, args.output)

if __name__ == "__main__" : 

    main()