import sys
import os
    

def bed_to_ntm(bed, outfile) : 
    """
    Convert a bed file to a number of times mapped file. 
    For each read divide the raw read counts by the number of locations in the genome where the read maps.

    Inputs : 
        - bed file
    
    Outputs : 
        - ntm dataframe 
    """
    
    ntms = {}
    with open(bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            seq = info[3]
            if seq in ntms.keys() : 
                ntms[seq] += 1
            else : 
                ntms[seq] = 1
    f.close()

    out = open(outfile, 'w')
    with open(bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            nfeilds = len(info)
            chrom = info[0]
            start = info[1]
            end = info[2]
            seq = info[3]
            ntm = float(info[4]) / ntms[seq]
            strand = info[5]
            mismatch = info[6]
            out.write(f"{chrom}\t{start}\t{end}\t{seq}\t{ntm}\t{strand}\t{ntms[seq]}\t{mismatch}\n")

    f.close()
    out.close()

def main() : 

    bed = sys.argv[1]
    outfile = sys.argv[2]

    bed_to_ntm(bed, outfile)

if __name__ == "__main__" : 

    main()

















