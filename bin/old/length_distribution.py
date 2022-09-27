
import sys
import os

def length_distribution(perfect_bed, outfile) :

    """
    A function to compute the length distribution of tailed reads 

    Input : tailed bed file

    Output : file with length distribution (length / tail / rpm / ratio)
    
    """

    if os.path.exists(outfile) : 
        os.remove(outfile)

    out_perfect = open(f"{outfile}_perfect_length_distribution.tsv", 'a') 

    # sum perfect
    perfect = 0
    with open(perfect_bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            rpm = float(info[4])
            perfect += rpm
    f.close()

    total = perfect

    # perfect
    perfect_len_dis = {}
    with open(perfect_bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            rpm = float(info[4])
            seq = str(info[3])
            seq_len = len(seq)

            key = f"{seq_len}\t*"
            if not key in perfect_len_dis.keys() : 
                perfect_len_dis[key] = rpm 
            else :
                perfect_len_dis[key] += rpm
    f.close()

    for key, value in perfect_len_dis.items() : 
        out_perfect.write(f"{key}\t{value}\t{100*(value/total)}\n")

def main() : 

    perfect_bed = sys.argv[1]
    outfile = sys.argv[2]

    length_distribution(perfect_bed, outfile)

if __name__ == "__main__" : 

    main()




                





