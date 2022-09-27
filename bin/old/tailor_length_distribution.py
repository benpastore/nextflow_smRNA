
import sys
import os

def length_distribution(tailed_bed, perfect_bed, outfile) :

    """
    A function to compute the length distribution of tailed reads 

    Input : tailed bed file

    Output : file with length distribution (length / tail / rpm / ratio)
    
    """

    if os.path.exists(outfile) : 
        os.remove(outfile)

    out_perfect = open(f"{outfile}_perfect_length_distribution.tsv", 'a') 
    out_tailed = open(f"{outfile}_tailed_length_distribution.tsv", 'a') 

    # sum perfect
    perfect = 0
    with open(perfect_bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            rpm = float(info[4])
            perfect += rpm
    f.close()

    # sum tailed 
    tailed = 0 
    with open(tailed_bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            rpm = float(info[4])
            tailed += rpm
    f.close()

    total = perfect + tailed

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

    # tailed 
    tailed_len_dis = {}
    with open(tailed_bed,'r') as f : 
        for line in f : 
            info = line.strip().split("\t") 
            rpm = float(info[4])
            seq = str(info[3])
            seqlen = len(seq)

            if len( set(str(info[6])) ) == 1 :
                tail = list(set(str(info[6])))[0]
            else : 
                tail = "other"

            key = f"{seqlen}\t{tail}\t"

            if not key in tailed_len_dis.keys() : 
                tailed_len_dis[key] = rpm 
            else : 
                tailed_len_dis[key] += rpm
    f.close()

    for key, value in tailed_len_dis.items() : 
        out_tailed.write(f"{key}\t{value}\t{100*(value/total)}\n")

    out_tailed.close()

def main() : 

    tailed_bed = sys.argv[1]
    perfect_bed = sys.argv[2]
    outfile = sys.argv[3]

    length_distribution(tailed_bed, perfect_bed, outfile)

if __name__ == "__main__" : 

    main()




                





