
import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO

def filter_artifact(tailor_bed, fasta, N, outfile) :

    """
    A function to filter reads where the length of the tail is not equal to the edit distance to the reference seq.

    **Also I will filter sequences based on 4 preceeding nucleotides before tail to remove possible sequencing errors. 
    (must be the same NT repeated N times)

    For example a 3 nt tail where there is only 1 mismatch.

    Method : iterate through bed file, for each bed entry compare the aligned seq to the reference. 
    If the number of mismatches to the reference equal to the length of the tail, keep, otherwise do not keep.

    Input : tailor bed file, output of tailor_sam_to_bed

    Output : cleaned bed file
    
    """
    
    if os.path.exists(outfile) : 
        os.remove(outfile)

    out = open(outfile, 'a')

    records = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))

    with open(tailor_bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")

            chrom  = info[0]
            start = int(info[1])
            end = int(info[2])
            seq = str(info[3])
            tail_length = int(info[7])

            seq_minus_tail = seq[:-tail_length]
            lastN = seq_minus_tail[-int(N):]

            if tail_length < 2 : 
                out.write(line)
            
            elif tail_length >= 2 and tail_length <= 4 : 
    
                ref = str(records[chrom][start : end + tail_length].seq)

                edit_distance = sum(1 for a, b in zip(seq, ref) if a != b)

                if edit_distance == tail_length :
                    if not len(set(lastN)) == 1 : 
                        out.write(line)
            
            elif tail_length > 4 : 

                # if reference is AATAA
                # and tail is     TTTTT
                # edit_distnace = 4
                # tail_length = 5
                # 5 - 4 = 1
                # the edit distance would be 4 but the tail length would be 5, 
                # yet it is likely that TTTTT is a real tail, even if the edit 
                # distance is not equal to the 

                # so we need function to say that if the difference between the edit distance and tail length 
                
                # The chance of a nucleotide is in a give position is 1/4,
                # thus I will allow a differnece in (tail_length - edit_distance) of 
                # less than or equal to the tail_length*(1/4), rounded to nearest whole number
                
                # also if the single 

                error_rate = 0.25
                allowance = round(tail_length*error_rate)

                ref = str(records[chrom][start : end + tail_length].seq)

                edit_distance = sum(1 for a, b in zip(seq, ref) if a != b)

                difference = tail_length - edit_distance

                if difference == 0 : 
                    out.write(line)

                elif difference < allowance :
                    if not len(set(lastN)) == 1 :
                        out.write(line)


def main() : 

    bed = sys.argv[1]
    fasta = sys.argv[2]
    outfile = sys.argv[3]

    N = 4

    filter_artifact(bed, fasta, N, outfile)

if __name__ == "__main__" : 

    main()