#!/usr/bin/python

import sys
import argparse

def calc_frequency(seq, count) : 
    
    count = count_characters(seq)
    
    most_frequent = max(count.values())
    
    frequency = ( most_frequent / len(seq) ) 

    return frequency

def count_characters(string) : 

    counts = {}
    counts['A'] = 0
    counts['G'] = 0
    counts['C'] = 0
    counts['T'] = 0
    counts['N'] = 0

    for i in string : 
        if i not in counts.keys() : # avoid non C, T, G, A, N characters
            return False
        counts[i] += 1

    return counts

def filter(file, output, threshold, fastq = False) :

    i = 0
    with open(file, 'r') as f : 
        for line in f : 
            i += 1
    f.close()

    if fastq : 
        n_line_per_entry = 4 # 4 lines per fastq entry ...
    else : 
        n_line_per_entry = 2 # 2 lines per fasta entry ...

    j = 0
    outlines = ''
    lines = []
    with open(file, 'r') as f : 
        for line in f : 
            info = line.strip()
            lines.append(info)

            if len(lines) == n_line_per_entry :
                seq = str(lines[1])
                count = count_characters(seq)

                if count : 
                    freq = calc_frequency(seq, count)
                    if freq < threshold : 
                        j += 1
                        if fastq : 
                            outlines += f"{lines[0]}\n{lines[1]}\n{lines[2]}\n{lines[3]}\n"
                        else : 
                            outlines += f"{lines[0]}\n{lines[1]}\n"
                lines = []
    f.close() 

    ou = open(output, 'w')
    ou.write(outlines)
    ou.close()

    print(f"Total entries = {i/4 if fastq else i/2}\nTotal entries passing filter {j}\nPercent entries filtered { 100 - 100*( j / (i/4 if fastq else i/2) ) }\n")

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    ar = parser.add_argument_group('required arguments')
    ar.add_argument("-i", "--input", type = str, required = True)
    ar.add_argument("-o", "--output", type = str, required = True)
    ar.add_argument("-t", "--threshold", type = float, required = False)
    ar.add_argument("-q", "--fastq", action = 'store_true', required = False)

    return parser.parse_args()

def main() : 

    args = get_args()

    if not args.threshold : 
        threshold = 0.70
    else : 
        threshold = args.threshold
    
    filter(file = args.input, output = args.output, threshold = threshold, fastq = args.fastq)

if __name__ == "__main__" : 

    main()