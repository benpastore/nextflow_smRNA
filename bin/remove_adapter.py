#!/usr/bin/python3 

def trim_3(sequence, trim = False) :

    if trim : 
        return sequence[0:len(sequence) - trim]
    
    return sequence

def trim_5(sequence, trim = False) :

    if trim : 
        return sequence[trim:len(sequence)]
    
    return sequence

def split_sequence(sequence, adapter) : 

    splt = sequence.split(adapter)
    if len(splt) > 1 :        
        return splt[0]
    
    else : 
        try_adapter = adapter[0:len(adapter)-1]
        
        if len(try_adapter) >= 3 :
            return split_sequence(sequence, try_adapter)
            
        else : 
            return False

def remove_adapter(file, adapter, outname, fastq = False, keep_untrimmed = True, min_length = 15, max_length = 40, trim_3p = False, trim_5p = False, uni_fasta = True) : 

    if uni_fasta : 
        seq_count = {}

    total = 0
    with open(file, 'r') as f : 
        for line in f : 
            total += 1
    f.close()

    if fastq : 
        n_line_per_entry = 4 # 4 lines per fastq entry ...
        total = total / 4
    else : 
        n_line_per_entry = 2 # 2 lines per fasta entry ...
        total = total / 2



    trimmed = 0
    missing = 0
    outlines = ''
    lines = []
    with open(file, 'r') as f : 
        for line in f :

            info = line.strip()
            lines.append(info)

            if len(lines) == n_line_per_entry :
                raw_seq = lines[1]
                
                seq = split_sequence(raw_seq, adapter)
                if seq == False : 
                    missing += 1
                    if keep_untrimmed : 
                        seq = raw_seq
                    else : 
                        lines = []
                        continue
                else : 
                    trimmed += 1

                if trim_3p : 
                    seq = trim_3(seq, trim_3p) 

                if trim_5p :
                    seq = trim_5(seq, trim_5p)
                
                if len(seq) >= min_length and len(seq) <= 40 :
                    if fastq : 
                        outlines += f"{lines[0]}\n{seq}\n{lines[2]}\n{lines[3]}\n"
                    else : 
                        outlines += f"{lines[0]}\n{seq}\n"

                    if uni_fasta : 
                        if not seq in seq_count.keys() : 
                            seq_count[seq] = 1
                        else : 
                            seq_count[seq] += 1

                lines = []
    f.close()

    output = open(outname, 'w')
    output.write(outlines)
    output.close()

    if uni_fasta :
        uni_lines = ''
        for k,v in seq_count.items() : 
            uni_lines += f">{k}:{v}\n{k}\n"
    
    uni_fasta_name = outname.replace(".fa", ".uni.fa").replace(".fq", ".uni.fa").replace(".fastq", ".uni.fa")
    uni_output = open(uni_fasta_name, 'w')
    uni_output.write(uni_lines)
    uni_output.close()

    print(f"Total number of reads : {total}")
    print(f"Total number of reads with adapter : {trimmed} ({round(100*(trimmed/total), 3)} %)")


def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    ar = parser.add_argument_group('required arguments')
    ar.add_argument("-i", "--input", type = str, required = True)
    ar.add_argument("-o", "--output", type = str, required = True)
    ar.add_argument("-a", "--adapter", type = str, required = True)
    ar.add_argument("-m", "--min_length", type = int, required = False)
    ar.add_argument("-M", "--max_length", type = int, required = False)
    ar.add_argument("-t", "--trim_3p", type = int, required = False)
    ar.add_argument("-T", "--trim_5p", type = int, required = False)
    ar.add_argument("-u", "--uni_fasta", action = 'store_true', required = False)
    ar.add_argument("-k", "--keep_untrimmed", action = 'store_true', required = False)
    ar.add_argument("-q", "--fastq", action = 'store_true', required = False)

    return parser.parse_args()

def main() : 

    args = get_args()

    if not args.min_length : 
        m = 15
    else : 
        m = args.min_length
    
    if not args.max_length : 
        M = 40
    else : 
        m = args.max_length

    remove_adapter(
        args.file, args.adapter, args.outname, fastq = args.fastq, keep_untrimmed = args.keep_untrimmed, 
        min_length = m, max_length = M, trim_3p = args.trim_3p, trim_5p = args.trim_5p, uni_fasta = args.uni_fasta
        )

if __name__ == "__main__" : 

    main()