import argparse

class Junc2bed() : 

    def __init__(self, infile, outfile) : 

        self._infile = infile
        self._outfile = outfile

    def ntm(self) : 

        self._ntms = {}
        with open(self._infile, 'r') as f : 
            for line in f : 
                info = line.strip().split("\t") 
                seq = str(info[3])

                if seq in self._ntms.keys() : 
                    self._ntms[seq] += 1
                else : 
                    self._ntms[seq] = 1
        f.close()

    def convert(self) :

        output = open(self._outfile, 'w')

        with open(self._infile, 'r') as f : 
            for line in f : 
                info = line.strip().split("\t")

                # 1-3 -> chrom, start, end
                # 4 -> seq
                # 5 -> ntm
                # 6 -> strand (all should be +)
                # 7 -> tail
                # 8 -> tail length
                # 9 -> number of locations mapped

                chrom = info[0]
                junction = chrom.split("|")
                read_start = int(info[1])
                read_end = int(info[2])
                seq = info[3]
                count = float(info[4])
                strand = info[5]
                tail = info[6]
                ntm = float(info[7])

                chrom, start, jpos, end = junction[0], int(junction[1]), junction[2].split('-'), int(junction[3])
                
                j_start = int(jpos[0])
                j_end = int(jpos[1])

                """
                start   j_start   j_end    end
                |      |          |       | 
                ########----------#########
                    G----------------->
                    |                  | 
                    read_start          read_end

                
                # entry number 1. 
                true start = start + read_start
                true end = start + j_end - 29
                """
                
                left = {}
                left['start'] = start+read_start-1
                left['end'] = j_start

                right = {}
                right['start'] = j_end
                right['end'] = j_end+read_start+len(seq)-29

                if not left['start'] > left['end'] : 
                    if not right['start'] > right['end'] : 
                        output.write(f"{chrom}\t{start+read_start-1}\t{j_start}\t{seq}\t{count}\t{strand}\t{tail}\t{ntm}\n")
                        output.write(f"{chrom}\t{j_end}\t{j_end+read_start+len(seq)-29}\t{seq}\t{count}\t{strand}\t{tail}\t{ntm}\n")

        output.close()

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument(
        "-i", 
        "--input", 
        type = str, 
        required = True
    )

    parser.add_argument(
        "-o",
        "--output",
        type = str,
        required=True
    )

    return parser.parse_args()

def main() : 

    args = get_args()

    run = Junc2bed(args.input, args.output)
    run.ntm()
    run.convert()

if __name__ == "__main__" : 

    main()

