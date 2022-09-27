
import sys
import os 

def process_bed(junction_bed, ref_length, outfile) : 

    """
    A function to convert a bedfile from junciton alignment to standard bedfile 

    Input : sam2bed output of junction.sam

    Output : standard bedfile, each alignment will get split into 2 bed entries
    
    """

    if os.path.exists(outfile) : 
        os.remove(outfile)

    out = open(outfile, 'a')

    with open(junction_bed, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            junction = info[0].split("|")

            chrom = junction[0]
            read_start = int(info[1])
            read_end = int(info[2])

            # split the junction 
            junc_right = {}
            junc_left = {}
            
            junc_right['start'] = int(junction[1])
            junc_right['end'] = int(junction[2].split("-")[0])

            junc_left['start'] = int(junction[2].split("-")[1])
            junc_left['end'] = int(junction[3])

            # split the alignment 
            right = {}
            left = {}

            right['start'] = junc_right['start'] + read_start
            right['end'] = junc_right['end']

            if right['start'] == right['end'] : 
                right['start'] = right['start'] - 1

            left['start'] = junc_left['start']
            left['end'] = junc_left['end'] - (ref_length - read_end)

            if left['start'] == left['end'] : 
                left['end'] = left['end'] + 1
            
            # deal with counting issue
            count = int(info[4])
            split_count = int(round(count/2)) if round(count/2) > 0 else 1

            if count%2 == 0 : 
                right_count = split_count
                left_count = split_count
            elif count == 1 : 
                right_count = split_count
                left_count = split_count
            else :
                if (right['end'] - right['start']) > (left['end'] - left['start']) : 
                    right_count = split_count
                    left_count = split_count - 1
                else :
                    right_count = split_count - 1
                    left_count = split_count
                    


            # require that the read spans the junction, split the alignment read count in half, round to nearest whole number
            if not right['start'] > right['end'] : 
                if not left['start'] > left['end'] :
                    entry1 = "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        chrom,
                        right['start'],
                        right['end'],
                        info[3],
                        right_count,
                        info[5]
                    )

                    entry2 = "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        chrom,
                        left['start'],
                        left['end'],
                        info[3],
                        left_count,
                        info[5]
                    )

                    out.write(entry1)
                    out.write(entry2)

    f.close()
    out.close() 

def main() : 

    junction_bed = sys.argv[1]
    outfile = sys.argv[2]
    ref_length = 58

    process_bed(junction_bed, ref_length, outfile)

if __name__ == "__main__" : 

    main()






            


