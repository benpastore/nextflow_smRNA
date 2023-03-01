#!/usr/bin/env python3

import pandas as pd

def filter_21U(ntm) : 

    ntag = 0
    chroms, starts, ends, seqs, counts, strands = [],[],[],[],[],[]
    with open(ntm, 'r') as f : 
        for line in f : 
            info = line.strip().split("\t")
            chrom = info[0]
            start = info[1]
            end = info[2]
            seq = str(info[3])
            count = float(info[4])
            strand = info[5]
            ntag += count
    f.close() 

    bed = pd.DataFrame({
        'chrom':chroms,
        'start':starts,
        'end':ends,
        'seq':seqs,
        'count':counts,
        'strand':strands
    })

    bed['rpm'] = bed.apply(lambda row: 1e6*(row['count']/ntag), axis = 1)
    bed = bed[['chrom','start','end','seq','rpm','strand']]

    # separate into sense / antisense (watson / crick) based on aligned strand
    sense = bed[ bed['strand'] == "+" ].reset_index(drop = True)
    anti = bed[ bed['strand'] == "-"].reset_index(drop == True)

    # aggregate 5' ends 
    sense_5p_agg = sense.groupby(by = ['chrom', 'start', 'strand'] )['rpm'].sum().reset_index()
    anti_5p_agg = anti.groupby(by = ['chrom', 'end', 'strand'] )['rpm'].sum().reset_index()
    
    # find most abundant end
    sense_3p = sense.sort_values('rpm', ascending=False).drop_duplicates(['chrom','start']).reset_index(drop = True)
    sense_3p = sense_3p[['chrom','start','end','seq']]
    
    anti_3p = anti.sort_values('rpm', ascending = False).drop_duplicates(['chrom','end']).reset_index(drop = True)
    anti_3p = anti_3p[['chrom','start','end','seq']]

    # merge 
    sense_collapsed = sense_5p_agg.merge(sense_3p, on = ['chrom','start'], how = 'outer')
    anti_collapsed = anti_5p_agg.merge(anti_3p, on = ['chrom','end'], how = 'outer')
    
    collapsed_df = sense_collapsed.append(anti_collapsed)
    collapsed_df = collapsed_df[['chrom','start','end','rpm','strand']]
    collapsed_df = collapsed_df[ (len(collapsed_df['seq']) == 21) & (collapsed_df['seq'][0] == "U") ]


def get_opts() : 
    
    parser = argparse.ArgumentParser(add_help = True)
    
    parser.add_argument("-i", "--input", type = str, required = True, help = "Directory with files")

    return parser.parse_args()

def main() :

    args = get_opts()
    ntm = args.input
    res = filter_21U(ntm)
    res.to_csv(ntm.replace(".bed",".21U.bed"), sep = "\t", header = False, index = False)
    
if __name__ == "__main__" : 

    main()