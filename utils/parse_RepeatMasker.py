import subprocess
import sys
import os

""" 
This script will work with WormBase gff3 files ONLY!!! 

The bed format needed for the pipeline is

1. chrom
2. start (base 0)
3. end
4. source (aka WormBase)
5. feature (gene, exon, utr)
6. strand
7. gene name
8. sequence name
9. public name / common name / locus (IF THIS IS A TRANSPOSON THIS FEILD IS USED FOR TRANSPOSON FAMILY)
10. gene biotype
11. class / pathways 

** If the gene is not annotated and only has a gene name (aka no locus name), make the locus name the gene name or sequence name
** If the gene is not annotated as beloning to class or biotype make the class or biotype unknown

"""


class RepeatMasker() : 
    
    def __init__(self, gff3, outprefix) : 
        self._gff3 = gff3
        self._outprefix = outprefix
        self._outfile_gff3 = f"{outprefix}.gff3"
        self._outfile_bed = f"{outprefix}.bed"
        
        if os.path.exists(self._outfile_gff3) : 
            os.remove(self._outfile_gff3)
        
        if os.path.exists(self._outfile_bed) : 
            os.remove(self._outfile_bed)
        
        self._out = open(self._outfile_gff3, 'w')
    
    def slim(self) : 

        if not os.path.exists(f"{self._outprefix}.RepeatMasker.gff3") : 
            if "*.gz" in self._gff3 : 
                cmd = f"""zcat {self._gff3} | awk '(match($2,"RepeatMasker"))' > {self._outprefix}.RepeatMasker.slim.gff3"""
                os.system(cmd)
            else : 
                cmd = f"""cat {self._gff3} | awk '(match($2,"RepeatMasker"))' > {self._outprefix}.RepeatMasker.slim.gff3"""
                os.system(cmd)
        else : 
            pass
            
        self._gff3 = f"{self._outprefix}.RepeatMasker.slim.gff3"
                    
    def parse_genes(self) : 
        
        self._header = ""
        self._attributes = {}
        self._biotypes = {}
        self._lines = ''
        
        with open(self._gff3, 'r') as f :
            for line in f :
                if not line.startswith("#") :
                    info = line.strip().split("\t")

                    source = str(info[1])
                    feature = str(info[2])

                    if source == "RepeatMasker" : 

                        attrs = dict(x.split("=") for x in  info[8].split(";"))

                        target = attrs['Target'].replace(" ","_")
                        sequence_name = target
                        locus = target
            
                        biotype = "RepeatMasker"
                       
                        # for each feature type have a count of biotypes
                        entry = (feature, biotype)
                        if entry not in self._biotypes.keys() :
                            self._biotypes[entry] = 1
                        else : 
                            self._biotypes[entry] += 1

                        classes_joined = 'repeat'

                        attrs_string = f"ID={target};Name={target};sequence_name={target};locus={target};biotype={biotype};Class={classes_joined}"
                        
                        self._lines += f"{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}\t{info[4]}\t{info[5]}\t{info[6]}\t{info[7]}\t{attrs_string}\n"

        f.close()

        for k,v in self._biotypes.items() : 
            self._header += f"# {k[0]} -> {k[1]} : ({v})\n"
        
        self._out.write(self._header)
        self._out.write(self._lines)

    def gff3_to_bed(self) : 
        
        out = open(self._outfile_bed, 'w')
        with open(self._outfile_gff3, 'r') as f : 
            for line in f : 
                if not line.startswith("#") :
                    info = line.strip().split("\t")
                    chrom = info[0]
                    start = int(info[3]) - 1
                    end = int(info[4])
                    strand = str(info[6])
                    
                    source = info[1]
                    feature = info[2]
                    
                    attrs = dict(x.split("=") for x in info[8].split(";"))
                    ID = attrs['ID']
                    name = attrs['Name']
                    seq_name = attrs['sequence_name']
                    locus = attrs['locus']
                    Class = attrs['Class']
                    biotype = attrs['biotype']
                    
                    out.write(f"{chrom}\t{start}\t{end}\t{source}\t{feature}\t{strand}\t{name}\t{seq_name}\t{locus}\t{biotype}\t{Class}\n") 

        f.close()
        out.close()

def run(gff, outfile) :

    x = RepeatMasker(gff, outfile)

    print(f">{gff}\n>{outfile}")

    print("Slim down gff3 file")
    x.slim()

    print("Parse genes")
    x.parse_genes()

    print("Convert to bed")
    x.gff3_to_bed()

def main() : 

    # celegans
    gff = "/fs/ess/PCON0160/ben/genomes/ce/WS279/c_elegans.PRJNA13758.WS279.annotations.gff3"
    outfile = "/fs/ess/PCON0160/ben/genomes/ce/WS279/ce_WS279_repeat_maskter"
    run(gff, outfile)

    # cbriggsae
    gff = "/fs/ess/PCON0160/ben/genomes/cb/WS279/c_briggsae.PRJNA10731.WS279.annotations.gff3"
    outfile = "/fs/ess/PCON0160/ben/genomes/cb/WS279/cb_WS279_repeat_masker"
    run(gff, outfile)

if __name__ == "__main__" : 

    main()

##########