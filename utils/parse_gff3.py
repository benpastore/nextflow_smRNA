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

###
def process_exons_utr(gff3, exon_gene_assoc, attributes_info, out, FEATURE, SOURCE) :

    if os.path.exists("tmp") :
        os.remove("tmp")

    bed = open("tmp", 'a')
    nlines = 0
    with open(gff3, 'r') as f :
        for line in f : 
            if not line.startswith("#") :
                info = line.strip().split("\t")

                source = str(info[1])
                feature = str(info[2])
                start = int(info[3]) - 1
                end = int(info[4])
                strand = str(info[6])
                chrom = info[0]

                if (feature == FEATURE and source == SOURCE) : 

                    attrs = dict(x.split("=") for x in  info[8].split(";"))
                    exon_string = attrs['Parent'].replace("Pseudogene:","").replace("Transcript:","").replace("transcript:","").split(".")

                    if "CB" in exon_string[0] : 
                        exon = exon_string[0]
                        exon_ID = exon[:-1] if exon[-1].isalpha() else exon
                    elif "ci" in gff3 : 
                        exon = f"{exon_string[0]}"
                        exon_ID = exon
                    else : 
                        exon = f"{exon_string[0]}.{exon_string[1]}"
                        exon_ID = exon[:-1] if exon[-1].isalpha() else exon

                    bed.write(f"{exon_ID}\t{start}\t{end}\t{chrom}\t.\t{strand}\n")
                    nlines += 1
    f.close()
    bed.close()

    if nlines > 0 : 
        cmd = "cat tmp | sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,6 -o distinct,distinct > merge"
        os.system(cmd)

        with open("merge", 'r') as f : 
            for line in f : 
                info = line.strip().split("\t")
                transcript = info[0].replace("transcript:", "")
                start = int(info[1]) + 1
                end = int(info[2])
                chrom = info[3]
                strand = info[4]
                try : 
                    gene = exon_gene_assoc[transcript]
                    attributes = attributes_info[gene]
                    out.write(f"{chrom}\t{SOURCE}\t{FEATURE}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
                except : 
                    attributes = f"ID={transcript};Name={transcript};sequence_name={transcript};locus={transcript};biotype=unknown;Class=unknown"
                    out.write(f"{chrom}\t{SOURCE}\t{FEATURE}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
        f.close()
    
        os.remove("merge")
        os.remove("tmp")


class process_gff3() : 
    
    def __init__(self, gff3, pathway_lists, outprefix) : 
        self._gff3 = gff3
        self._outprefix = outprefix
        
        self._pathway_list = []
        if pathway_lists is not None : 
            with open(pathway_lists, 'r') as f : 
                for line in f : 
                    info = line.strip().split("\t")
                    clas = info[0]
                    path = info[1]
                    mes = info[2]
                    
                    cur = [clas, path, mes]
                    self._pathway_list.append(cur)
            f.close() 
        else : 
            self._pathway_list = []
            
        self._outfile_gff3 = f"{outprefix}.gff3"
        self._outfile_bed = f"{outprefix}.bed"
        
        if os.path.exists(self._outfile_gff3) : 
            os.remove(self._outfile_gff3)
        
        if os.path.exists(self._outfile_bed) : 
            os.remove(self._outfile_bed)
        
        self._out = open(self._outfile_gff3, 'w')
    
    def slim(self) : 

        if not os.path.exists(f"{self._outprefix}.slim.gff3") : 
            if ".gz" in self._gff3 : 
                print("Gzipped")
                cmd = f"""zcat {self._gff3} | awk '(match($2,"WormBase") || match($2,"WormBase_imported"))' | sed -e 's/WormBase_imported/WormBase/g' > {self._outprefix}.slim.gff3"""
                os.system(cmd)
            else : 
                cmd = f"""cat {self._gff3} | awk '(match($2,"WormBase") || match($2,"WormBase_imported"))' | sed -e 's/WormBase_imported/WormBase/g' > {self._outprefix}.slim.gff3"""
                os.system(cmd)
        else : 
            pass
            
        self._gff3 = f"{self._outprefix}.slim.gff3"

    def parse_pathways(self) :
                
        self._header = '##GFF3\n'
        if len(self._pathway_list) > 0 :
        
            self._pathways = {}
            for entry in self._pathway_list : 
                
                pathway = entry[0]
                file = entry[1]
                message = entry[2]
                
                genes = []
                with open(file, 'r') as f :
                    for line in f : 
                        info = line.strip().split("\t") 
                        genes.append(str(info[0]))
                f.close() 

                self._header += f"# {message} ({len(genes)})\n"
                if pathway in self._pathways.keys() : 
                    self._pathways[pathway] += genes
                else : 
                    self._pathways[pathway] = genes
                    
        else : 
            self._pathways = {}
                    
    def parse_genes(self) : 
        
        self._attributes = {}
        self._biotypes = {}
        self._lines = ''
        
        with open(self._gff3, 'r') as f :
            for line in f :
                if not line.startswith("#") :
                    info = line.strip().split("\t")

                    source = str(info[1])
                    feature = str(info[2])

                    if (feature == "gene" and source == "WormBase") or (feature == "transposable_element" and source == "WormBase_transposon" and "Predicted" not in info[8]) or (feature == "miRNA" and source == "WormBase") or (feature == "pre_miRNA" and source == "WormBase") : 

                        attrs = dict(x.split("=") for x in  info[8].split(";"))

                        id = attrs['ID'].replace("gene:", "") if 'ID' in attrs.keys() else 'unknown'
                        name = attrs['Name'] if 'Name' in attrs.keys() else 'unknown'
                        sequence_name = attrs['sequence_name'] if 'sequence_name' in attrs.keys() else attrs['Name']
                        locus = attrs['locus'] if 'locus' in attrs.keys() else attrs['Name']

                        if feature == "transposable_element" : 
                            family = attrs['Family']
                        
                        if feature == 'miRNA' : 
                            if "matureType" in attrs.keys() : 
                                locus = attrs['matureType'] 
                            else : 
                                pass

                        if 'biotype' in attrs.keys() :
                            biotype = attrs['biotype']
                        elif feature == "transposable_element" : 
                            biotype = "transposable_element"
                        elif feature == "miRNA" : 
                            biotype = "miRNA"
                        elif feature == "pre_miRNA" :
                            biotype = "pre_miRNA"
                        else : 
                            biotype = "unknown"

                        # for each feature type have a count of biotypes
                        entry = (feature, biotype)
                        if entry not in self._biotypes.keys() :
                            self._biotypes[entry] = 1
                        else : 
                            self._biotypes[entry] += 1

                        # assign pathways
                        if len(self._pathways) > 0 and self._pathways is not None:
                            classes = []
                            for k,v in self._pathways.items() : 
                                if name in v : 
                                    classes.append(k)
                                elif sequence_name in v : 
                                    classes.append(k)
                                elif locus in v : 
                                    classes.append(k)
                        else : 
                            classes = []

                        classes_joined = 'none' if len(classes) == 0 else ",".join(classes)

                        if not feature == "transposable_element" : 
                            attrs_string = f"ID={id};Name={name};sequence_name={sequence_name};locus={locus};biotype={biotype};Class={classes_joined}"
                        else : 
                            attrs_string = f"ID={id};Name={name};sequence_name={sequence_name};locus={locus};Family={family};biotype={biotype};Class={classes_joined}"
                        self._attributes[name] = attrs_string

                        self._lines += f"{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}\t{info[4]}\t{info[5]}\t{info[6]}\t{info[7]}\t{attrs_string}\n"

        f.close()

        for k,v in self._biotypes.items() : 
            self._header += f"# {k[0]} -> {k[1]} : ({v})\n"
        
        self._out.write(self._header)
        self._out.write(self._lines)
 
    def assoc_exon_with_gene(self) : 
                
        self._exon_gene_assoc = {}
        for k,v in self._attributes.items() : 
            name = k
            attrs = dict(x.split("=") for x in  v.split(";"))
            sequence_name = attrs['sequence_name'].replace("transcript:","")
            
            self._exon_gene_assoc[sequence_name] = name

    def parse_exons(self) :

        process_exons_utr(self._gff3, self._exon_gene_assoc, self._attributes, self._out, 'five_prime_UTR', 'WormBase')
        process_exons_utr(self._gff3, self._exon_gene_assoc, self._attributes, self._out, 'three_prime_UTR', 'WormBase')
        process_exons_utr(self._gff3, self._exon_gene_assoc, self._attributes, self._out, 'exon', 'WormBase')
        process_exons_utr(self._gff3, self._exon_gene_assoc, self._attributes, self._out, 'intron', 'WormBase')
        process_exons_utr(self._gff3, self._exon_gene_assoc, self._attributes, self._out, 'CDS', 'WormBase')
        self._out.close()
        
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
                    
                    if "Transposon:" in ID : 
                        Family = attrs['Family']
                        out.write(f"{chrom}\t{start}\t{end}\t{source}\t{feature}\t{strand}\t{name}\t{seq_name}\t{Family}\t{biotype}\t{Class}\n") 
                    else : 
                        out.write(f"{chrom}\t{start}\t{end}\t{source}\t{feature}\t{strand}\t{name}\t{seq_name}\t{locus}\t{biotype}\t{Class}\n") 

        f.close()
        out.close()

def run(gff, pathways, outfile) :

    x = process_gff3(gff, pathways, outfile)

    print(f">{gff}\n>{pathways}\n>{outfile}")

    print("Slim down gff3 file")
    x.slim()

    print("Parse pathways")
    x.parse_pathways()

    print("Parse genes")
    x.parse_genes()

    print("Associate exons with gene")
    x.assoc_exon_with_gene()

    print("Parse exons, UTR")
    x.parse_exons()

    print("Convert to bed")
    x.gff3_to_bed()
    
def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument("-gff3", type = str, required = True, help="gff3 file")
    parser.add_argument("-o", type = str, required = True, help="Output file")
    parser.add_argument("-classes", type = str, required = False, help="pathways file")
    
    return parser.parse_args()

def main() : 
    
    args = get_args()
    gff = args.gff3
    outfile = args.o
    pathways = args.classes
    
    run(gff, pathways, outfile)

if __name__ == "__main__" : 

    main()

##########
