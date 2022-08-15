from sys import argv
from collections import namedtuple
from collections import defaultdict
import os
#features in gff
#gnl|Prokka|IIJKPODI_1	prokka	gene	55	558	.	-	.	ID=IIJKPODI_00001_gene;locus_tag=IIJKPODI_00001
#gnl|Prokka|IIJKPODI_1	Prodigal:002006	CDS	55	558	.	-	0	ID=IIJKPODI_00001;Parent=IIJKPODI_00001_gene;inference=ab initio prediction:Prodigal:002006;locus_tag=IIJKPODI_00001;product=hypothetical protein;protein_id=gnl|Prokka|IIJKPODI_00001
#faa:IIJKPODI_00001

def create_dir(working_dir):
    try:
        outpath = working_dir+"rm_system_annotation/" #output directory
        os.mkdir(outpath)
        print ("Path is created:"+outpath)
    except Exception:
        print ("Path exists:"+outpath)
    return outpath
def uncolored(d_color):
    for i in d_color:
        if d_color[i] == 0:
            return i
    return False
def color_func(lt, color, d_pair, d_color):
    if d_color[lt]==0:
        d_color[lt]=color
        for lt_nei in d_pair[lt]:
            color_func(lt_nei, color, d_pair, d_color)




def rm_annotate(rm_file, features_file, outpath):
    rm = open(rm_file, "r")
    features = open(features_file, "r")
    name = rm_file.split("/")[-1]
    out = open(outpath+name+".systems", "w")
    out.write(f'genome_faa\tprotein_id\trm_type#protein_type\tchrom\tstart\tstop\tstrand\n')

    cds = namedtuple('cds',['chrom', 'start', 'stop','strand','desc','locus_tag', 'prot'])
    d_cds = defaultdict(cds)
    d_lt_pos = defaultdict(int)
    d_lt = {}

    i=0
    for line in features:
        line = line.strip().split("\t")
        if len(line)>3:
            if line[2]=="CDS":
                chrom = line[0]
                start = line[3]
                stop = line[4]
                strand = line[6]
                desc = line[8]
                desc_lst = line[8].split(";")
                for item in desc_lst:
                    if "ID=" in item:
                        prot_id = item.split("=")[1]
                        if "cds-" in prot_id:
                            prot_id = prot_id[4:]
                    if "locus_tag=" in item:
                        locus_tag = item.split("=")[1]
                if locus_tag in d_lt:
                    d_lt[locus_tag].append(prot_id)
                else:
                    d_lt[locus_tag]=[]
                    d_lt[locus_tag].append(prot_id)
                d_cds[i] = cds(chrom,start,stop,strand,desc,locus_tag, prot_id)
                d_lt_pos[locus_tag] = i
                i+=1
        
    wps = []
    d_rm = {}

    for line in rm:
        line = line.strip()
        sline = line.split()
        wp = sline[1]
        wps.append(wp)
        d_rm[wp] = line 
    d_pair = defaultdict(set)
    
    for lt in d_lt:
        for prot in d_lt[lt]: 
            if prot in wps:
                k = d_lt_pos[lt]
                for z in range(k-4,k+5):
                    try:
                        if d_cds[z].prot in wps:
                            d_pair[lt].add(d_cds[z].locus_tag)
                            d_pair[d_cds[z].locus_tag].add(lt)
                    except:
                        pass
    
    d_color = {key:0 for key in d_pair}
    unc = uncolored(d_color)
    color = 1
    
    while unc:
        color_func(unc, color, d_pair, d_color)
        color+=1
        unc = uncolored(d_color)
    
    for key in d_color:
        for prot in d_lt[key]:
            if prot in d_rm:
                new_line = d_rm[prot].replace("_protein.faa.0.01.blast", "")
                rm_id = d_cds[d_lt_pos[key]].chrom+"-"+str(d_color[key])
                out.write(f'{new_line}\t{rm_id}\t{key}\t{d_cds[d_lt_pos[key]].chrom}\t{d_cds[d_lt_pos[key]].start}\t{d_cds[d_lt_pos[key]].stop}\t{d_cds[d_lt_pos[key]].strand}\n')
        
def main(rm, features, working_dir):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    rm_annotate(rm, features, outpath)
    print ("rm predictions are ready")

if __name__ == "__main__":
    rm = argv[1]
    features = argv[2]
    working_dir = argv[3] #output directory
    main(rm, features, working_dir) 

