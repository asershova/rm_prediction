from sys import argv
from collections import namedtuple
from collections import defaultdict
import os
from rm_predictor_types import BlastData
from typing import Dict, NamedTuple, List, Tuple, Set
#features in gff
#gnl|Prokka|IIJKPODI_1	prokka	gene	55	558	.	-	.	ID=IIJKPODI_00001_gene;locus_tag=IIJKPODI_00001
#gnl|Prokka|IIJKPODI_1	Prodigal:002006	CDS	55	558	.	-	0	ID=IIJKPODI_00001;Parent=IIJKPODI_00001_gene;inference=ab initio prediction:Prodigal:002006;locus_tag=IIJKPODI_00001;product=hypothetical protein;protein_id=gnl|Prokka|IIJKPODI_00001
#faa:IIJKPODI_00001

def create_dir(working_dir):
    try:
        outpath = working_dir+"rm_system_with_RE/" #output directory
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

def read_rm_dict(rm_file: str ) -> Dict[str, Dict[str,   List[str]]]:
                                     # type      protein genome      
    #rm
#GCA_000003955.1_ASM395v1_protein.faa.0.01.blast EEL90457.1      Type_II#M       25.301  63.5    Type_II#M       1.0     Domains
    rm = open(rm_file, "r")
    res = defaultdict(dict)
    for line in rm:
        line = line.strip()
        sline = line.split()
        genome = sline[0]
        wp = sline[1]
        ptype = sline[2]
        res[ptype][wp]=sline
    return res

class CDS(NamedTuple):
    chrom: str
    start: str
    stop: str 
    strand: str
    desc: str
    locus_tag: str
    prot: str

def cds_read_dict(cds_file: str) -> Tuple[Dict[int, CDS], Dict[str, int], Dict[str, str]]:
    features = open(cds_file, "r")
    d_cds = defaultdict(CDS)
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
                prot_id = ""
                locus_tag = ""
                for item in desc_lst:
                    if "ID=" in item:
                        prot_id = item.split("=")[1]
                        if "cds-" in prot_id:
                            prot_id = prot_id[4:]
                    if "locus_tag=" in item:
                        locus_tag = item.split("=")[1]
                if prot_id == "":
                    print(f'Warning: no prot_id for line: {line}')
                if locus_tag == "":
                    print(f'Warning: no locus_tag for line: {line}')
                if locus_tag in d_lt:
                    print (f'Warning:{locus_tag} is duplicated')
                else:
                    d_lt[locus_tag]=prot_id
                d_cds[i] = CDS(chrom,start,stop,strand,desc,locus_tag, prot_id)
                d_lt_pos[locus_tag] = i
                i+=1
    return (d_cds, d_lt_pos, d_lt)

def read_blast_dict(blast_file): #-> str, Dict[str, BlastData ]]
                        # type     protein   genome
    #GCA_000003955.1_ASM395v1_protein.faa.0.01.blast Type_II#R       EEL90114.1      HpyAII  30.108  93      54      2       61      153     237     318      0.004   36.6
    blast = open(blast_file, "r")
    blast_dict = defaultdict(dict)
    for line in blast:
        if line[0]!='#':
            sline = line.strip().split("\t")
            #p_type = rebase_dict[line[2]].prot_type
            query = sline[2]
            p_type = sline[1]
            blast_dict[p_type][query] = BlastData(
                genome = sline[0],
                p_type = sline[1],
                query = sline[2],
                rm_prot = sline[3],
                bit_score = float(sline[13]),
                aln = float(sline[4])
            )
    return blast_dict

def get_existing_rm_proteins(rm_dict: Dict[str, Dict[str, str]]) -> Set[str]:
    res = set()
    for ptype, pnames in rm_dict.items():
        for pname in pnames:
            res.add(pname)
    return res

def re_search(rm_file, features_file, blast_file):
    rm_dict = read_rm_dict(rm_file)
    existing_rm = get_existing_rm_proteins(rm_dict)
    blast_dict = read_blast_dict(blast_file)
    #name = rm_file.split("/")[-1]
    #out = open(outpath+name+".added_R", "w")
    #out.write(f'genome_faa\tprotein_id\trm_type#protein_type\tchrom\tstart\tstop\tstrand\n')
    d_cds, d_lt_pos, d_lt = cds_read_dict(features_file)
    mt_types = { 
            'Orphan_M#M': 'Type_II#R',
            'Type_III#M': 'Type_III#R',
            'Type_II#M': 'Type_II#R',
            'Type_I#M': 'Type_I#R'
            }
    for lt, prot in d_lt.items():
        for m_type, r_type in mt_types.items():
            if prot in rm_dict[m_type]:
                k = d_lt_pos[lt]
                for z in range(k-4,k+5):
                    
                    try:                       
                        zprot = d_cds[z].prot
                        if zprot in blast_dict[r_type]:
                            if zprot not in existing_rm:
                                rm_dict[r_type][zprot] = [blast_dict[r_type][zprot].genome, zprot, r_type, str(blast_dict[r_type][zprot].aln), str(blast_dict[r_type][zprot].bit_score), "no", "-1", "Neighbour"]
                    except:
                        pass
    return rm_dict
    
def write_dict(rm_file, rm_dict, outpath):
    name = rm_file.split("/")[-1]
    out = open(outpath+name+".added_R", "w")
    for rm_type in rm_dict:
        for prot in rm_dict[rm_type]:
            line = "\t".join(rm_dict[rm_type][prot])
            out.write(f'{line}\n')
        
def main(rm, features, blast_file, working_dir):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    rm_new = re_search(rm, features,blast_file)
    write_dict(rm, rm_new, outpath)
    print ("new REs are added")

if __name__ == "__main__":
    rm = argv[1]
    features = argv[2] #gff file
    blast_file = argv[3]
    working_dir = argv[4] #output directory
    main(rm, features, blast_file, working_dir) 

