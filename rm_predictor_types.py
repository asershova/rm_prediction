from sys import argv
from collections import defaultdict
from collections import namedtuple
import os
from typing import Dict, NamedTuple

def create_dir(working_dir):
    try:
        outpath = working_dir+"predicted_rm_proteins/" #output directory
        os.mkdir(outpath)
        print ("Path is created")
    except Exception:
        print (outpath+" exists")
    return outpath

class BlastData(NamedTuple):
    genome: str
    p_type: str
    query: str
    rm_prot: str
    bit_score: float
    aln: float




def make_prediction2(cdd_dict: Dict[str, float], blast_dict: Dict[str, BlastData], thresholds: Dict[str, NamedTuple]):
    
    #thrs = namedtuple('thrs',['dom_score','bit'])
    #thresholds = defaultdict(thrs)
    #thresholds = {'Orphan_M#M':thrs(0.75,50),'Type_IIG#RM':thrs(0.5,40), "Type_III#M":thrs(0.5,165),\
    #        "Type_III#R":thrs(0.5,50), "Type_II#M":thrs(0.75,50), "Type_IIM#R":thrs(0.6,30),"Type_II#R":thrs(0.75,45),\
    #        "Type_II#V":thrs(0,120), "Type_I#M":thrs(0.5,175), "Type_I#R":thrs(0.5,40), "Type_I#S":thrs(0.5,75),\
    #        "Type_IV#R":thrs(0.4,95), "-#C":thrs(1,0),\
    #        "Type_II#H":thrs(1,0),"Type_II#S":thrs(1,0),\
    #            "Type_I#RM":thrs(1,0), "Homing#R":thrs(1,0) }
    clean_cdd = {}
    
    for ptype in cdd_dict:
        if ptype in blast_dict:
            if ptype in thresholds:
                if cdd_dict[ptype] > thresholds[ptype].dom_score:
                    if blast_dict[ptype].bit_score >= thresholds[ptype].bit:
                        clean_cdd[ptype]=blast_dict[ptype].bit_score
    if clean_cdd == {}:
        return ("0")
    else:
        best = 0
        for ptype in clean_cdd:
            if clean_cdd[ptype]>best:
                best = clean_cdd[ptype]
                best_ptype = ptype
    return best_ptype



#def make_rebase_dict(): #rebase file should be at the current directory
#    rebase_gs = open("rebase_007.only_gold.tsv", "r")
#    rebase_anno = namedtuple('rebase_anno', ['prot_type','site'])
#    rebase_dict = defaultdict(rebase_anno)
#    for line in rebase_gs:
#        line = line.strip().split("\t")
#        rebase_dict[line[0]] = rebase_anno(line[2].replace(" ","_")+"#"+ line[3], line[4])
#    return rebase_dict
def get_thresholds(thr_file):
    f = open(thr_file, "r")
    thrs = namedtuple('thrs',['dom_score','bit'])
    thresholds = defaultdict(thrs)
    for line in f:
        line = line.strip().split()
        thresholds[line[0]] = thrs(float(line[1]), float(line[2]))
        
    #thresholds = {'Orphan_M#M':thrs(0.75,50),'Type_IIG#RM':thrs(0.5,40), "Type_III#M":thrs(0.5,165),\
    #        "Type_III#R":thrs(0.5,50), "Type_II#M":thrs(0.75,50), "Type_IIM#R":thrs(0.6,30),"Type_II#R":thrs(0.75,45),\
    #        "Type_II#V":thrs(0,120), "Type_I#M":thrs(0.5,175), "Type_I#R":thrs(0.5,40), "Type_I#S":thrs(0.5,75),\
    #        "Type_IV#R":thrs(0.4,95), "-#C":thrs(1,0),\
    #        "Type_II#H":thrs(1,0),"Type_II#S":thrs(1,0),\
    #            "Type_I#RM":thrs(1,0), "Homing#R":thrs(1,0) }
    return thresholds
def write_results(cdd_file, blast_file,m,mode,outpath, thresholds_file):
    cdd = open(cdd_file, "r")
    blast = open(blast_file, "r")
    out_name = outpath+blast_file.split("/")[-1]+"."+m+"."+mode+".predicted_rm_types"
    out = open(out_name, "w")
    thresholds = get_thresholds(thresholds_file)
    cdd_dict = defaultdict(dict)
    #rebase_dict = make_rebase_dict()
    blast_dict = defaultdict(dict)

    for line in cdd:
        line = line.strip().split("\t")
        type_name = line[-1].split(".")[0]
        cdd_dict[line[0]][type_name] = float(line[-2])
    
    for line in blast:
        if line[0]!='#':
            sline = line.strip().split("\t")
            #p_type = rebase_dict[line[2]].prot_type
            query = sline[2]
            p_type = sline[1]
            blast_dict[query][p_type] = BlastData(
                genome = sline[0],
                p_type = sline[1],
                query = sline[2],
                rm_prot = sline[3],
                bit_score = float(sline[13]),
                aln = float(sline[4])
            )

    for query in blast_dict:
        if query not in cdd_dict:
            for p_type in blast_dict[query]:
                bd = blast_dict[query][p_type]
                if bd.aln > 95:
                    out.write(f'{bd.genome}\t{bd.query}\t{p_type}\t{bd.aln}\t{bd.bit_score}\tno\t-1\tBlast\n')
        else:
            
            best_type = make_prediction2(cdd_dict[query], blast_dict[query], thresholds)
             
            if best_type !="0":
                bd = blast_dict[query][best_type]
                dom_score = cdd_dict[query][best_type]
                print(bd.genome+" predicted")
                out.write(f'{bd.genome}\t{bd.query}\t{best_type}\t{bd.aln}\t{bd.bit_score}\t{best_type}\t{dom_score}\tDomains\n')
    
           
def main(working_dir, cdd_file, blast_file, m, mode, thresholds_file):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    write_results(cdd_file, blast_file, m, mode, outpath, thresholds_file)
    print ("rm systems are predicted")

if __name__ == "__main__":
    working_dir = argv[1]
    cdd_file = argv[2]
    blast_file = argv[3]
    mode = argv[4]
    m = argv[5]
    thresholds_file = argv[6]
#score_threshold = float(argv[3])
#bit-score_threshold = float(argv[4])
    main(working_dir, cdd_file, blast_file, m, mode, thresholds_file)
