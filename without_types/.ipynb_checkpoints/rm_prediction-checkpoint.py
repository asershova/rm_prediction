from sys import argv
from collections import defaultdict
from collections import namedtuple

cdd = open(argv[1], "r")
blast = open(argv[2], 'r')
#score_threshold = float(argv[3])
#bit-score_threshold = float(argv[4])

out_name = argv[2]+".predicted_rm"
out = open(out_name , "w")

def make_prediction(prot_type, dom_s, bit_s):
    thrs = namedtuple('thrs',['dom_score','bit'])
    thresholds = defaultdict(thrs)
    thresholds = {'Orphan_M#M':thrs(0,50),'Type_IIG#RM':thrs(0,40), "Type_III#M":thrs(0,165),\
            "Type_III#R":thrs(0,50), "Type_II#M":thrs(0,50), "Type_IIM#R":thrs(0,30),"Type_II#R":thrs(0,45),\
            "Type_II#V":thrs(0,120), "Type_I#M":thrs(0,175), "Type_I#R":thrs(0,40), "Type_I#S":thrs(0,75),\
            "Type_IV#R":thrs(0,95), "-#C":thrs(1,0), "Homing#R":thrs(1,0), "Orphan_M#M":thrs(0,50)}
    if dom_s > thresholds[prot_type].dom_score:
        if bit_s >= thresholds[prot_type].bit:
            return True
    return False
def make_rebase_dict(): #rebase file should be at the current directory
    rebase_gs = open("rebase_007.only_gold.tsv", "r")
    rebase_anno = namedtuple('rebase_anno', ['prot_type','site'])
    rebase_dict = defaultdict(rebase_anno)
    for line in rebase_gs:
        line = line.strip().split("\t")
        rebase_dict[line[0]] = rebase_anno(line[2].replace(" ","_")+"#"+ line[3], line[4])
    return rebase_dict

cdd_dict = defaultdict(float)

for line in cdd:
    line = line.strip().split("\t")
    cdd_dict[line[0]] = float(line[-1])
rebase_dict = make_rebase_dict()
for line in blast:
    if line[0]!='#':
        line = line.strip().split("\t")
        genome = line[0]
        prot = line[1]
        bit_score = float(line[12])
        p_type = rebase_dict[line[2]].prot_type
        if make_prediction(p_type, cdd_dict[prot],bit_score):
            print ("yes")
            out.write(f'{genome}\t{prot}\t{p_type}\n')
        

