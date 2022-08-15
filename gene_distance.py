from sys import argv
from collections import namedtuple
from collections import defaultdict
rm = open(argv[1], "r")
features = open(argv[2], "r")
outdir = argv[3] #output directory
name = argv[1].split("/")[-1]
out = open(outdir+"/"+name+".systems", "w")
cds = namedtuple('cds',['genome','chrom', 'start', 'stop','strand','prot','locus_tag'])
d_cds = defaultdict(cds)
d_lt_pos = defaultdict(int)
d_lt = defaultdict(str)

i=0
for line in features:
    line = line.strip().split("\t")
    if line[0]=="CDS":
        genome = line[2]
        chrom = line[6]
        start = line[7]
        stop = line[8]
        strand = line[9]
        prot = line[10]
        locus_tag = line[16]
        d_cds[i] = cds(genome,chrom,start,stop,strand,prot,locus_tag)
        d_lt_pos[locus_tag] = i
        d_lt[locus_tag] = prot

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
print (d_rm)
for lt in d_lt:
    if d_lt[lt] in wps:
        k = d_lt_pos[lt]
        for z in range(k-4,k+5):
            try:
                if d_cds[z].prot in wps:
                    d_pair[lt].add(d_cds[z].locus_tag)
                    d_pair[d_cds[z].locus_tag].add(lt)
            except:
                pass
print (d_pair)
d_color = {key:0 for key in d_pair}
def uncolored(d_color):
    for i in d_color:
        if d_color[i] == 0:
            return i
    return False
unc = uncolored(d_color)
color = 1
def color_func(lt, color, d_pair, d_color):
    if d_color[lt]==0:
        d_color[lt]=color
        for lt_nei in d_pair[lt]:
            color_func(lt_nei, color, d_pair, d_color)

    
while unc:
    color_func(unc, color, d_pair, d_color)
    color+=1
    unc = uncolored(d_color)

for key in d_color:
    assembly_id = d_rm[d_lt[key]].split("\t")[0]
    rm_id = assembly_id+"-"+str(d_color[key])
    print (d_rm[d_lt[key]])
    out.write(f'{d_rm[d_lt[key]]}\t{rm_id}\t{key}\t{d_cds[d_lt_pos[key]].chrom}\t{d_cds[d_lt_pos[key]].start}\t{d_cds[d_lt_pos[key]].stop}\
            \t{d_cds[d_lt_pos[key]].strand}\n')
