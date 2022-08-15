from sys import argv
from collections import defaultdict
#>Cluster 0
#0       3227aa, >WP_001070504.1... *
#1       3227aa, >WP_001070504.1... at 100.00%
#2       3227aa, >WP_001070504.1... at 100.00%
#3       3227aa, >WP_001070504.1... at 100.00%
#4       2250aa, >WP_115596477.1... at 99.96%
#>Cluster 1
#0       2210aa, >WP_151532485.1... *
clstr = argv[1]
out1 = open(argv[1]+".table", "w")
def parse_clusters(filename):
    f = open(filename, "r")
    seq_dict = defaultdict(str)
    for line in f:
        line = line.strip()
        if line[0] == ">":
            #print(str(wp)+","+str(rebase))
            clstr_id = line.split(" ")[1]
        else:
            line = line.split(" ")
            prot_id = line[1][1:-3]
            if prot_id in seq_dict and seq_dict[prot_id]!=clstr_id:
                print ("The same protein in different clusters!: "+prot_id, clstr_id, seq_dict[prot_id])
            seq_dict[prot_id] = clstr_id
    return seq_dict
seq_dict = parse_clusters(clstr)
for prot in seq_dict:
    out1.write(f'{prot}\t{seq_dict[prot]}\n')
print ("The result is in "+argv[1]+".table")
