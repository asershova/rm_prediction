from sys import argv
in_file = open("proteins_faa/faa_list", "r")
besthits_file = open("blast_best_hits_anna/" + argv[1], "r")
out = open("seqret_anna/" + argv[1] + ".seqret", "w")
#create two dicts
rm_dict = {}
genome_dict = {}
#rm_prediction_file
#read line by line, strip and split lines
for line in besthits_file:
    sline = line.strip().split()
    genome = sline[0]
    protein_id = sline[1]
    rm_dict[protein_id] = genome
#faa_file
#read line by line, strip and split lines
for line in in_file:
    sline = line.strip().split()
    sfilename = sline[0].split("_")
    short_name = sfilename[0]+"_"+sfilename[1]
    genome_dict[short_name] = sline[0]
for key in rm_dict:
    out.write(genome_dict[rm_dict[key]] + ":" + key + "\n")
