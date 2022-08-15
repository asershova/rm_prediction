from sys import argv
from collections import defaultdict, namedtuple
import csv
import os

#‘qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
#0   1   2   3   4   5   6   7   8   9   10  11
def create_dir(working_dir):
    directory = working_dir+"blast_best_hits/"
    try:
       os.mkdir(directory)
    except:
       pass
    return directory
def best_blast_hits(blast_dir, outdir):
    entries = os.scandir(blast_dir)
    rebase_dict = make_rebase_dict()
    for entry in entries:
        blast_dict = defaultdict(list)
        blast = open(f'{blast_dir}{entry.name}', "r")
        out1 = open(f'{outdir}{entry.name}.best_hits', "w")
        for line in blast:
            sline = line.strip().split()
            query = sline[0]
            subject = sline[1]
            aln_len = sline[3]
            bitscore = sline[11]
            rm_type = rebase_dict[subject].prot_type
            key = (query, rm_type)
            if key in blast_dict:
                if float(blast_dict[key][11]) < float(bitscore):
                    blast_dict[key] = sline
                else:
                    continue
            else:
                blast_dict[key] = sline
        
        for i in blast_dict:
            blast_string = '\t'.join(blast_dict[i])
            genome_id = entry.name
            out1.write(f'{genome_id}\t{i[1]}\t{blast_string}\n')
def make_rebase_dict(): #rebase file should be at the current directory
    #rebase_gs = open("rebase_007.only_gold.tsv", "r")
    rebase_gs = open("rebase_007/rebase_007.tsv", "r")
    rebase_anno = namedtuple('rebase_anno', ['prot_type','site'])
    rebase_dict = defaultdict(rebase_anno)
    for line in rebase_gs:
        line = line.strip().split("\t")
        rebase_dict[line[0]] = rebase_anno(line[2].replace(" ","_")+"#"+ line[3], line[4])
    return rebase_dict
def main(blast_dir, working_dir):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    best_blast_hits(blast_dir, outpath)
    print ("best hits results are ready")

if __name__ == "__main__":
    blast_dir = argv[1]
    working_dir = argv[2]
    main(blast_dir, working_dir)
