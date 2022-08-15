#from sys import argv
import argparse
from collections import defaultdict
from typing import NamedTuple, List, Set
from std_asn_out_objects_set import Domain, Query, line_to_dom, read_asn_out 
from comparison_lib import * 
import os
from sys import argv
import pdb

def create_dir(working_dir):
    try:
        outpath = working_dir+"domain_comparison/" #output directory
        os.mkdir(outpath)
        print ("Path is created:"+outpath)
    except Exception:
        print ("Path exists:"+outpath)
    return outpath
def compare_domains(doms_dir, mode, modified, outdir):
    
    for i in os.scandir(doms_dir):
        name = i.name.split("/")[-1]
        if name.split(".")[-1] == "out":
            out = outdir+name+"."+str(modified)+"."+mode
            dom_dict = read_asn_out(doms_dir+i.name)
            gold_dicts = {}
            gold = "gold_list"
            #gold_path = "/mnt/md0/anna/rm_systems/scripts/"
            gold_file = open(gold, "r")
            for line in gold_file:
                line = line.strip().split()
                f = line[0]
                if int(modified) >0:
                    gold_dict = gold_dict_edit(read_asn_out(f),f)
                else:
                    gold_dict = read_asn_out(f)
                gold_dicts[line[0]]=gold_dict
            for rm_type, gold_dict in gold_dicts.items():
                d_match = make_d_match(dom_dict, gold_dict, mode)
                d_best_score = get_best_score(dom_dict, gold_dict, d_match, mode)
                print_best_scores(d_best_score, rm_type, out) #from gold_dicts of different rm-types
def main(working_dir, asn_dir, mode, modified):
    outpath = create_dir(working_dir)
    print ("See results in "+outpath)
    compare_domains(asn_dir,mode,modified,outpath)
    print ("Comparison results are ready")

if __name__ == "__main__":
    asn_dir = argv[1]
    working_dir = argv[2]
    mode = argv[3] #family or superfamily
    modified = argv[4] #0 or 1 if change domains
    main(working_dir, asn_dir, mode, modified)
