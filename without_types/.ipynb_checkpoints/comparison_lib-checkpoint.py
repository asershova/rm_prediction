#from sys import argv
import argparse
from collections import defaultdict
from typing import NamedTuple, List, Set
from std_asn_out_objects_set import Domain, Query, line_to_dom, read_asn_out 
import csv

def gold_dict_edit(gold_dict):
    gold_dict["49#Query_1043"].best_domains[0].append(Domain(hit_type='Non-specific', PSSM_ID='236618', start=1, stop=348, e_value=3.87155e-142, bitscore=409.539, accession='COG4268', short_name='McrC', incomplete='N', superfamily_PSSM_ID='382678', true_superfamily='cl01737'))
    gold_dict["49#Query_1043"].best_domains[0].append(Domain(hit_type='Non-specific', PSSM_ID='370821', start=13, stop=304, e_value=4.8546e-20, bitscore=88.8876, accession='pfam10117', short_name='McrBC', incomplete= '-', superfamily_PSSM_ID='382678', true_superfamily='cl01737'))
    gold_dict["10#Query_214"].best_domains[0].append(Domain(hit_type='Non-specific', PSSM_ID='372641', start=28, stop=169, e_value=1.25851e-17, bitscore=79.6608, accession='pfam13589', short_name='HATPase_c_3', incomplete= '-', superfamily_PSSM_ID='381791', true_superfamily='cl00075'))
    gold_dict["57#Query_1216"].best_domains[0].append(Domain(hit_type='Non-specific', PSSM_ID='372641', start=33, stop=157, e_value=1.36997e-14, bitscore=70.8013, accession='pfam13589', short_name='HATPase_c_3', incomplete= '-', superfamily_PSSM_ID='381791', true_superfamily='cl00075'))
    return gold_dict

class Match():
    query_list:List[int]
    subj_list:List[int]
    non_zero_ratio: float
    def __init__(self, prot, gs, mode, bad_doms=set()):
        prot_list = get_family_list(prot, mode)
        gs_list = get_family_list(gs, mode)
        self.query_list = [0 for n in prot_list]
        self.subj_list = [0 for k in gs_list]
        for j in range(len(gs_list)):
            gs_list_fin = gs_list[j]-bad_doms
            for i in range(len(prot_list)):
                if len(prot_list[i].intersection(gs_list_fin)) > 0:
                    self.query_list[i] +=1
                    self.subj_list[j] +=1
        self.non_zero_ratio = self.non_zero()
    def non_zero (self) -> int:
        countq = 0
        counts = 0
        for i in self.query_list:
            if i !=0:
                countq+=1
        for j in self.subj_list:
            if j !=0:
                counts+=1
        return (countq+counts)/(len(self.query_list)+len(self.subj_list))
def get_family_list (query, mode):
    if mode == "superfamily":
        return query.superfamilies()
    if mode == "family":
        return query.families()
   
def make_d_match (dom_dict, gold_dict, mode, bad_doms = set()):
    d_match = defaultdict(dict)
    for prot in dom_dict:
        for rm in gold_dict:
            match = Match(dom_dict[prot],gold_dict[rm], mode, bad_doms)
            if match.non_zero_ratio > 0:
                d_match[prot][rm] = match
    return d_match

class Best_Score (NamedTuple):
    domains:str
    best_rm_name:str
    best_rm:str
    max_score:float


def get_best_score (dom_dict, gold_dict, d_match, mode):
    d_best_score = {}
    for prot in dom_dict.keys():
        domains = get_family_list(dom_dict[prot], mode)
        if prot in d_match:
            max_score = 0
            best_rm = "unk"
            best_rm_name = "unk"
            for rm in d_match[prot]:
                if d_match[prot][rm].non_zero_ratio > max_score:
                    max_score = d_match[prot][rm].non_zero_ratio
                    best_rm = rm
                    best_rm_name = gold_dict[best_rm].name
            d_best_score[(dom_dict[prot].name,prot)] = Best_Score(
                domains = ",".join(y for x in domains for y in x),
                best_rm_name = best_rm_name, 
                best_rm = best_rm, 
                max_score = max_score)    
        else:
            d_best_score[(dom_dict[prot].name,prot)] = Best_Score(
                domains = ",".join(y for x in domains for y in x),
                best_rm_name = "no", 
                best_rm = "no", 
                max_score = 0)
    return d_best_score

def print_best_score (d_best_score, out):
    w = csv.writer(open(out, "w"), delimiter = "\t")
    for key, val in d_best_score.items():
        w.writerow([key[0].split()[0], val.domains, val.best_rm_name, val.max_score])
        

        



