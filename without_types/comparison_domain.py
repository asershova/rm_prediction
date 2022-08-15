#from sys import argv
import argparse
from collections import defaultdict
from typing import NamedTuple, List, Set
from std_asn_out_objects_set import Domain, Query, line_to_dom, read_asn_out 
from comparison_lib import * 
parser = argparse.ArgumentParser(description='Compare domain composition of proteins from two files')
parser.add_argument("file1", help="asn file with cdd domain results for query")
parser.add_argument("file2", help="asn file with cdd domain results for comparison file,\
        like rebase gold standard")
parser.add_argument("mode", help="comparison mode: can be family or superfamily")
parser.add_argument("--modified", help="add several domains in gold_standard data", action="store_true")
args = parser.parse_args()

doms = args.file1
gold = args.file2
mode = args.mode
if args.modified:
    out = args.file1+".modified."+mode
else:
    out = args.file1+"."+mode
dom_dict = read_asn_out(doms)

if args.modified:
    gold_dict = gold_dict_edit(read_asn_out(gold))
else:
    gold_dict = read_asn_out(gold)


d_match = make_d_match(dom_dict, gold_dict, mode)
d_best_score = get_best_score(dom_dict, gold_dict, d_match, mode)
print_best_score(d_best_score, out)



