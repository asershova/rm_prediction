from sys import argv
from collections import defaultdict
from typing import NamedTuple, List, Set, Dict

class Domain(NamedTuple):
    hit_type: str
    PSSM_ID: str
    start: int
    stop: int
    e_value: float
    bitscore: float
    accession: str
    short_name: str
    incomplete: str
    superfamily_PSSM_ID: str
    true_superfamily: str
class Query():
    query_id: str
    name: str
    domains: List[Domain]
    best_domains: List[Domain]
    
    def __init__(self, session: str, line: List[str]):  
        self.query_id = session + "#" + line[1]
        self.name = line[4]
        self.domains = []
        self.best_domains = []    

    def get_best_domains (self)-> List[Domain]:
        list_of_best = []
        for dom in self.domains:
            if dom.hit_type == "Specific":
                i = 0
                for dom_list in list_of_best:
                    for best_dom in dom_list:
                        #if not ((dom.stop < best_dom.start and dom.start != best_dom.stop)\
                        #       or (best_dom.stop < dom.start and dom.stop != best_dom.start)):
                        if min(dom_intersect(dom, best_dom)[1], dom_intersect(dom, best_dom)[2])  > 0.75:
                            dom_list.append(dom)
                            i = 1
                            break
                if i == 0:
                    new_list = []
                    new_list.append(dom)
                    list_of_best.append(new_list)
        if len(list_of_best) > 0:
            return list_of_best
        for dom in self.domains:
            if dom.hit_type == "Non-specific":
                i = 0
                if len(list_of_best)>0:
                    for dom_list in list_of_best:
                        for best_dom in dom_list:
                            #if not ((dom.stop < best_dom.start and dom.start != best_dom.stop)\
                            #       or (best_dom.stop < dom.start and dom.stop != best_dom.start)):
                            if min(dom_intersect(dom, best_dom)[1], dom_intersect(dom, best_dom)[2])  > 0.75:
                                i = 1
                                if dom.e_value < best_dom.e_value:
                                    dom_list.remove(best_dom)
                                    dom_list.append(dom)
                                                               
                    if i == 0:
                        new_list = []
                        new_list.append(dom)
                        list_of_best.append(new_list)
                else:
                    new_list = []
                    new_list.append(dom)
                    list_of_best.append(new_list)

        return list_of_best
       
                
    def superfamilies (self) -> List[Domain]:
        sp_dom_set = []
        for doms in self.best_domains:
            dom_set = set()
            for dom in doms:
                dom_set.add(dom.true_superfamily)
            sp_dom_set.append(dom_set)
        return sp_dom_set
    def families (self) -> Set[str]:
        sp_dom_set = []
        for doms in self.best_domains:
            set_dom = set()
            for dom in doms:
                if isinstance(dom, Domain):
                    set_dom.add(dom.short_name)
                else:
                    import pdb 
                    pdb.set_trace()
            sp_dom_set.append(set_dom)
        return sp_dom_set

def line_to_dom (line:List[str], family_dict)->Domain:
    #1       Query_1 Superfamily     390081  12      283     7.806e-137      387.154 cl23779 MethyltransfD12 -       -
    return Domain(hit_type = line[2], PSSM_ID = line[3],
            start = int(line[4]), stop = int(line[5]),
            e_value = float(line[6]), bitscore = float(line[7]),
            accession = line[8], short_name = line[9],
            incomplete = line[10], superfamily_PSSM_ID = line[11], true_superfamily = family_dict[line[8]])
family_links = open("family_superfamily_links", "r")
family_dict = defaultdict()
for line in family_links:
    line = line.strip().split()
    if len(line)>0:
        family_dict[line[0]] = line[2]

def read_asn_out (filename) -> Dict[str, Query]:
    d_id = defaultdict()
    f = open(filename, "r")
    session = 0
    for line in f:
        if line[0]!="#" and len(line) > 0:
            line = line.strip().split("\t")
            if line[0] == "SESSION":
                session = line[1]            
            if line[0] == "QUERY":
                current_query=Query(session, line)
                d_id[current_query.query_id] = current_query
            
            if line[0] == session:
                dom = line_to_dom(line, family_dict)
                current_query.domains.append(dom)
        
    for query in d_id:
        d_id[query].best_domains = d_id[query].get_best_domains()
        
    return d_id
def make_all_dom_dict (d_id):
    all_dom_dict = defaultdict(int)
    for query in d_id:
        for dom in d_id[query].domains:
            all_dom_dict[dom.short_name]+=1
    return all_dom_dict 
def dom_intersect(dom1, dom2):
    max_start = max(dom1.start, dom2.start)
    min_stop = min(dom1.stop, dom2.stop)
    if max_start < min_stop:
        interlen = min_stop - max_start
        d1_interlen = interlen/(dom1.stop - dom1.start)
        d2_interlen = interlen/(dom2.stop - dom2.start)
    else:
        return 0, 0, 0
    return interlen, d1_interlen, d2_interlen

