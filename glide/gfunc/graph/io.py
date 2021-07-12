import numpy as np
import time
import sys
import mygene
import pkg_resources
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_ncbi_associations

def parse_go_label_file(fname):
    """Parses a GO label file.

    Outputs two dicts:
    - One associating GO labels to proteins
    - One associating proteins to GO labels
    """
    with open(fname, "r") as f:
        rows = f.readlines()
        go_to_proteins = {}
        proteins_to_go = {}
        for row in rows:
            words = row.split()
            golabel, proteins = words[0], words[1:]
            go_to_proteins[golabel] = proteins
            for protein in proteins:
                if protein in proteins_to_go:
                    proteins_to_go[protein].append(golabel)
                else:
                    proteins_to_go[protein] = [golabel]

        return go_to_proteins, proteins_to_go

def parse_graph_file(fname, weighted = True, normalize = True):
    """Parses a graph represented as an adjacency list. Works on either
    directed or undirected graphs.

    Outputs a triple:
    - An edgelist for the weighted graph G
    - A list mapping node indices to names
    - A dictionary mapping node names to node indices
    """
    with open(fname, "r") as f:
        edgelist, node_map = [], {}

        counter = 0
        edges = f.readlines()
        for edge in edges:
            edge = edge.rstrip().lstrip()
            if edge == "":
                continue
            edge = edge.split()
            if weighted:
                u, v, weight = edge[0], edge[1], float(edge[2])
            else:
                try:
                    u, v         = edge[0], edge[1]
                except:
                    print(f"{edge}")
                weight       = 1
            for x in [u, v]:
                if x not in node_map:
                    node_map[x] = counter
                    counter += 1

            edgelist.append((node_map[u], node_map[v], weight))

        n, m = len(node_map), len(edgelist)
        node_list = np.empty(n, dtype=object)
        for name, index in node_map.items():
            node_list[index] = name
        
        if normalize:
            mx = edgelist[0][2]
            for e in edgelist:
                if mx < e[2]:
                    mx = e[2]
            if mx > 1:
                for i in range(len(edgelist)):
                    p, q, w = edgelist[i]
                    edgelist[i] = (p, q, w / mx)
        return edgelist, node_list, node_map

def write_graph_to_file(filename, edgelist, node_map=None):
    """Given an output filename and a list of edges, writes the contents
    of the list to the file.
    """
    list_str = ""
    for edge in lst:
        p, q, wt = edge
        if node_map is not None:
            p = node_map[p]
            q = node_map[q]
            
        e_str = str(p) + " " + str(q) + " " + str(wt) + "\n"
        list_str += e_str

    list_str = list_str.lstrip().rstrip()
    
    with open(filename, "w") as fp:
        fp.write(list_str)

class GoTool:
    def __init__(self):
        gdagfile = pkg_resources.resource_filename('glide', 'data/go-basic.obo.dat')
        # self.godag = get_godag('go-basic.obo', optional_attrs='relationship')
        self.godag = get_godag(gdagfile, optional_attrs='relationship')
    def get_labels(self, filters = None):
        if filters == None:
            return list(self.godag.keys())
        go_terms   = []
        for k in self.godag.keys():
            k_val = self.godag[k]
            if "max_level" in filters:
                if k_val.level > filters["max_level"]:
                    continue
            if "min_level" in filters:
                if k_val.level < filters["min_level"]:
                    continue
            if "namespace" in filters:
                if k_val.namespace != filters["namespace"]:
                    continue
            go_terms.append(k)
        return go_terms
    def get_parents(self, go_id, filters = None):    
        gosubdag = GoSubDag(go_id, self.godag, relationships=True, prt=False)
        nts      = gosubdag.get_vals("id")
        parents  = []
        for n in nts:
            nObj = gosubdag.go2nt[n]
            ns   = nObj.NS
            if ns == "BP":
                ns = "biological_process"
            elif ns == "MF":
                ns = "molecular_function"
            else:
                ns = "cellular_component"
            level= nObj.level
            
            if (("max_level" in filters and filters["max_level"] < level) or
                ("min_level" in filters and filters["min_level"] > level) or
                ("namespace" in filters and filters["namespace"] != ns)):
                continue
            parents.append(nObj.id)
        return parents
        
def get_go_labels_and_parents(filter_protein, filter_label, filter_parent, entrez_labels, anno_map = lambda x : x):
    g2gofile = pkg_resources.resource_filename('glide', 'data/gene2go.dat')
    objanno = Gene2GoReader(g2gofile, taxids=[9606])
    go2geneids_human = objanno.get_id2gos(namespace=filter_protein["namespace"], 
                                          go2geneids=True)
    mg = mygene.MyGeneInfo()
    gt          = GoTool()
    labels      = gt.get_labels(filter_label)
    print("Read gene2go.dat File")
    labels_dict = {}
    f_labels    = []
    for key in labels:
        if key not in go2geneids_human:
            continue
        assoc_genes   = go2geneids_human[key]
        f_assoc_genes = [] 
        for a in assoc_genes:
            if a in entrez_labels:
                f_assoc_genes.append(anno_map(a))
        if len(f_assoc_genes) > filter_protein["lower_bound"]:
            labels_dict[key] = f_assoc_genes
            f_labels.append(key)
    print(f"Number of Labels: {len(f_labels)}")
    parent_dict = {}
    for l in f_labels:
       parent_dict[l] = gt.get_parents(l, filter_parent)
    return f_labels, labels_dict, parent_dict
