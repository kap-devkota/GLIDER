"""
Before running this file, first go to the denoise folder and do:
         pip install . 
in order to install the denoise package into the system.
"""
import argparse
import random
import numpy as np
import json
import denoise.graph.io as graph_io
from denoise.algorithms.dsd import denoise 
from denoise.algorithms.dsd import computations
from denoise.graph import operations
from denoise import scoring
from denoise import predict
import ctypes
import scipy.spatial.distance as spatial
import pandas as pd
import pickle

def parse_args():
    """
    Function for parsing the command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--network", help="A Tab-delimited network file")
    parser.add_argument("--output", help  = "The output GLIDE list")
    parser.add_argument("--glide_loc", default = 'cw_normalized', choices=['l3', 'cw' , 'cw_normalized'], help = "Type of Local Score")
    parser.add_argument("-v", action="store_true", help="Verbose mode")
    
    """
    DSD arguments
    """
    parser.add_argument("--return-dsd-emb", 
                        action = "store_true", 
                        default = False, 
                        help = "If set to True, only returns the DSD embedding, else returns the GLIDE matrix")
    
    parser.add_argument("--return-dsd-dist", 
                        action = "store_true", 
                        default = False, 
                        help = "If set to True, bypasses the --return-dsd-emb command to return the pairwise distance matrix from the dsd embedding")
    parser.add_argument("--dsd-dist-norm", 
                        choices = ["l1", "l2"], 
                        default = "l2", 
                        help = "Only used in conjunction with the --return-dsd-dist argument. Decides whether to use the `l1` or `l2` norm while computing distance")
    parser.add_argument("--normalized", 
                        action = "store_true", 
                        default = False, 
                        help = "If set to false, returns the classic cDSD, else returns normalized cDSD embedding.")
    
    """
    GLIDE arguments
    """
    parser.add_argument("--gamma", default = 1.0, type = float, help = "GLIDE gamma parameter")
    parser.add_argument("--alpha", default = 0.1, type = float, help = "GLIDE alpha parameter")
    parser.add_argument("--beta", default = 1000, type = float, help = "GLIDE beta parameter")
    parser.add_argument("--delta", default = 1, type = "float", help = "GLIDE delta parameter")
    
    return parser.parse_args()


def main():
    args = parser.parse_args()
    print(args)
    return
    
    def log(msg):
        """
        If verbose, print
        """
        if args.v:
            print(msg)

    network  = args.network

    ## Parsing the network
    log("Parsing network")
    dfnet    = pd.read_csv(network, sep = "\t", header = None)
    nodelist = list(set(dfnet[0]).union(dfnet[1]))
    nodemap  = {k:i for i, k in enumerate(nodelist)}
    
    
    edgelist, node_list, node_map = graph_io.parse_graph_file(network)
    A            = operations.densify(edgelist)
    n, _         = A.shape

    
    ## Getting the GLIDE Parameters
    gamma, alpha, beta, delta, gn = args.gamma, args.alpha, args.beta, args.delta, args.normalized
    
    log("Computing DSD embedding")
    X = computations.compute_embedding(edgelist, 
                                       lm = lam, 
                                       is_normalized = gn)
    params = {"alpha" : alp, "beta" : beta, "delta" : delta, "loc" : args.glide_loc}
    log("Predicting GLIDE links")
    ranked_edgelist = denoise.glide_predict_links(edgelist, 
                                                  X, 
                                                  params = params)
    log("Writing GLIDE list")

    r_nodemap = {}
    for key in node_map:
        r_nodemap[node_map[key]] = key

    with open(f"{args.output}", "w") as of:
        str_ = ""
        for r in ranked_edgelist:
            p, q, w  = r
            p        = r_nodemap[p]
            q        = r_nodemap[q]
            str_    += f"{p} {q} {w}\n"
        of.write(str_)

if __name__ == "__main__":
    main()
        


