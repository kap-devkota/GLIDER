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

verbose              = False

def log(msg):
    if verbose:
        print(msg)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--network", help="The network file.")
    parser.add_argument("--output", help  = "The output GLIDE list")
    parser.add_argument("--glide_loc", default = 'cw_normalized', choices=['l3', 'cw' , 'cw_normalized'], help = "Type of Local Score")
    parser.add_argument("-v", action="store_true", help="Verbose mode")
    parser.add_argument("--glide_params", default="1:0.1:1000:1:N", help = "GLIDE Params: format {gamma:alpha:beta:delta}. If UDSED then N, else (DSED) 'Y'")
    args = parser.parse_args()

    global verbose
    verbose = args.v
    network  = args.network

    log("Parsing network")
    edgelist, node_list, node_map = graph_io.parse_graph_file(network)
    A            = operations.densify(edgelist)
    n, _         = A.shape

    lam, alp, beta, delta, normalized = args.glide_params.split(":")
    lam = float(lam)
    alp = float(alp)
    beta = float(beta)
    delta = float(delta)
    glide_normalized = True if normalized == "Y" else False 
    log("Computing DSD embedding")

    X = computations.compute_embedding(edgelist, 
                                       lm = lam, 
                                       is_normalized = glide_normalized)
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
        


