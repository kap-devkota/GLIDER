"""
Before running this file, first go to the denoise folder and do:
         pip install . 
in order to install the glidetool package into the system.
"""
import argparse
import random
import numpy as np
import json
import networkx as nx

import glidetools.algorithm.dsd as dsd
import glidetools.algorithm.glide as glide
import glidetools.construct_graph as cg

from scipy.spatial.distance import squareform, pdist
import pandas as pd
import pickle

def parse_args():
    """
    Function for parsing the command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--network", help="A Tab-delimited network file")
    parser.add_argument("--output", help  = "The output URL. If the output is a matrix, it is always saved in a pickle format along with the name-to-index mapping dictionary")
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
    parser.add_argument("--reduced-dims",
                       default = -1, type = int,
                       help = "If set to a positive value, the output is a reduced normalized DSD with reduced dimensions given by --reduced_dims")
    parser.add_argument("--gamma", default = 1.0, type = float, help = "DSD gamma parameter")
    
    """
    GLIDE arguments
    """
    parser.add_argument("--get-glide-neighbors", default = False, action = "store_true", help = "If set to true, --get_glide_neighbors returns glide neighbors instead of glide matrix")
    parser.add_argument("--glide-neighbors-k", default = 10, type = int, help = "If --get_glide_neighbors is set to true, the code uses --glide_neighbors to decide on the number of neighbors")
    parser.add_argument("--neighbors-return-format", default = "dataframe", choices = ["dataframe",
                                                                                      "graph",
                                                                                      "map"],
                       help = ("This parameter decides the output format for the GLIDE neighbors. If `dataframe` is selected, the code returns output as a panda DataFrame." + 
                       "If `graph` is selected, the code returns output as a networkx graph, otherwise the output is returned as a simple dictionary {NODE: LIST[NODE]}, where LIST[NODE]" + 
                       "is the list of neighbors for the particular node")) 
    parser.add_argument("--alpha", default = 0.1, type = float, help = "GLIDE alpha parameter")
    parser.add_argument("--beta", default = 1000, type = float, help = "GLIDE beta parameter")
    parser.add_argument("--delta", default = 1, type = float, help = "GLIDE delta parameter")
    parser.add_argument("--local", default = "cw", choices = ["cw", "l3"], help = "The local parameter for GLIDE")
    parser.add_argument("--normalize-local", default = True, action = "store_false", help = "If set to False, the local measures are not normalized")
    parser.add_argument("--weighted-local", default = True, action = "store_false", help = ("If set to False, the adjacency matrix is converted to a unweighted form (setting every non-zero elements to 1)" + 
                       "before applying local measures"))
    parser.add_argument("--scale-local", default = False, action = "store_true", help = "If set to True, scales the local measures by their max value before returning it")
    
    return parser.parse_args()


def main():
    args = parse_args()
    def log(msg):
        """
        If verbose, print
        """
        if args.v:
            print(msg)

    network  = args.network

    ## Parsing the network
    log("Parsing network...")
    dfnet    = pd.read_csv(network, delim_whitespace = True, header = None)
    assert(len(dfnet.columns) > 1) and (len(dfnet.columns) <= 3)
    if len(dfnet.columns) == 2:
        dfnet[2] = 1 # set weight = 1
    nodelist = list(set(dfnet[0]).union(dfnet[1]))
    nodemap  = {k:i for i, k in enumerate(nodelist)}
    n        = len(nodemap)
    
    log("Constructing the adjacency matrix...")
    ## Generate the adjacency matrix
    A   = np.zeros((n, n))
    for p, q, w in dfnet.values:
        p_, q_ = [nodemap[r] for r in [p, q]]
        A[p_, q_] = w
        A[q_, p_] = w
        
    # If return dsd matrix instead of the GLIDE Matrices
    if args.return_dsd_emb or args.return_dsd_dist:
        log("Generating DSD matrix...")
        if args.reduced_dims < 0:
            R = dsd.compute_dsd_embedding(A, gamma = args.gamma, is_normalized = args.normalized)
            if args.return_dsd_dist:
                R = squareform(pdist(R))
        else:
            R = dsd.compute_dsd_reduced_embedding(A, dims = args.reduced_dims)
            if args.return_dsd_dist:
                R = squareform(pdist(np.real(R)))
        with open(args.output, "wb") as f:
            pickle.dump([R, nodemap], f)
        return
                
    
    ## Returning GLIDE outputs
    gamma, alpha, beta, delta, gn, local = args.gamma, args.alpha, args.beta, args.delta, args.normalized, args.local
    kwargs = {"weighted": args.weighted_local, "rescale":args.scale_local, "normalize": args.normalize_local, "annotate": nodemap}
    
    if args.get_glide_neighbors:
        log("Getting GLIDE neighbors...")
        kwargs["output_type"] = args.neighbors_return_format
        Out = cg.get_glide_neighbors(A, k = args.glide_neighbors_k, gamma = gamma, 
                                     alpha = alpha, beta = beta, delta = delta, normalize_dsd = gn,
                                     local = local, **kwargs)
        if isinstance(Out, pd.DataFrame):
            Out.to_csv(args.output, sep = "\t", index = None)
        elif isinstance(Out, nx.Graph):
            nx.write_graphml_lxml(Out, args.output)
        else: # Dump it in a pickle format
            with open(args.output, "wb") as op:
                pickle.dump(Out, op)
        return
    
    log("Computing GLIDE matrix...")
    # Returning the GLIDE matrix as a whole
    G = glide.glide(A, gamma = gamma, alpha = alpha, beta = beta, normalize_dsd = gn, local = local, **kwargs)
    with open(args.output, "wb") as f:
        pickle.dump([G, nodemap], f)
    return
            

if __name__ == "__main__":
    main()
        


