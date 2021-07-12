from .computations import compute_embedding, compute_degree_mat, compute_laplacian, compute_X_normalized
from denoise.graph.operations import densify
from numpy.linalg  import norm

import numpy as np
import itertools
import scipy.spatial.distance as spatial
import ctypes

def rank_edges(edgelist, X):
    """Ranks the edges in the edgelist according to their L2 distance from
    each other in the embedding.
    """
    existing_edgelist = []
    best_known_list    = []
    for ed in edgelist:
        p, q, wt = ed
        _norm    = norm(X[p] - X[q])
        best_known_list.append((p, q, wt, _norm))
    return sorted(best_known_list, key = lambda l : l[3])

    
def predict_links(X, metric="euclidean"):
    """Predicts the most likely links in a graph given an embedding X
    of a graph.

    Returns a ranked list of (edges, distances) sorted from closest to
    furthest.
    """
    distances = spatial.pdist(X, metric=metric)
    n = X.shape[0]
    edges = itertools.combinations(range(n), 2) # generator expression doesn't actualize list :)
    edges_and_distances = list(zip(edges, distances))
    edges_and_distances.sort(key=lambda x: x[1])
    return edges_and_distances


###################################################################################################################################################
def create_edge_dict(edgelist):
    """
    Creates an edge dictionary with the edge `(p, q)` as the key, and weight `w` as the value.
    @param  edgelist -> A list with elements of form `(p, q, w)`
    @return edgedict -> A dictionary with key `(p, q)` and value `w`.
    """
    edgedict             = {}
    for (p, q, w) in edgelist:
        edgedict[(p, q)] = w
    return edgedict

def create_neighborhood_dict(edgelist):
    """
    Create a dictionary with nodes as key and a list of neighborhood nodes as the value
    @param edgelist          -> A list with elements of form `(p, q, w)`
    @param neighborhood_dict -> A dictionary with key `p` and value, a set `{p1, p2, p3, ...}`
    """
    ndict                = {}
    for ed in edgelist:
        p, q, _          = ed
        if p not in ndict:
            ndict[p]     = set()
        if q not in ndict:
            ndict[q]     = set()
        ndict[p].add(q)
        ndict[q].add(p)
    return ndict

def compute_cw_score(p, q, edgedict, ndict, params = None):
    """
    Computes the common weighted score between p and q
    @param p        -> A node of the graph
    @param q        -> Another node in the graph
    @param edgedict -> A dictionary with key `(p, q)` and value `w`.
    @param ndict    -> A dictionary with key `p` and the value a set `{p1, p2, ...}`
    @param params   -> Should always be none here
    @return         -> A real value representing the score
    """
    if (len(ndict[p]) > len(ndict[q])):
        temp  = p
        p     = q
        q     = temp            
    score     = 0
    for elem in ndict[p]:
        if elem in ndict[q]:
            p_elem  = edgedict[(p, elem)] if (p, elem) in edgedict else edgedict[(elem, p)]
            q_elem  = edgedict[(q, elem)] if (q, elem) in edgedict else edgedict[(elem, q)]
            score  += p_elem + q_elem
    return score

def compute_cw_score_normalized(p, q, edgedict, ndict, params = None):
    """
    Computes the common weighted normalized score between p and q
    @param p        -> A node of the graph
    @param q        -> Another node in the graph
    @param edgedict -> A dictionary with key `(p, q)` and value `w`.
    @param ndict    -> A dictionary with key `p` and the value a set `{p1, p2, ...}`
    @param params   -> Should always be none here
    @return         -> A real value representing the score
    """
    if (len(ndict[p]) > len(ndict[q])):
        temp  = p
        p     = q
        q     = temp            
    score     = 0
    for elem in ndict[p]:
        if elem in ndict[q]:
            p_elem  = edgedict[(p, elem)] if (p, elem) in edgedict else edgedict[(elem, p)]
            q_elem  = edgedict[(q, elem)] if (q, elem) in edgedict else edgedict[(elem, q)]
            score  += p_elem + q_elem
    degrees  = params["deg"]
    return score / np.sqrt(degrees[p] * degrees[q])

def compute_l3_unweighted_mat(A):
    A_u  = np.where(A>0, 1, 0)
    d, _ = A_u.shape 
    e    = np.ones((d, 1))
    deg  = A_u @ e
    ideg = np.where(deg > 0, 1 / deg, 0)
    sdeg = np.diag(np.sqrt(ideg).flatten())
    A1   = sdeg @ A_u @ sdeg
    return A1

def compute_l3_weighted_mat(A):
    d, _ = A.shape
    e    = np.ones((d, 1))
    deg  = A @ e
    ideg = np.where(deg > 0, 1 / deg, 0)
    sdeg = np.diag(np.sqrt(ideg).flatten())
    A1   = sdeg @ A @ sdeg        
    return A1

def compute_l3_score_mat(p, q, edgedict, ndict, params = None):
    L3 = params["l3"]
    return L3[p, q]

def compute_degree_vec(edgelist):
    A   = densify(edgelist)
    e   = np.ones((A.shape[0], 1))
    deg = A @ e
    return deg.flatten()

###################################################################################################################################################

def glide_predict_links(edgelist, X, params={}):
    """Predicts the most likely links in a graph given an embedding X
    of a graph.
    Returns a ranked list of (edges, distances) sorted from closest to 
    furthest.
    @param edgelist -> A list with elements of type `(p, q, wt)`
    @param X        -> A nxk embedding matrix
    @param params   -> A dictionary with entries 
    {
        alpha       => real number
        beta        => real number
        delta       => real number
        loc         => String, can be `cw` for common weighted, `l3` for l3 local scoring
    
        ### To enable ctypes, the following entries should be there ###
        
        ctypes_on   => True  # This key should only be added if ctypes is on (dont add this
                           # if ctypes is not added)
        so_location => String location of the .so dynamic library
    }
    """
    edgedict      = create_edge_dict(edgelist)
    ndict         = create_neighborhood_dict(edgelist)
    params_        = {}
        
    # Embedding
    pairwise_dist = spatial.squareform(spatial.pdist(X))
    N             = X.shape[0]

    # Choose Scoring Method
    if "method" not in params or params["method"] == "sum":
        method        = "sum"
        alpha         = params["alpha"]
        beta          = params["beta"]
        delta         = params["delta"]
    elif params["method"] == "prod":
        method        = "prod"
        p_delta       = params["delta"] if "delta" in params else 0.1

    local_metric  = params["loc"]
    if local_metric == "l3_u" or local_metric == "l3":
        A         = densify(edgelist)
        L3        = compute_l3_unweighted_mat(A)
        params_["l3"] = L3
        local_metric  = compute_l3_score_mat
    elif local_metric == "l3_w":
        A         = densify(edgelist)
        L3        = compute_l3_weighted_mat(A)
        params_["l3"] = L3
        local_metric  = compute_l3_score_mat  
    elif local_metric == "cw":
        local_metric = compute_cw_score
    elif local_metric == "cw_normalized":
        params_["deg"]  = compute_degree_vec(edgelist)
        local_metric    = compute_cw_score_normalized
    else:
        raise Exception("[x] The local scoring metric is not available.")
    
    edgelist_with_scores = []
    for i in range(N):
        for j in range(i):
            local_score = local_metric(i, j, edgedict, ndict, params_)
            dsed_dist   = pairwise_dist[i, j]
            if method  == "sum":
                glide_score = (np.exp(alpha / (1 + beta * dsed_dist)) * local_score
                               + delta * 1 / dsed_dist)
            elif method == "prod":
                glide_score = (p_delta + local_metric(i, j, edgedict, ndict, params_)) / dsed_dist
            edgelist_with_scores.append((i, j, glide_score))
    return sorted(edgelist_with_scores, key = lambda l : l[2], reverse = True)    
    

def glide_remove_links(edgelist, params):
    edgedict      = create_edge_dict(edgelist)
    ndict         = create_neighborhood_dict(edgelist)
    params_ = {}        
    
    # Compute Embedding
    A_            = densify(edgelist)
    D_            = compute_degree_mat(A_)
    X             = compute_X_normalized(A_, D_)
    N             = X.shape[0]
    alpha         = params["alpha"]
    local_metric  = params["loc"]
    beta          = params["beta"]
    delta         = params["delta"]
    no_removed    = int(len(edgelist) * (params["removed_frac"]))
    if local_metric == "l3_u" or local_metric == "l3":
        A         = densify(edgelist)
        L3        = compute_l3_unweighted_mat(A)
        params_["l3"] = L3
        local_metric  = compute_l3_score_mat
    elif local_metric == "l3_w":
        A         = densify(edgelist)
        L3        = compute_l3_weighted_mat(A)
        params_["l3"] = L3
        local_metric  = compute_l3_score_mat
    elif local_metric == "cw":
        local_metric  = compute_cw_score
    elif local_metric == "cw_normalized":
        params_["deg"]  = compute_degree_vec(edgelist)
        local_metric    = compute_cw_score_normalized
    elif "local_metric_function" in params:
        local_metric    = params["local_metric_function"]
    else:
        raise Exception("[x] The local scoring metric is not available.")
    updated_edgelist = []
    for ed in edgelist:
        p, q, w     = ed
        p           = int(p)
        q           = int(q)
        w           = float(w)
        local_score = local_metric(p, q, edgedict, ndict, params_)
        dsed_dist   = norm(X[p] - X[q])
        glide_score = (np.exp(alpha / (1 + beta * dsed_dist)) * local_score
                       + delta * 1 / dsed_dist)
        updated_edgelist.append((p, q, glide_score, w))
    updated_edgelist = sorted(updated_edgelist, 
                              key = lambda l: l[2], 
                              reverse = True) 
    ret_edgelist     = return_connected_edgelist(updated_edgelist, no_removed)
    return ret_edgelist

def return_connected_edgelist(edgelist, no_removed):
    p, q, r, w  = edgelist.pop(0)
    n_added     = {}
    n_added[p]  = True
    n_added[q]  = True
    tree  = [(p, q, w)]
    extra = []
    while(True):
        edgelist_top = []
        if len(edgelist) == 0:
            break
        while(len(edgelist) > 0): 
            p, q, r, w = edgelist.pop(0)
            if (p in n_added and q in n_added):
                extra.append((p, q, w))
            elif (p in n_added or q in n_added):
                tree.append((p, q, w))
                for m in [p, q]:
                    if m not in n_added:
                        n_added[m] = True
                break
            else:
                edgelist_top.append((p, q, r, w))
        edgelist = edgelist_top + edgelist
    if no_removed > len(extra):
        return tree
    else:
        ext   = extra[: -no_removed]
        return tree + ext

        
        
    
