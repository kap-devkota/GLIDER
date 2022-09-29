import glidetools.algorithm.dsd as dsd
import glidetools.algorithm.local as lc
from scipy.spatial.distance import pdist, squareform
import numpy as np


"""
Code for GLIDE computation
"""
def glide(A, 
          alpha = 0.1,
          beta  = 1000,
          delta = 1,
          gamma = 1,
          normalize_dsd = False,
          local = "",
          **kwargs):
    
    """
    Given an adjacency matrix, returns a GLIDE matrix representing pairwise GLIDE score
    
    parameters:
    ----------
    A : A n x n numpy matrix
    alpha, beta, delta, gamma => glide parameters: For more information, see :
    normalize_dsd: If set to True, generates the normalized version of DSD embedding
    
    
    Additional kwargs arguments:
    ----------------------------
    You can also provide your own local and global functions for GLIDE
    
    localf: a function that takes in adjacency matrix and returns the local pairwise score
    globalf: a function that takes in adjacency matrix and returns the global pairwise score
    
    G_emb: if the global embedding is precomputed, returns that embedding. Note: It is an embedding, NOT a distance 
    L_sc: If the local scores are precomputed, returns that score
    
    normalize: If set to true, returns the normalized local score
    weighted:  If set to true, returns the weighted_score
    """
    assert (A is not None) and (len(A.shape) == 2) and (A.shape[0] == A.shape[1]) and (local in {"","cw", "l3"})
    if local == "":
        assert ("local" in kwargs) or ("L_sc" in kwargs)

    n, _ = A.shape
    # local score
    L  = None
    Ge = None
    Gd = None
    
    # In precomputed
    if ("G_emb" in kwargs) or ("L_sc" in kwargs):
        if  "G_emb" in kwargs:
            Ge = kwargs["G_emb"]
        if "L_sc" in kwargs:
            L = kwargs["L_sc"]
    
    # If the local and/or global functions are provided
    elif ("localf" in kwargs) or ("globalf" in kwargs):
        if "localf" in kwargs:
            local = kwargs["localf"]
            L = local(A, **kwargs)
        if "globalf" in kwargs:
            global_ = kwargs["globalf"]
            Ge = global_(A, **kwargs)
    
    if L == None:
        lfunc = lc.cw if local == "cw" else lc.l3
        L = lfunc(A, **kwargs)
    
    if Ge == None:
        Ge = dsd.compute_dsd_embedding(A, gamma = gamma, is_normalized = normalize_dsd)
    
    Gd = squareform(pdist(Ge)) + np.eye(n, n)
    
    # Compute the GLIDE score
    R = np.exp(alpha / (1 + beta * Gd)) * L + (delta / Gd)
    return R
    
