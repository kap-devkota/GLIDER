from gfunc.algorithm.glide import glide_predict_links
from gfunc.algorithm.computations import compute_embedding


def glide_mat(edgelist, 
              is_annotated = True,
              lamb         = 1,
              is_normalized= False,
              glide_alph   = 0.1,
              glide_beta   = 1000,
              glide_delta  = 1,
              glide_loc    = "cw_normalized"):
    """
    edgelist: list of edges in the format [(p, q, w), ... ]
    is_annotated: True if the graph nodes in the edgelist have string names, False otherwise
    lamb        : Lambda parameter for DSD: Default  = 1
    is_normalized: Is the DSD embedding normalized, default = False
    """
    if is_annotated:
        e_list  = []
        nodemap = {}
        count   = 0
        for e in edgelist:
            p, q, w = e
            for m in [p, q]:
                if m not in nodemap:
                    nodemap[m] = count
                    count     += 1
            e_list.append((nodemap[p], nodemap[q], w))
        edgelist = e_list

    X = compute_embedding(edgelist, lm = lamb, is_normalized = is_normalized)
    gmat = glide_predict_links(edgelist, 
                               X, 
                               params = {
                                   "alpha": glide_alph,
                                   "beta" : glide_beta,
                                   "delta": glide_delta,
                                   "loc"  : glide_loc
                               })
    if is_annotated:
        return gmat, nodemap
    return gmat
