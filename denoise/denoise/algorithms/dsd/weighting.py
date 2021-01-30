import networkx as nx
import numpy as np

def constant_reweight(new_edges, params):
    """
    Takes in a list of `new_edges` of elements `(p, q, wt)` and returns 
    `(p, q, params["wt"])`
    @param new_edges : A list of form `(p, q, _)`, with p and q the nodes
    @param params    : A dictionary with one key "wt", which is a constant 
                       weight to be applied to every edges in `new_edges`.
    @return          : A list with all edges in `new_edges` reweighted to 
                       params["wt"]
    """
    ret_list    = []
    try:
        cons_wt = params["wt"]
    except:
        raise Exception("[x] `params` Dictionary Illegal!")
    for ed in new_edges:
        p, q, _ = ed
        ret_list.append((p, q, cons_wt))
    return ret_list


def linear_reweight(new_edges, params):
    """
    Perform linear reweighting on the `new_edges`.
    @param new_edges: A list of form `(p, q, _)`, with p and q as nodes.
    @param params   : A dictionary with keys; 
    {
        "max_weight"  => Initial (max) weight. Should be greater than 0.
        "spacing"     => After how many edges should the weight be reduced.
                         If spacing = 10, the same weight appled to the 10 
                         consecutive edges, before a linear operation is 
                         applied.
        "red_val"     => The value to be reduced after a spacing. If `red_val` is 0.05,
                         and `spacing` is 10, the weight is reduced by 0.05 after
                         10 consecutive runs.
        "min_weight"  => What is the minimum weight given to the edges. Should 
                         be greater than 0 and less than `max_weight`.
    }
    @return          : A list with updated weight edges.
    """
    ret_list    = []
    try:
        max_w   = params["max_weight"]
        spacing = params["spacing"]
        spacing = int(spacing)
        if (spacing < 1):
            spacing = 1
        red_val = params["red_val"]
        min_w   = params["min_weight"]
    except:
        raise Exception("[x] `params` Dictionary Illegal!")
    curr_w      = max_w
    count       = 1
    for ed in new_edges:
        p, q, _ = ed
        if (curr_w > min_w) and (count % spacing == 0):
            curr_w = curr_w - red_val
        ret_list.append((p, q, curr_w))
        count   += 1
    return ret_list
            
        

def geometric_reweight(new_edges, params):
    """
    Perform geometric reweighting on the `new_edges`.
    @param new_edges: A list of form `(p, q, _)`, with p and q as nodes.
    @param params   : A dictionary with keys; 
    {
        "max_weight"  => Initial (max) weight. Should be greater than 0.
        "spacing"     => After how many edges should the weight be reduced.
                         If spacing = 10, the same weight appled to the 10 
                         consecutive edges, before a linear operation is 
                         applied.
        "factor"      => > 1
                         ; If scale is 2 and factor is 4, after k update,
                         ; the new weight becomes `new_weight` = 2 * `max` * 1 / (4)^k 
        "min_weight"  => What is the minimum weight given to the edges. Should 
                         be greater than 0 and less than `max_weight`.
    }
    @return          : A list with updated weight edges.
    """
    ret_list    = []
    try:
        max_w   = params["max_weight"]
        spacing = params["spacing"]
        spacing = int(spacing)
        if (spacing < 1):
            spacing = 1
        factor  = params["factor"]
        min_w   = params["min_weight"]
    except:
        raise Exception("[x] `params` Dictionary Illegal!")
    curr_w      = max_w
    count       = 1
    for ed in new_edges:
        p, q, _ = ed
        if (curr_w > min_w) and (count % spacing == 0):
            curr_w = curr_w / (factor)
        ret_list.append((p, q, curr_w))
        count   += 1
    return ret_list        

def inverse_reweight(new_edges, params):
    """
    Perform geometric reweighting on the `new_edges`.
    @param new_edges: A list of form `(p, q, _)`, with p and q as nodes.
    @param params   : A dictionary with keys; 
    {
        "max_weight"  => Initial (max) weight. Should be greater than 0.
        "spacing"     => After how many edges should the weight be reduced.
                         If spacing = 10, the same weight appled to the 10 
                         consecutive edges, before a linear operation is 
                         applied.
        "scale"       => > 0
                         ; If scale is 2, after k update,
                         ; the new weight becomes `new_weight` = scale * `max_wt` * 1 / k 
        "min_weight"  => What is the minimum weight given to the edges. Should 
                         be greater than 0 and less than `max_weight`.
    }
    @return          : A list with updated weight edges.
    """
    ret_list    = []
    try:
        max_w   = params["max_weight"]
        spacing = params["spacing"]
        spacing = int(spacing)
        if (spacing < 1):
            spacing = 1
        scale   = params["scale"]
        factor  = params["factor"]
        min_w   = params["min_weight"]
    except:
        raise Exception("[x] `params` Dictionary Illegal!")
    curr_w      = max_w
    count       = 1
    for ed in new_edges:
        p, q, _ = ed
        if (curr_w > min_w) and (count % spacing == 0):
            curr_w = scale * max_w / (count / spacing)
        ret_list.append((p, q, curr_w))
        count   += 1
    return ret_list        


def dijkstra_reweight(new_edges, A_mat):
    """
    Function that performs dijkstra reweighting.
    """
    e_inv_f  = np.vectorize(lambda l : 0.0 if l == 0.0 else 1.0 / l)
    A_mat_i  = e_inv_f(A_mat)
    G        = nx.from_numpy_matrix(A_mat_i)
    r_d      = nx.floyd_warshall(G)
    ret_list = []
    for e in new_edges:
        p, q, _ = ed
        wt      = 1 / (r_d[p][q])
        ret_list.append((p, q, wt))
    return ret_list
