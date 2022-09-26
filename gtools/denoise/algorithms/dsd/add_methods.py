import numpy as np
import json
import networkx as nx

def compute_degree(A, weighted = True):
    e   = np.ones((A.shape[0], 1))
    if weighted:
        deg = A @ e
    else:
        A_  = np.where(A > 0, 1, 0)
        deg = A_ @ e
    return deg.flatten()

def compute_clustering(A, sc):
    G  = nx.from_numpy_matrix(A)
    cc = nx.clustering(G)
    av = 0
    for c in cc:
        av += cc[c]
    av = av / (len(cc) *  sc)
    return cc, av

def weight_dijkstra(A, p, q,
                    f_s_to_d,
                    f_d_to_s,
                    smat):
    ps = f_s_to_d(p)
    qs = f_d_to_s(q)
    rs = qs
    r  = q
    wt = 0
    while rs != ps:
        preds = smat[ps, rs]
        pred  = f_s_to_d(preds)
        wt   += 1 / A[pred, r]
        r     = pred
        rs    = preds
    return 1 / wt

def add_node_Dst(A, nA, ranked_edgelist, l_to_i_src, SDIST_MAT, SDIST_JS, nodes_to_look, no_edges_added, return_added_edgelist = False):
    with open(SDIST_JS, "r") as sp:
        l_to_i_dist = json.load(sp)
    i_to_l_src = {}
    for k in l_to_i_src:
        i_to_l_src[l_to_i_src[k]] = k
    i_to_l_dist = {}
    for k in l_to_i_dist:
        i_to_l_dist[l_to_i_dist[k]] = k
    def id_s_to_d(p):
        return l_to_i_dist[i_to_l_src[p]]
    def id_d_to_s(p):
        return l_to_i_src[i_to_l_dist[p]]
    smat          = np.load(SDIST_MAT)
    n_nodes       = A.shape[0]
    A_added_edges = A.copy()
    added_dict    = {}
    n_ed          = len(ranked_edgelist)
    if return_added_edgelist:
        added_edgelist = []
    for n in nodes_to_look: 
        added_dict[n] = no_edges_added
    counter       = 0
    n_added       = 0
    while len(added_dict) != 0:
        u, v, score = ranked_edgelist[counter]
        counter    += 1
        if counter >= n_ed:
            break
        accepted    = u in added_dict or v in added_dict
        if not accepted:
            continue
        if A_added_edges[u, v] == 0:
            for m in [u, v]:
                if m in added_dict:
                    added_dict[m]    -= 1
                    if added_dict[m] == 0:
                        del added_dict[m]
            n_added            += 1
            wt                  = weight_dijkstra(A, u, v,
                                                  id_s_to_d,
                                                  id_d_to_s,
                                                  smat)
            A_added_edges[u, v] = wt
            A_added_edges[v, u] = wt
            if return_added_edgelist:
                added_edgelist.append((u, v))
    print(f"Percentage added : {n_added/ nA}")
    if return_added_edgelist:
        return A_added_edges, added_edgelist
    return A_added_edges

def add_nodebased_Dst(A, nA, ranked_edgelist, l_to_i_src, SDIST_MAT, SDIST_JS, n_add = 10, d_thres = 25, d_min= 5, between_both = True):
    with open(SDIST_JS, "r") as sp:
        l_to_i_dist = json.load(sp)
    i_to_l_src = {}
    for k in l_to_i_src:
        i_to_l_src[l_to_i_src[k]] = k
    i_to_l_dist = {}
    for k in l_to_i_dist:
        i_to_l_dist[l_to_i_dist[k]] = k
    def id_s_to_d(p):
        return l_to_i_dist[i_to_l_src[p]]
    def id_d_to_s(p):
        return l_to_i_src[i_to_l_dist[p]]
    smat          = np.load(SDIST_MAT)
    A_added_edges = A.copy()
    deg_vec       = compute_degree(A, weighted = False)
    added_dict    = {}
    n_ed          = len(ranked_edgelist)
    for i in range(deg_vec.shape[0]): 
        if deg_vec[i] >= d_min and deg_vec[i] <= d_thres:
            added_dict[i] = n_add
    counter       = 0
    n_added       = 0
    while len(added_dict) != 0:
        u, v, score = ranked_edgelist[counter]
        counter    += 1
        if counter >= n_ed:
            break
        accept      = u not in added_dict or v not in added_dict
        if between_both:
            accept  = u not in added_dict and v not in added_dict
        if accept:
            continue
        if A_added_edges[u, v] == 0:
            for m in [u, v]:
                if m in added_dict:
                    added_dict[m] -= 1
                    if added_dict[m] == 0:
                        del added_dict[m]
            n_added += 1
            ed                  = weight_dijkstra(A, u, v,
                                                  id_s_to_d,
                                                  id_d_to_s,
                                                  smat)
            A_added_edges[u, v] = ed
            A_added_edges[v, u] = ed
    print(f"Percentage added : {n_added/ nA}")
    return A_added_edges


def add_nodebased_Res(A, nA, ranked_edgelist, l_to_i_src, INVL_MAT, INVL_JS, n_add = 10, d_thres = 10, between_both = True):
    with open(INVL_JS, "r") as ip:
        l_to_i_invl = json.load(ip)
    i_to_l_src = {}
    for k in l_to_i_src:
        i_to_l_src[l_to_i_src[k]] = k
    i_to_l_invl = {}
    for k in l_to_i_invl:
        i_to_l_invl[l_to_i_invl[k]] = k
    def id_s_to_i(p):
        return l_to_i_invl[i_to_l_src[p]]
    def id_i_to_s(p):
        return l_to_i_src[i_to_l_invl[p]]
    def e_f(i):
        e_i = np.zeros((A.shape[0], 1))
        e_i[i, 0] = 1
        return e_i
    imat        = np.load(INVL_MAT)
    A_added_edges      = A.copy()
    if d_thres <= 0:
        deg_vec, d_thres = compute_clustering(A, -d_thres)
    else:
        deg_vec       = compute_degree(A, weighted = False)
    added_dict    = {}
    for i in range(deg_vec.shape[0]): 
        if deg_vec[i] <= d_thres:
            added_dict[i] = n_add
    counter       = 0
    n_added       = 0
    while len(added_dict) != 0:
        u, v, score = ranked_edgelist[counter]
        counter    += 1
        accept      = u not in added_dict or v not in added_dict
        if between_both:
            accept  = u not in added_dict and v not in added_dict
        if accept:
            continue
        if A_added_edges[u, v] == 0:
            for m in [u, v]:
                if m in added_dict:
                    added_dict[m] -= 1
                    if added_dict[m] == 0:
                        del added_dict[m]
            n_added += 1
            u1   = id_s_to_i(u)
            v1   = id_s_to_i(v)
            X_uv = e_f(u1) - e_f(v1)
            Re   = (X_uv.T @ imat) @ X_uv
            A_added_edges[v, u] = 1.0 / Re[0,0]
            A_added_edges[u, v] = 1.0 / Re[0,0]
        counter += 1
    print(f"Percentage Added {n_added/nA}")
    return A_added_edges
    

def add_global_Dst(A, nA, ranked_edgelist, l_to_i_src, SDIST_MAT, SDIST_JS, p_add = 0.1, return_added_edgelist = False):
    with open(SDIST_JS, "r") as sp:
        l_to_i_dist = json.load(sp)
    i_to_l_src = {}
    for k in l_to_i_src:
        i_to_l_src[l_to_i_src[k]] = k
    i_to_l_dist = {}
    for k in l_to_i_dist:
        i_to_l_dist[l_to_i_dist[k]] = k
    def id_s_to_d(p):
        return l_to_i_dist[i_to_l_src[p]]
    def id_d_to_s(p):
        return l_to_i_src[i_to_l_dist[p]]
    smat          = np.load(SDIST_MAT)
    A_added_edges = A.copy()
    n_ed          = len(ranked_edgelist)
    counter, n_added = 0, 0
    if return_added_edgelist:
        added_edgelist = []
    while n_added <= p_add * nA:
        u, v, score = ranked_edgelist[counter]
        counter    += 1
        if counter >= n_ed:
            break
        if A_added_edges[u, v] == 0:
            wt                  = weight_dijkstra(A, u, v,
                                                  id_s_to_d,
                                                  id_d_to_s,
                                                  smat)
            A_added_edges[u, v] = wt
            A_added_edges[v, u] = wt
            n_added            += 1
            if return_added_edgelist:
                added_edgelist.append((u, v))
    if return_added_edgelist:
        return A_added_edges, added_edgelist
    return A_added_edges


def add_global_Res(A, nA, ranked_edgelist, l_to_i_src, INVL_MAT, INVL_JS, p_add = 0.1):
    with open(INVL_JS, "r") as ip:
        l_to_i_invl = json.load(ip)
    i_to_l_src = {}
    for k in l_to_i_src:
        i_to_l_src[l_to_i_src[k]] = k
    i_to_l_invl = {}
    for k in l_to_i_invl:
        i_to_l_invl[l_to_i_invl[k]] = k
    def id_s_to_i(p):
        return l_to_i_invl[i_to_l_src[p]]
    def id_i_to_s(p):
        return l_to_i_src[i_to_l_invl[p]]
    def e_f(i):
        e_i = np.zeros((A.shape[0], 1))
        e_i[i, 0] = 1
        return e_i
    imat          = np.load(INVL_MAT)
    A_added_edges = A.copy()
    n_added, counter = 0, 0
    while n_added < nA * p_add:
        u, v, score = ranked_edgelist[counter]
        counter    += 1
        if A_added_edges[u, v] == 0:
            u1   = id_s_to_i(u)
            v1   = id_s_to_i(v)
            X_uv = e_f(u1) - e_f(v1)
            Re   = (X_uv.T @ imat) @ X_uv
            A_added_edges[u, v] = 1.0 / Re[0,0]
            A_added_edges[v, u] = 1.0 / Re[0,0]
            n_added += 1
    return A_added_edges
    
def add_nodebased_DSD(A, nA, ranked_edgelist, dsd_dist, min_deg = 5, max_deg = 20, n_add = 25, between_both = True):    
    A_added_edges      = A.copy()
    deg_vec       = compute_degree(A, weighted = False)
    added_dict    = {}
    for i in range(deg_vec.shape[0]): 
        if deg_vec[i] >= min_deg and deg_vec[i] <= max_deg:
            added_dict[i] = n_add
    counter       = 0
    n_ranked      = len(ranked_edgelist)
    n_added       = 0
    while len(added_dict) != 0:
        u, v, score = ranked_edgelist[counter]
        counter    += 1
        if counter >= n_ranked:
            break
        accept      = u not in added_dict or v not in added_dict
        if between_both:
            accept  = u not in added_dict and v not in added_dict
        if accept:
            continue
        if A_added_edges[u, v] == 0:
            for m in [u, v]:
                if m in added_dict:
                    added_dict[m] -= 1
                    if added_dict[m] == 0:
                        del added_dict[m]
            n_added += 1
            A_added_edges[u, v] = 1.0 / dsd_dist[u, v] if dsd_dist[u, v] > 1.0 else 1
            A_added_edges[v, u] = A_added_edges[u, v]
    print(f"Percentage added {n_added / nA}")
    return A_added_edges


def add_node_tri_Dst(A, nA, ranked_edgelist, l_to_i_src, SDIST_MAT, SDIST_JS, n_add = 10, d_thres = 25, d_min= 5):
    with open(SDIST_JS, "r") as sp:
        l_to_i_dist = json.load(sp)
    i_to_l_src = {}
    for k in l_to_i_src:
        i_to_l_src[l_to_i_src[k]] = k
    i_to_l_dist = {}
    for k in l_to_i_dist:
        i_to_l_dist[l_to_i_dist[k]] = k
    def id_s_to_d(p):
        return l_to_i_dist[i_to_l_src[p]]
    def id_d_to_s(p):
        return l_to_i_src[i_to_l_dist[p]]
    smat          = np.load(SDIST_MAT)
    A_added_edges = A.copy()
    deg_vec       = compute_degree(A, weighted = False)
    added_dict    = {}
    n_ed          = len(ranked_edgelist)
    for i in range(deg_vec.shape[0]): 
        if deg_vec[i] >= d_min and deg_vec[i] <= d_thres:
            added_dict[i] = {"neighbors" : [], "len" : 0}  
    counter       = 0
    n_a           = 0
    while n_a < len(added_dict) * n_add:
        u, v, score = ranked_edgelist[counter]
        counter    += 1
        if counter >= n_ed:
            break 
        n_accept  = u not in added_dict or v not in added_dict
        if n_accept:
            continue
        if A_added_edges[u, v] == 0:
            if added_dict[u]["len"] < n_add:
                added_dict[u]["len"] += 1
                n_a += 1
                added_dict[u]["neighbors"] += [v]
            if added_dict[v]["len"] < n_add:
                added_dict[v]["len"] += 1
                n_a += 1
                added_dict[v]["neighbors"] += [u]

    edges_to_create = {}
    n_added = 0
    for u in added_dict:
        neigh = added_dict[u]["neighbors"]
        for n in neigh:
            if (n, u) not in edges_to_create and u != n:
                edges_to_create[(u, n)] = True
                n_added += 1
                for nn in added_dict[n]["neighbors"]:
                    if (nn, u) not in edges_to_create and u != nn:
                        edges_to_create[(u, nn)] = True
                        n_added += 1
    for key in edges_to_create:
        u, v                = key
        A_added_edges[u, v] = weight_dijkstra(A, u, v,
                                                  id_s_to_d,
                                                  id_d_to_s,
                                                  smat)
        A_added_edges[v, u] = A_added_edges[u, v]
    print(f"Percentage Added {n_added/nA}")
    return A_added_edges
