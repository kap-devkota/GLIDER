import glidetools.algorithm.glide as glide 
import pandas as pd
import networkx as nx
import numpy as np

def get_glide_neighbors(network, k = 10, **kwargs):
    """
        Given a network, returns the top glide neighbors
        
        parameters:
        -----------
            network: Can be either dataframe or an adjacency matrix. 
                       If the network is a dataframe. columns should be named ["p", "q"] or ["p", "q", "w"]
            k: The top `k` glide neighbors
            
        Additional parameters:
        ---------------------
            kwargs represents the additional parameters for glide
    """
    assert (isinstance(network, np.ndarray)) or (isinstance(network, pd.DataFrame))
    if isinstance(network, pd.DataFrame):
        assert (len(network.columns) >= 2) and (len(network.columns) <= 3)
        if len(network.columns) == 2:
            network["w"] = 1
    
        nodes = list(set(network["p"]).union(network["q"]))
    
        # Forward and reverse map
        nmap  = {k: i for k, i in enumerate(nodes)}
        rnmap = {i: k for k, i in nmap.items()}
    
        n  = len(nodes)
        A  = np.zeros((n, n))
        for p, q, w in network.values:
            p_, q_ = [nmap[p], nmap[q]]
            A[p_, q_] = w
            A[q_, p_] = w 
    else:
        A    = network
        n, _ = A.shape
        nodes = range(n)
        if "annotate" in kwargs and kwargs["annotate"] is not None: 
            nmap  = kwargs["annotate"]
            nodes = list(nmap.keys())
            rnmap = {i:k for k, i in nmap.items()}
        else:
            nmap  = {i:i for i in range(n)}
            rnmap = nmap
    
    # Get the GLIDE matrix
    G = glide.glide(A, **kwargs)
    np.fill_diagonal(G, 0)
    
    # Sort the rows by the GLIDE scores and get the "k" largest entries
    G1 = np.argsort(-G, axis = 1)[:, :k]
    
    # Return the networkx graph 
    if ("output_type" in kwargs) and (kwargs["output_type"] == "graph"):
        edgelist = []
        for i in range(n):
            for j in G1[i]:
                edgelist.append((rnmap[i], rnmap[j], G[i, j]))
        Gnet = nx.Graph()
        Gnet.add_weighted_edges_from(edgelist)
        return Gnet
    
    # Return the dictionary of nodes and their nearest GLIDE neighbors
    neighbors = {}
    for i, node in enumerate(nodes):
        neighbors[node] = [rnmap[j] for j in G1[i]]

    # Return the pandas dataframe
    if ("output_type" in kwargs) and (kwargs["output_type"] == "dataframe"):
        ret = pd.DataFrame(neighbors).T.reset_index()
        ret.columns = ["Node", *[f"Neighbor-{i}" for i in range(k)]]
        return ret
    
    return neighbors




