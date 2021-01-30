import numpy as np
import json
import networkx as nx

def sparsify(A, directed = False, node_map = None):
    """Given an adjacency matrix as a numpy matrix, returns the
    sparsified form of the matrix (or adjacency list).
    """
    dim = A.shape[0]
    edgelist = []
    for i in range(dim):
        for j in range(i + 1):
            
            if i <= j:
                continue

            if A[i, j] != 0:
                edgelist.append((i, j, A[i, j]))

            if directed:
                if A[j, i] != 0:
                    edgelist.append((j, i, A[j, i]))
    if node_map is not None:
        r_node_map = {}
        for k in node_map:
            r_node_map[node_map[k]] = k
        a_edgelist = []
        for (p, q, w) in edgelist:
            a_edgelist.append(r_node_map[p], r_node_map[q], w)
        return a_edgelist, dim
    else:
        return edgelist, dim

def densify(edgelist, dim = None, directed = False):
    """Given an adjacency list for the graph, computes the adjacency
    matrix.
    """
    if dim is None:
        dim = get_dim(edgelist)

    A = np.zeros((dim, dim), dtype = np.double)
    
    for edge in edgelist:
        p, q, wt = edge
        A[p, q] = wt

        if not directed:
            A[q, p] = wt

    return A

def get_dim(edgelist):
    """Given an adjacency list for a graph, returns the number of nodes in
    the graph.
    """
    node_dict = {}
    node_count = 0

    for edge in edgelist:
        p, q = edge[ :2]

        if p not in node_dict:
            node_dict[p] = True
            node_count += 1

        if q not in node_dict:
            node_dict[q] = True
            node_count += 1

    return node_count

def add_random_numbering(edgelist):
    """Adds a random number to each of the edges in the adjacency list,
    for randomization.
    """
    no_edges = len(edgelist)
    perm = np.random.permutation(no_edges)
    _edgelist = []
    for i in range(no_edges):
        _edgelist.append((edgelist[i][0], edgelist[i][1], edgelist[i][2], perm[i]))

    return _edgelist

def compute_reduced_graph(edgelist, dims = None, p_reduced = 0.1):
    """Given an adjacency list `edgelist`, splits `edgelist` into two. One of the
    subgraph should be connected, have all the nodes in `edgelist` and
    should contain `(1 - p_reduced) * edgelist` number of edges.
    """
    if dims == None:
        _edgelist = add_random_numbering(edgelist)

    no_edges = len(edgelist)
    node_dict = {}
    nodes_added = 0

    rA = []

    _edgelist = sorted(_edgelist, key = lambda x: x[3])

    # Automatically add the first edge

    p, q, weight, _ = _edgelist[0]
    rA.append((p, q, weight))

    node_dict[p] = True
    node_dict[q] = True
    nodes_added = 2
    extras = []

    # Remove the first element
    _edgelist = _edgelist[1: ]

    refresh = 0
    connected = True

    while (len(_edgelist) != 0):
        
        elem = _edgelist.pop(0)
        p, q, wt, _ = elem

        if p not in node_dict and q not in node_dict:
            _edgelist.append(elem)
            refresh += 1
            if refresh >= len(_edgelist):
                connected = False
                break
        else:
            refresh = 0
            if p in node_dict and q not in node_dict:
                print((p, q, wt))
                rA.append((p, q, wt))
                node_dict[q] = True
                nodes_added += 1
            elif p not in node_dict and q in node_dict:
                print((p, q, wt))
                rA.append((p, q, wt))
                node_dict[p] = True
                nodes_added += 1
            else:
                extras.append((p, q, wt))

    if connected == False:
        return None

    add_length = int((1 - p_reduced) * no_edges - nodes_added + 1) 
    if add_length > len(extras):
        add_length = len(extras)
            
    rA = rA + extras[: add_length]
    if add_length == extras:
        rrA = None
    else:
        rrA = extras[add_length: ]

    return rA, rrA

def get_connected_components(edgelist):
    """This function takes in an adjacency list and returns one of the
    connected components from the list
    """
    _edgelist = add_random_numbering(edgelist)
    node_dict = {}
    nodes_added = 0

    cc = []
    _edgelist = sorted(_edgelist, key = lambda x: x[3])

    # Automatically add the first edge
    p, q, weight, _ = _edgelist[0]
    cc.append((p, q, weight))
    node_dict[p] = True
    node_dict[q] = True
    nodes_added = 2

    # Remove the first element
    _edgelist = _edgelist[1: ]
    refresh = 0

    while (len(_edgelist) != 0):
        elem = _edgelist.pop(0)
        p, q, wt, _ = elem
        if p not in node_dict and q not in node_dict:
            _edgelist.append(elem)
            refresh += 1
            if refresh >= len(_edgelist):
                break
        else:
            refresh = 0
            cc.append((p, q, wt))
            if p not in node_dict:
                node_dict[p] = True
                nodes_added += 1
            if q not in node_dict:
                node_dict[q] = True
                nodes_added += 1
    return cc
        

def compute_graph_reduced_nodes(edgelist, no_nodes):
    """
    Given an adjacency list, returns a connected subgraph from the adjacency list containing all the edge 
    connections, with randomly selected nodes of size `no_nodes`.
    """
    _edgelist = add_random_numbering(edgelist)
    no_edges = len(edgelist)
    node_dict = {}
    nodes_added = 0
    rA = []
    _edgelist = sorted(_edgelist, key = lambda x: x[3])

    # Automatically add the first edge
    p, q, weight, _ = _edgelist[0]
    rA.append((p, q, weight))

    node_dict[p] = True
    node_dict[q] = True
    nodes_added = 2
    extras = []

    # Remove the first element
    _edgelist = _edgelist[1: ]
    refresh = 0
    while (len(_edgelist) != 0):
        elem = _edgelist.pop(0)
        p, q, wt, _ = elem
        
        if nodes_added < no_nodes:
            if p not in node_dict and q not in node_dict:
                _edgelist.append(elem)
                refresh += 1
                if refresh >= len(_edgelist):
                    break
            elif p in node_dict and q not in node_dict:
                rA.append((p, q, wt))
                print((p, q, wt))
                node_dict[q] = True
                nodes_added += 1
            elif p not in node_dict and q in node_dict:
                rA.append((p, q, wt))
                print((p, q, wt))
                node_dict[p] = True
                nodes_added += 1
            else:
                print((p, q, wt))
                rA.append((p, q, wt))
        else:
            if p in node_dict and q in node_dict:
                rA.append((p, q, wt))
    return rA


def compute_clustering_coeff(edgelist):
    """Create a datastructure that is a dictionary of dictionaries,
    indexed by node label, and uses it to compute the clustering
    coeffecient.
    """
    edge_dict = {}
    node_dict = {}    
    for e in edgelist:
        p, q, wt = e
        if p == q:
            continue
        if p not in node_dict:
            node_dict[p] = {}
            node_dict[p]["neighbors"] = {}
            node_dict[p]["n_list"]    = []
            node_dict[p]["degree"]    = 0
        
        if q not in node_dict:
            node_dict[q] = {}
            node_dict[q]["neighbors"] = {}
            node_dict[q]["n_list"]    = []
            node_dict[q]["degree"]    = 0

        if (p, q) not in edge_dict and (q, p) not in edge_dict:
            edge_dict[(p, q)]             = True
            node_dict[p]["neighbors"][q]  = True
            node_dict[q]["neighbors"][p]  = True
            
            node_dict[p]["n_list"].append(q)
            node_dict[q]["n_list"].append(p)
            
            node_dict[p]["degree"]        += 1
            node_dict[q]["degree"]        += 1
    
    c_coeff  = 0.0
    n_deg_2  = 0
    for n in node_dict:
        deg_n = node_dict[n]["degree"]
        if deg_n == 1:
            continue
        n_deg_2  += 1
        
        node_dict[n]["delta"] = 0
        for i in range(deg_n):
            for j in range(i):
                if node_dict[n]["n_list"][i] in node_dict[node_dict[n]["n_list"][j]]["neighbors"]:
                    node_dict[n]["delta"] += 1
                        
        r_n      = deg_n * (deg_n - 1) / 2
        c_coeff += node_dict[n]["delta"] / float(r_n)
                    
    return c_coeff / n_deg_2
                
############################DIAMETER CALCULATION#######################################
def compute_shortest_separation(adjacency_dict, node, computed_dict = None):
    """
    Takes in an adjacency dict, which is a dictionary of dictionaries, a node as a key to the dict
    and a already computed dict, returns the dictionary whose key is the tuple of two nodes, and 
    the value is the distance between nodes.
    """
    nodes = list(adjacency_dict.keys())
    len_nodes = len(nodes)
    accounted_for_nodes = {node : True}

    return_dict = {(node, node) : 0}
    nodes_at_curr_depth = {node: True}

    next_depth = 1

    while(not (len(accounted_for_nodes) == len_nodes)):
        nodes_at_next_depth = {}

        for nd in nodes_at_curr_depth:
            for _nd in adjacency_dict[nd]:

                if _nd not in accounted_for_nodes:
                    nodes_at_next_depth[_nd] = True
                    accounted_for_nodes[_nd] = True
                    
                    return_dict[(node, _nd)] = next_depth

        nodes_at_curr_depth = nodes_at_next_depth
        next_depth += 1

    if computed_dict is not None:
        computed_dict.update(return_dict)
        return computed_dict

    return return_dict

def compute_graph_separation_dict(adjacency_dict):
    """Given a dictionary that represents two nodes with edges between
    them, returns the dictionary whose key is the tuple of two nodes,
    and the value is the distance between nodes.
    """
    nodes = list(adjacency_dict.keys())

    computed_dict = None

    for node in nodes:
        if computed_dict == None:
            computed_dict = compute_shortest_separation(adjacency_dict, node, computed_dict = None)
        else:
            computed_dict = compute_shortest_separation(adjacency_dict, node, computed_dict = computed_dict)

    return computed_dict

def convert_to_adjacency_dict(edge_connections):
    """Given an edge list, returns the adjacency dict `adict` where
    `adict[p][q]` is true if p and q are nodes and p and q have edges
    between them.
    """
    adjacency_dict = {}
    
    for edge in edge_connections:
        e1, e2 = edge[ : 2]
        
        if e1 not in adjacency_dict:
            adjacency_dict[e1] = {}

        if e2 not in adjacency_dict:
            adjacency_dict[e2] = {}

        adjacency_dict[e1][e2] = True
        adjacency_dict[e2][e1] = True

    return adjacency_dict

def compute_graph_distances(edge_connections):
    adjdict = convert_to_adjacency_dict(edge_connections)
    return compute_graph_separation_dict(adjdict)


def compute_shortest_paths(edgelist):
    """ 
    """
    G               = nx.Graph()
    G.add_weighted_edges_from(edgelist)
    predecessors, _ = nx.floyd_warshall_predecessor_and_distance(G)
    n_nodes         = len(G)
    p_matrix        = np.zeros((n_nodes, n_nodes), dtype = np.int)
    for src in predecessors:
        for intr in predecessors[src]:
            p_matrix[src, intr] = predecessors[src][intr]
    return p_matrix
