import numpy as np

def add_edges(A, ranked_edgelist, alpha):
    """Denoises a graph by adding edges to it.

    Inputs:
    - A: an adjacency matrix representing the graph to add edges to.
    - ranked_edgelist: a list of edges of the form (u, v, w) to add to the
    graph ranked from most to least likely.
    - alpha: a parameter in (0, 1] denoting the proportion of edges to add to
    the graph.

    Returns: a new adjacency matrix representing the graph
    with a portion of the edges added to it.
    """
    A_added_edges = A.copy()
    num_added, counter = 0, 0
    while num_added < len(ranked_edgelist) * alpha:
        u, v, w = ranked_edgelist[counter]
        if A_added_edges[u, v] == 0:
            A_added_edges[u, v] = w
            num_added += 1
        counter += 1
    return A_added_edges

def lsr_convert(ranked_edgelist, edgelist):
    """Converts a ranked edgelist of scores to a ranked edgelist of
    weights using a least squares linear regression model

    Inputs:
    - ranked_edgelist: a list of edges of the form (u, v, s) ranked from
    most to least likely. 
    - edgelist: the original graphs edgelist

    Returns: a ranked edgelist of the form (u, v, w) where the weights
    are predicted from the scores.
    """

    n = len(edgelist)
    scores_matrix = np.ones((n, 2))
    weights_vector = np.ones((n, 1))

    edgedict = {}
    for (p, q, w) in edgelist:
        edgedict[(p, q)] = w

    counter = 0
    for rank in range(len(ranked_edgelist)):
        (u, v, score) = ranked_edgelist[rank]

        if (u, v) in edgedict or (v, u) in edgedict:
            weight = edgedict[(u, v)] if (u, v) in edgedict else edgedict[(v, u)]
            scores_matrix[counter] = [1, score]
            weights_vector[counter] = weight
            counter += 1

    w = np.linalg.inv(scores_matrix.T @ scores_matrix) @ scores_matrix.T @ weights_vector

    new_ranked_edgelist = []
    for rank in range(len(ranked_edgelist)):
        (u, v, score) = ranked_edgelist[rank]

        input_v = np.array([1, score])
        predicted_weight = (input_v @ w)[0]
        new_ranked_edgelist.append((u, v, predicted_weight))

    return new_ranked_edgelist
    

    
