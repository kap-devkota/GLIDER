"""
TODO: Incorporate isolated vertices into 
"""

import numpy as np
from scipy.sparse import diags
def network_enhancement(W, k=20, alpha=0.9, order=2):
    """
    Performs the network enhancement algorithm.

    Inputs:
    - W the the adjacency matrix of the input network of size n x n
    - k the number of neighbors.
    - alpha is the regularization parameter.
    - order determines the extent of diffusion. Typical values are 0.5,1,2.
    Outputs:
    - An adjacency matrix of the denoised network.
    """

    # TODO: Maybe remove nodes with no neighbors? No edges to denoise.
    DD = np.sum(W, axis=0)
    W = transition_matrix(W)
    W = (W + W.T) / 2
    n, _ = W.shape
    
    P = nearest_neighbors_graph(W, min(k, n - 1))
    P = P + np.identity(n) + np.diag(np.sum(P, axis=-1))
    P = transition_fields(P)

    D, U = np.linalg.eig(P)
    D = np.diag(D)
    D = np.divide((1 - alpha) * D, 1 - alpha * np.linalg.matrix_power(D, order))

    W = U @ D @ U.T

    W_div = np.reshape(1 - np.diag(W), (n, 1))
    W_div = np.tile(W_div, (1, n))
    W = np.multiply(W, 1 - np.identity(n)) # Destroys diagonal
    W = np.divide(W, W_div) # adds it back in

    D = diags(DD)
    W = D @ W # Rescale by initial weights
    W = (W + W.T) / 2 # Convert undirected -> directed

    return W

"""
A subroutine used by the Network Enhancement algorithm.
"""
def transition_fields(W):
    n, _ = W.shape
    W = W * n # This seems to be a no-op... a little confused here.
    W = transition_matrix(W)

    w = np.sqrt(np.sum(W, axis=0))
    w = np.tile(w, (n, 1))

    W = np.divide(W, w)
    W = W @ W.T
    return W

"""
Creates the k-nearest-neighbors graph.

The undirected output graph G is defined to have the same set of
vertices and an edge (u, v) if the weight on (u, v) is one of the top
k heaviest weights coming out of u.

Inputs:
  - An adjacency matrix for a graph W
  - The number of nearest neighbors to capture knn
Outputs:
  - An adjacency matrix for the KNN graph.
"""
def nearest_neighbors_graph(W, knn):
    n, _ = W.shape

    indices = np.argsort(-1 * W)
    indicator_mat = np.zeros((n, n))

    for i in range(n):
        for k in range(knn):
            j = indices[i][k]
            indicator_mat[i, j] = 1
        
    knn_mat = np.multiply(W, indicator_mat)
    knn_mat = (knn_mat + knn_mat.T) / 2
    return knn_mat

"""
Returns the transition matrix P = D^-1W
"""
def transition_matrix(W):
    D = np.diag(1 / np.sum(W, axis=-1))
    return D@W
