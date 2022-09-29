import numpy as np
from numpy import linalg as LA


def compute_dsd_embedding(A, t = -1, gamma = 1, is_normalized = True):
    """
    Code for computing the cDSD and normalized cDSD embedding
    NOTE: This does not return the distance matrix. Additional pairwise distance computations
          need to be done in order to produce the distance matrix after the embedding is obtained
          
    parameters:
        A: The square adjacency matrix representation of a network
        t: The number of random walks for the DSE computation. Setting it to `-1` sets this to infinity
        is_normalized: If set to true, normalizes the embedding at the end
        
    output:
        A n x n output embedding. NOTE: since the dimension of X is n x n, when n is large, it is more space-efficient to use the reduced dimension representation of
        x, using the function `compute_dsd_reduced_embedding` function for `gamma` == 1
    """
    n, _ = A.shape
    e    = np.ones((n, 1))
    
    # compute the degree vector and the Transition matrix `P`
    d = A @ e
    P = A/d

    # Compute scaling `W`, which is a rank one matrix used for removing the component from P with eigenvalue 1
    W     = (1 / np.sum(d)) * np.dot(e, (e * d).T)

    # This resulting `P'` has all eigenvalues less than 1
    P1 = gamma * (P - W)
    
    # X = (I - P')^{-1}, which is the cDSD embedding
    X  = np.linalg.pinv(np.eye(n, n) - P1)

    if t > 0:
        """
        This is computationally inefficient for larger t
        """
        P1_t = np.eye(n, n) - np.linalg.matrix_power(P1, t)
        # (I - P')^{-1}(I - (P')^t) = I + P' + P'^2 + ... + P'^{t-1}
        X    = np.matmul(X, P1_t)
    
    if is_normalized == False:
        return X
    
    # If normalized is true, normalize with the square root of the steady state vector
    return (X / np.sqrt(d).T)

def compute_dsd_reduced_embedding(A, dims = 50):
    """
    Performs the dimension reduction on the normalized DSD embedding, returning a
    reduced dimension matrix.
    
    parameters:
        A: A n x n numpy matrix representing the adjacency matrix 
        dims(d): The reduced dimension 
        
    output:
        A n x d dimensional embedding of the network.
    """
    
    n, _  = A.shape
    # Get the normalized adjacency matrix N = D^{-1/2} A D^{-1/2}, L = I - N
    d1_2 = np.sqrt(A @ np.ones((n, 1)))
    N    = (A / d1_2) / d1_2.T       
    L    = np.eye(n, n) - N
    
    # Get the eigendecomposition of the laplacian
    lvals, X_ls  = LA.eig(L + np.eye(n, n))
    
    # Get the smallest `dims` eigencomponents
    l_ids = np.argsort(lvals)[1 : dims + 1]
    l_r   = lvals[l_ids] - 1
    
    X_lr  = X_ls[:, l_ids] / d1_2
    
    return np.real(X_lr) * np.real(l_r).reshape(1, -1)
