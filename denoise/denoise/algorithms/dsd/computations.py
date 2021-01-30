import numpy as np
from denoise.graph.operations import densify
from numpy import linalg as LA

def compute_degree_mat(A):
    e = np.ones((A.shape[0], 1))
    deg = np.matmul(A, e)
    return np.diag(deg.flatten())

def compute_laplacian(A):
    D = np.diag(deg.flatten())
    return D - A

def compute_pinverse_diagonal(diag):
    i_diag = np.zeros((diag.shape[0], diag.shape[0]))
    for i in range(diag.shape[0]):
        di = diag[i, i]
        if di != 0.0:
            i_diag[i, i] = 1 / float(di)

    return i_diag

def compute_X_normalized(A, D, t = -1, lm = 1, is_normalized = True):
    D_i = compute_pinverse_diagonal(D)
    P = np.matmul(D_i, A)
    Identity = np.identity(A.shape[0])
    e = np.ones((A.shape[0], 1))

    # Compute W
    scale = np.matmul(e.T, np.matmul(D, e))[0, 0]
    W = np.multiply(1 / scale, np.matmul(e, np.matmul(e.T, D)))

    up_P = np.multiply(lm, P - W)
    X_ = Identity - up_P
    X_i = np.linalg.pinv(X_)

    if t > 0:
        LP_t = Identity - np.linalg.matrix_power(up_P, t)
        X_i = np.matmul(X_i, LP_t)
    
    if is_normalized == False:
        return X_i
    
    # Normalize with steady state
    SS = np.sqrt(np.matmul(D, e))
    SS = compute_pinverse_diagonal(np.diag(SS.flatten()))

    return np.matmul(X_i, SS)

def compute_sim_matrix(embedding_mat, meth="rbf", params = None):
    """
    Given a embedding matrix, returns the similarity matrix
    @param embedding_mat -> A  (nxn) numpy matrix containing embedding
    @param meth          -> A  method by which the similarity metric computed
                            Currently only rbf kernel implemented
    @param params        -> parameters required for @meth
    @return sim          -> Similarity matrix computed from the embedding
    """
    prod    = embedding_mat @ embedding_mat.T
    diag    = np.diag(prod).reshape(-1, 1)
    e       = np.ones(diag.shape)
    Diag    = diag @ e.T + e @ diag.T
    l2_     = Diag - 2 * prod
    sim     = None
    if meth == "rbf":
        sim = np.exp(-1 / (2 * params["rbf"] ** 2) * l2_)
    return sim

# TODO: This function feels out of place. It doesn't do anything
# interesting and clutters the codebase.
def compute_embedding(edge_list, lm = 1, is_normalized = True):
    A     = densify(edge_list)
    D     = compute_degree_mat(A)
    X     = compute_X_normalized(A, D, lm = lm, is_normalized = is_normalized)
    return X


def compute_reduced_X_embedding(A, D, dims = 50):
    n            = A.shape[0]
    D_p          = np.sqrt(D)
    D_n          = np.sqrt(compute_pinverse_diagonal(D))
    N            = np.matmul(np.matmul(D_n, A), D_n)
    L            = np.identity(n) - N
    spec, X_spec = LA.eig(L + np.identity(n))
    index_chosen = np.argsort(spec)[1 : dims + 1]
    spec_mat     = spec[index_chosen] - 1
    spec_mat     = np.diag(spec_mat)
    X_spec_r     = X_spec[:, index_chosen]
    X_spec_r     = np.matmul(D_n, X_spec_r)
    reduced_mat  = np.matmul(X_spec_r, spec_mat)
    return reduced_mat
