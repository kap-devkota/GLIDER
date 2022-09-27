import numpy as np

def l3(A, **kwargs):
    """
        Implementation of: 
                Kovacs, Istvan et al. Network-based prediction of protein interactions. Nature Communications. 10. 10.1038/s41467-019-09177-y.  

        Given an adjacency matrix `A` matrix representing the L3 score between all pairwise nodes
        parameters:
            A - n x n matrix representing the adjacency matrix
            **kwargs - will check if it contains `weighted` key, if set to `True`, return the weighted version of L3
    """
    assert isinstance(A, np.ndarray) and (len(A.shape) == 2) and (A.shape[0] == A.shape[1])
    
    weighted = False if ("weighted" not in kwargs) else kwargs["weighted"]
    n, _     = A.shape
    if not weighted:
        A = np.where(A > 0, 1, 0)
    
    d1_2    = np.sqrt(A @ np.ones((n, 1)))
    
    # a_{p, q} = \sum_{u, v} \frac{a_{p, u} a_{u,v} a_{v, q}}{\sqrt{k_u * k_v}} 
    A_   =  (A / d1_2) / d1_2.T
    
    R = A @ A_ @ A
    
    if "rescale" in kwargs:
        R = R / np.max(R)
    
    return R

def cw(A, **kwargs):
    """
    Given an adjacency matrix `A` matrix representing the cw score between all pairwise nodes
        parameters:
            A - n x n matrix representing the adjacency matrix
            **kwargs - will check if it contains `weighted` key, if set to `True`, return the weighted version of cw
    """
    assert isinstance(A, np.ndarray) and (len(A.shape) == 2) and (A.shape[0] == A.shape[1])
    
    weighted   = False if ("weighted" not in kwargs) else kwargs["weighted"]
    normalized = True if ("normalize" not in kwargs) else kwargs["normalize"]
    n, _ = A.shape
    if not weighted: 
        A  = np.where(A > 0, 1, 0)
    A1 = np.where(A > 0, 1, 0)
    
    # CW(p, q) = \sum_{r} a_{p, r} + a_{q, r}
    C  = A @ A1
    C  = C + C.T
    
    if normalized:
        d1_2 = np.sqrt(A @ np.ones((n, 1)))
        C    =(C / d1_2) / d1_2.T
    
    if "rescale" in kwargs:
        C = C / np.max(C)

    return C
    