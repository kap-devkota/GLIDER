from sklearn import random_projection
from sklearn.random_projection import johnson_lindenstrauss_min_dim
import numpy as np

def r_projection(input_data, no_components = None, e = 0.1):
    if no_components == None:
        no_components = johnson_lindenstrauss_min_dim(n_samples = input_data.shape[0], eps = e)

    projected_data = random_projection.GaussianRandomProjection(n_components = no_components).fit_transform(input_data)

    return projected_data
