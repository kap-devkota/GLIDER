# GLIDETOOLS: A python based package for computing Diffusion State Distance and GLIDE

## Licensing

This code is copyrighted under the MIT License. 

## Publications

**GLIDE:** Devkota, Kapil, James M. Murphy, and Lenore J. Cowen. "GLIDE: combining local methods and diffusion state embeddings to predict missing interactions in biological networks." Bioinformatics 36.Supplement_1 (2020): i464-i473.

**GLIDER:** Devkota, K., Schmidt, H., Werenski, M., Murphy, J.M., Erden, M., Arsenescu, V. and Cowen, L.J., 2022. GLIDER: Function Prediction from GLIDE-based Neigborhoods. Bioinformatics.

**DSD:** Cao, Mengfei, et al. "New directions for diffusion-based network prediction of protein function: incorporating pathways with confidence." Bioinformatics 30.12 (2014): i219-i227. 

## Dependencies

This package requires the following dependencies:

1. numpy
2. scipy
3. pandas
4. matplotlib
5. json
6. networkx

## How to install

You can install the package using `pip`. 

```
pip install glidetools
```

You can also go to the glidetools repository at `https://github.com/kap-devkota/GLIDER` and clone the latest version under the `main branch`.
After you enter the `glidetools` folder, run 

```
python -m pip install glidetools
```

## Package Description

### Computing DSD matrix

This can be done by using the function  `glidetools.algorithm.dsd:compute_dsd_embedding` 

```
def compute_dsd_embedding(A, 
                        t = -1, 
                        gamma = 1, 
                        is_normalized = True)
```

Where,
- A: a numpy adjacency matrix (N x N)
- t: The number of random walks to get the DSD matrix. Setting `t` to a negative value implies `t` is infinity.
- gamma: Set it to 1 to get the default cDSD embedding
- is_normalized: If set to True, a normalized form of cDSD (by the steady state vector) is returned

This function returns a (N x N) cDSD embedding. **Note:** The output is an embedding, not a distance. To compute the cDSD distance,
do the following

```
from scipy.spatial.distance import squareform, pdist
squareform(pdist(X))
```

Where, `X` is the output from the `compute_dsd_embedding` function.


### Computing the GLIDE Matrix

This can be done by using the function  `glidetools.algorithm.glide:glide`

```
def glide(A, 
          alpha = 0.1,
          beta  = 1000,
          delta = 1,
          gamma = 1,
          normalize_dsd = False,
          local = "",
          **kwargs)
```

Where,
- A : A N x N numpy matrix
- alpha, beta, delta, gamma => glide parameters: For more information, see :
- normalize_dsd: If set to True, generates the normalized version of DSD embedding
- local: Can be either `cw`(common weighted) or [`l3`](https://www.nature.com/articles/s41467-019-09177-y).  
    
    
You can also provide your own local and global functions for GLIDE
    
- localf: a custom function that takes in adjacency matrix and returns the local pairwise score
- globalf: a custom function that takes in adjacency matrix and returns the global pairwise score

### Using `glide_compute`

If you have installed the pip package, you can the entrypoint `glide_compute` to obtain both the DSD and GLIDE outputs. 

```
usage: glide-compute [-h] [--network NETWORK] [--output OUTPUT] [-v] [--return-dsd-emb] [--return-dsd-dist] [--dsd-dist-norm {l1,l2}] [--normalized] [--reduced-dims REDUCED_DIMS] [--gamma GAMMA] [--get-glide-neighbors]
                     [--glide-neighbors-k GLIDE_NEIGHBORS_K] [--neighbors-return-format {dataframe,graph,map}] [--alpha ALPHA] [--beta BETA] [--delta DELTA] [--local {cw,l3}] [--normalize-local] [--weighted-local] [--scale-local]

optional arguments:
  -h, --help            show this help message and exit
  --network NETWORK     A Tab-delimited network file
  --output OUTPUT       The output URL. If the output is a matrix, it is always saved in a pickle format along with the name-to-index mapping dictionary
  -v                    Verbose mode
  --return-dsd-emb      If set to True, only returns the DSD embedding, else returns the GLIDE matrix
  --return-dsd-dist     If set to True, bypasses the --return-dsd-emb command to return the pairwise distance matrix from the dsd embedding
  --dsd-dist-norm {l1,l2}
                        Only used in conjunction with the --return-dsd-dist argument. Decides whether to use the `l1` or `l2` norm while computing distance
  --normalized          If set to false, returns the classic cDSD, else returns normalized cDSD embedding.
  --reduced-dims REDUCED_DIMS
                        If set to a positive value, the output is a reduced normalized DSD with reduced dimensions given by --reduced_dims
  --gamma GAMMA         DSD gamma parameter
  --get-glide-neighbors
                        If set to true, --get_glide_neighbors returns glide neighbors instead of glide matrix
  --glide-neighbors-k GLIDE_NEIGHBORS_K
                        If --get_glide_neighbors is set to true, the code uses --glide_neighbors to decide on the number of neighbors
  --neighbors-return-format {dataframe,graph,map}
                        This parameter decides the output format for the GLIDE neighbors. If `dataframe` is selected, the code returns output as a panda DataFrame.If `graph` is selected, the code returns output as a networkx graph,
                        otherwise the output is returned as a simple dictionary {NODE: LIST[NODE]}, where LIST[NODE]is the list of neighbors for the particular node
  --alpha ALPHA         GLIDE alpha parameter
  --beta BETA           GLIDE beta parameter
  --delta DELTA         GLIDE delta parameter
  --local {cw,l3}       The local parameter for GLIDE
  --normalize-local     If set to False, the local measures are not normalized
  --weighted-local      If set to False, the adjacency matrix is converted to a unweighted form (setting every non-zero elements to 1)before applying local measures
  --scale-local         If set to True, scales the local measures by their max value before returning it
```



