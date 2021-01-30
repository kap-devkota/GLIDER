
# Table of Contents

1.  [Project Plan :: DREAM function prediction pipeline <code>[75%]</code>](#orgd073853)
    1.  [Parse DREAM networks](#orgcbb03ce)
    2.  [Parse GO labels <code>[33%]</code>](#orgdb4d83f)
        1.  [Decide on which labels to keep](#org4c97453)
        2.  [Load raw labels from FA file](#orgaba3c54)
    3.  [Implement network denoising algorithm(s) <code>[100%]</code>](#orgff0cdee)
        1.  [Implement DSD denoising algorithm](#org0632f35)
        2.  [Implement identity denoising algorithm](#orgc1576ac)
        3.  [Implement NE denoising algorithm <code>[66%]</code>](#org5192c54)
    4.  [Implement scoring strategy](#orgf1e45dc)
        1.  [k-fold cross validation algorithm](#org1b41fbd)
    5.  [Create a program that takes in a graph, embedding -> functional labelling](#orga3fff57)
        1.  [: Created a program that takes in graph, produces embedding and returns the similarity matrix. Example in \`DREAM Network Testbench\`](#orge2fa776)
    6.  [(Henri) Added the embedding and similarity matrix producing code in \`DREAM Network Testbench\` notebook. Need to compute a    function that takes in similarity matrix and produces a function label score.](#org6757c3d)
    7.  [(Henri) Test function prediction against the network enhancement](#org7eccfed)
    8.  [(Henri/Kapil) Set up SLURM repository for collaboration](#orgcc9c8d7)
    9.  [(Henri) Try different voting methods (KNN, WMV)](#org7fec321)
    10. [(Henri) Try other kernels besides RBF kernel](#org8bd9433)
    11. [(Kapil) Create small datasets for us to work on.](#orgcbd7af0)
    12. [(Kapil) Write up GLIDE code.](#org45b311c)


<a id="orgd073853"></a>

# TODO Project Plan :: DREAM function prediction pipeline <code>[75%]</code>


<a id="orgcbb03ce"></a>

## DONE Parse DREAM networks

-   DREAM networks consist of adjacency lists with proteins listed as
    nodes.
-   Need to decide how I'm going to load the networks for processing.
-   I should keep the node names with their indices so it is easy to
    attach labels for scoring.
-   Ensure everything is **deanonymized**. Keep it simple.


<a id="orgdb4d83f"></a>

## DONE Parse GO labels <code>[33%]</code>

We need to load the labels for the nodes. There should be a standard
way to keep the labels with the nodes to ensure that scoring is easy.
which labels are kept is another problem of concern. Currently,
I use all labels for the nodes and the labels are stored in the format

GO<sub>LABEL1</sub> PROTEIN1 PROTEIN2 PROTEIN3 ..
GO<sub>LABEL2</sub> PROTEIN1 PROTEIN2 PROTEIN3 ..


<a id="org4c97453"></a>

### TODO Decide on which labels to keep


<a id="orgaba3c54"></a>

### DONE Load raw labels from FA file


<a id="orgff0cdee"></a>

## DONE Implement network denoising algorithm(s) <code>[100%]</code>


<a id="org0632f35"></a>

### DONE Implement DSD denoising algorithm


<a id="orgc1576ac"></a>

### DONE Implement identity denoising algorithm


<a id="org5192c54"></a>

### DONE Implement NE denoising algorithm <code>[66%]</code>

1.  DONE Reproduce scoring on raw butterfly network.

2.  DONE Port over denoising algorithm.

3.  TODO Implement zero handling strategy.


<a id="orgf1e45dc"></a>

## DONE Implement scoring strategy


<a id="org1b41fbd"></a>

### DONE k-fold cross validation algorithm


<a id="orga3fff57"></a>

## DONE Create a program that takes in a graph, embedding -> functional labelling

Does this based on the l2 distance induced by the embedding. DSD
embedding. [Embedding Code](https://github.com/kap-devkota/Trimming_Functional/blob/master/src/Utils/dse_computations.py): DSD embedding compute<sub>X</sub><sub>normalized</sub>(&#x2026;).


<a id="orge2fa776"></a>

### DONE : Created a program that takes in graph, produces embedding and returns the similarity matrix. Example in \`DREAM Network Testbench\`


<a id="org6757c3d"></a>

## DONE (Henri) Added the embedding and similarity matrix producing code in \`DREAM Network Testbench\` notebook. Need to compute a    function that takes in similarity matrix and produces a function label score.


<a id="org7eccfed"></a>

## TODO (Henri) Test function prediction against the network enhancement


<a id="orgcc9c8d7"></a>

## TODO (Henri/Kapil) Set up SLURM repository for collaboration

What I would like is a repository somewhere on the cluster that is
extremely organized for all of our slurm scripts/python scripts that
test various parameters. It should **not** have any core code (like your
GLIDE method).

There should be top-level module that imports all of our code and then
we can use that in all of our one-off python scripts. Then, there
should be a place for the SLURM code that calls the python scripts.

For the data directory, it would be useful to have some sort of
organization so things don't become too chaotic. Maybe a README
in each directory describing what everything is. Also for
very temporary stuff we can use a /temp directory


<a id="org7fec321"></a>

## DONE (Henri) Try different voting methods (KNN, WMV)


<a id="org8bd9433"></a>

## TODO (Henri) Try other kernels besides RBF kernel


<a id="orgcbd7af0"></a>

## DONE (Kapil) Create small datasets for us to work on.


<a id="org45b311c"></a>

## DONE (Kapil) Write up GLIDE code.

