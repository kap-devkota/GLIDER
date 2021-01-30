
# Table of Contents

1.  [High-Level Description of Project](#org2cb932a)
2.  [Running Code Locally](#org1d093ac)
    1.  [Dependencies](#orgcbaa4f2)
    2.  [Running project locally](#org4fd94e9)
    3.  [Running project remotely on cluster](#org15c2993)


<a id="org2cb932a"></a>

# High-Level Description of Project

The high-level goal of our project is to build an algorithm to `reduce
noise` in biological networks to improve the performance of various
`information mining tasks` on the denoised network.

COPIED: We want to build an algorithm for improving the
signal-to-noise ratio of undirected, weighted networks.

Some of the `information mining tasks` considered are listed as
follows:

-   Network clustering
-   Function prediction
-   Low dimensional embedding
-   Link prediction

Evidently, the performance metrics differ between tasks. We will first
consider the task of function prediction.


<a id="org1d093ac"></a>

# Running Code Locally

Although not required to use the code, the project uses Jupyter
Notebooks to visualize, tinker, and play with the code. It is highly
recommended, but if you just want to use the code in some pipeline
then you only need to install the correct packages.


<a id="orgcbaa4f2"></a>

## Dependencies

The project has a few dependencies:

-   python 3.7+
-   matplotlib
-   scipy
-   numpy
-   networkit
-   pandas

The easiest way to set this all up is with conda, which you can
install through miniconda or anaconda. Then, run the following
commands:

    conda create --name denoise # Creates a new conda environment named denoise
    conda activate denoise # Activates the environment
    conda install jupyterlab matplotlib scipy numpy pandas # Installs packages into environment
    conda install -c conda-forge networkit 


<a id="org4fd94e9"></a>

## Running project locally

For running Jupyter notebook locally, start at the root of the directory
of the project and run the following commands:

    conda activate denoise # Activates environment with jupyter notebook
    jupyter notebook # Runs notebook locally

Everything should now work as expected!


<a id="org15c2993"></a>

## Running project remotely on cluster

Sometimes it is necessary to run the project on the cluster. To do
this, we run the notebook remotely on the cluster and use an SSH
tunnel to forward the notebook to your laptop.

[Jupyter notebook remote tutorial](https://amber-md.github.io/pytraj/latest/tutorials/remote_jupyter_notebook)

    ssh -N -f -L localhost:[LOCAL]:localhost:[REMOTE] username@your_remote_host_name
    ssh -N -f -L localhost:8888:localhost:4000 hschmi02@login.cluster.tufts.edu

