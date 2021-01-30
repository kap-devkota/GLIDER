'''
Reference implementation of node2vec. 
Author: Aditya Grover
For more details, refer to the paper:
node2vec: Scalable Feature Learning for Networks
Aditya Grover and Jure Leskovec 
Knowledge Discovery and Data Mining (KDD), 2016
'''

import argparse
import numpy as np
import networkx as nx
import denoise.algorithms.n2vec.node2vec as nv
from gensim.models import Word2Vec
from datetime import datetime
import os

def default_args(args = {}):
	keys   = ["input", "dimensions", "walk-length", "num-walks", "window-size", "iter", "workers", "p", "q", "weighted", "unweighted", "directed", "undirected"]
	values = ["karate.edgelist",
		  50,
		  80,
		  10,
		  10,
		  1,
		  8,
		  1,
		  1,
		  True,
		  False,
		  False,
		  True]
	for i in range(len(keys)):
		key   = keys[i]
		value = values[i]
		if key not in args:
			args[key] = value
	return args
		
def read_graph(args):
	'''
	Reads the input network in networkx.
	'''
	if args["weighted"]:
		G = nx.read_edgelist(args["input"], nodetype=int, data=(('weight',float),), create_using=nx.DiGraph())
	else:
		G = nx.read_edgelist(args["input"], nodetype=int, create_using=nx.DiGraph())
		for edge in G.edges():
			G[edge[0]][edge[1]]['weight'] = 1

	if not args["directed"]:
		G = G.to_undirected()

	return G

def learn_embeddings(walks, args):
	'''
	Learn embeddings by optimizing the Skipgram objective using SGD.
	args -- A dictionary
	    dimensions - 
	'''
	walks = [list(map(str, walk)) for walk in walks]
	print(args)
	model = Word2Vec(walks, size=args["dimensions"], window=args["window-size"], min_count=0, sg=1, workers=args["workers"], iter=args["iter"])
	print("Here")
	fname = datetime.now().strftime("%d-%m-%Y-%H-%M-%S.txt")
	model.wv.save_word2vec_format(fname)
	with open(fname, "r") as fp:
		emb = []
		id  = 0
		node_map = {}
		for line in fp:
			words                   = line.rstrip().lstrip().split()
			node_map[int(words[0])] = id
			id                     += 1
			vec                     = [float(f) for f in words[1: ]]
			emb.append(vec)
	s_list = []
	for i in range(len(node_map)):
		s_list.append(node_map[i])
	emb = np.array(emb)
	return emb[s_list]


def compute_embedding(edgelist, args):
	'''
	Pipeline for representational learning for all nodes in a graph.
	'''
	args = default_args(args)
	nx_G = nx.Graph()
	nx_G.add_weighted_edges_from(edgelist)
	G = nv.Graph(nx_G, args["directed"], args["p"], args["q"])
	G.preprocess_transition_probs()
	walks = G.simulate_walks(args["num-walks"], args["walk-length"])
	return learn_embeddings(walks, args)
