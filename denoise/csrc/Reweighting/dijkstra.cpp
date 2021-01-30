#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <map>
#include <queue>
#include <bits/stdc++.h>
#include "dijkstra.h"

/*
  @param q         <int> : Id of the neighbor node `q`
  @param weight <double> : Weight of `p` and `q`
*/
struct N_Info {
  int q;
  double weight;
  N_Info() {
  }
  N_Info(int q,
	 double weight) {
    this->q      = q;
    this->weight = weight;
    return;
  }
};

/*
  @param p             <int> : Id of the node in the network
  @param degree        <int> : degree of the node `p`
  @param neighbors <N_Info*> : Information about the neighbours of `p` 
*/
struct Neighbors {
  int p;
  int degree;
  N_Info* neighbors;
};

/*
  @param no_nodes         <int> : Number of nodes in the network
  @param no_edges         <int> : Number of edges in the network
  @param neighbors <Neighbors*> : All the nodes in the networks and their neighbors
*/
struct Graph {
  int no_nodes;
  int no_edges;
  Neighbors* neighbors;
};

/*
  Function that, given a graph `g` and a id of `source` node, returns the predicted distance
  from `source` to all the other nodes in `g`.
  @param g     <Graph *> : A graph `g`
  @param source    <int> : A source node
  @return pred <double*> : Predicted weight from source to all the other nodes in the network
*/
double* Graph_predict_source_weights_dijkstra(Graph * g,
					      int source) {
  int* prev     = new int[g->no_nodes];
  bool* traced  = new bool[g->no_nodes];
  double* pred  = new double[g->no_nodes];
  std::queue<int> to_check;

  /* Take all the node connected to the source and assign them
     a default weight value */
  for(int i = 0; i < g->no_nodes; i ++) {
    if(i == source) {
      traced[source]        = true;
      Neighbors s_neighbors = g->neighbors[source];
      for (int j = 0; j < s_neighbors.degree; j ++) {
	to_check.push(s_neighbors.neighbors[j].q);
	pred[s_neighbors.neighbors[j].q] = s_neighbors.neighbors[j].weight;
      }
      pred[source]          = 0; /* No self loop */
    }
    else {
      traced[i] = false;
    }
  }

  /* Pefrorm BFS to find the weights */
  while(to_check.size() != 0) {
    int curr                 = to_check.front();
    to_check.pop();
    Neighbors curr_neighbors = g->neighbors[curr];
    traced[curr]             = true;
    double c_weight = pred[curr];
    for(int j = 0; j < curr_neighbors.degree; j ++) {
      N_Info n_curr = curr_neighbors.neighbors[j];
      if (traced[n_curr.q] == true) {
	continue;
      }
      double n_weight = n_curr.weight;
      pred[n_curr.q] = n_weight * c_weight / (n_weight + c_weight);
      to_check.push(n_curr.q);
    }
  }
  return pred;
}


/* Function to initialize graph given adjacency matrix and the number of nodes 
   @param adj_matrix <double **> : Adjacency matrix of the nodes
   @param no_nodes         <int> : The number of nodes in the network
   @return g            <Graph*> : A graph structure representing a graph
*/
Graph* Graph_init(double ** adj_matrix,
		  int no_nodes) {
  Graph* g    = new Graph();
  g->no_nodes = no_nodes;
  g->no_edges = 0;
  g->neighbors = new Neighbors[no_nodes];
  for(int i = 0; i < no_nodes; i ++) {
    g->neighbors[i].p = i;
    int degree = 0;
    for (int j = 0; j < no_nodes; j ++) {
      if (adj_matrix[i][j] == 0) continue;
      degree ++;
    }
    g->no_edges              += degree;
    g->neighbors[i].neighbors = new N_Info[degree];
    g->neighbors[i].degree    = degree;
    for(int j = 0, d = 0; d < degree; j ++) {
      if (adj_matrix[i][j] == 0.0) continue;
      g->neighbors[i].neighbors[d].q      = j;
      g->neighbors[i].neighbors[d].weight = adj_matrix[i][j];
      d += 1;
    }
  }
  return g;
}

/*
  Function to destroy the graph
  @param g <Graph*> : Graph structure
*/
void Graph_destroy(Graph* g) {
  Neighbors * ns = g->neighbors;
  for(int i = 0; i < g->no_nodes; i ++) {
    N_Info* nns  = ns[i].neighbors;
    delete [] nns; 
  }
  delete [] ns;
  delete g;
  return;
}

/*
  Structure that stores key and value.
 */
struct Entry {
  int key, value;
  Entry() {}
  Entry(int key,
	int value) {
    this->key   = key;
    this->value = value;
    return;
  }
};

/* 
   Structure that stores an edge (p, q, wt) 
*/
struct Edge {
  int p, q;
  double wt;
  Edge() {}
  Edge(int p,
       int q) {
    this->p = p;
    this->q = q;
    this->wt = -1; /* Uninitialized */
  }
  Edge(int p,
       int q,
       double wt) {
    this->p = p;
    this->q = q;
    this->wt = wt; 
  }
  void set(int p,
	   int q) {
    this->p = p;
    this->q = q;
    this->wt = -1; /* Uninitialized */
  }
  void set(int p,
	   int q,
	   int wt) {
    this->p = p;
    this->q = q;
    this->wt = wt;
  }
};


/*
  Function for sorting the Entries
 */
bool compare_entries(Entry a,
		     Entry b) {
  return a.value > b.value; 
}


/*
 */
int * vertex_pos(Edge * edges,
		 int no_nodes,
		 int no_edges,
		 int& no_nodes_in_edges) {
  Entry * entry = new Entry[no_nodes];
  for(int i = 0; i < no_nodes; i ++) {
    entry[i].key   = i;
    entry[i].value = 0;
  }
  for(int i = 0; i < no_edges; i ++) {
    int p = edges[i].p;
    int q = edges[i].q;
    entry[p].value ++;
    entry[q].value ++;
  }
  std::sort(entry, entry + no_nodes, compare_entries);
  no_nodes_in_edges = 0;
  int* pos = new int[no_nodes];
  for(int i = 0; i < no_nodes; i ++) {
    pos[entry[i].key] = i;
    if(entry[i].value > 0) {
      no_nodes_in_edges ++;
    }
  }
  return pos;
}

typedef std::vector<std::pair<int, Edge*>> V_Edge, * V_Edge_P;
typedef std::pair<int, V_Edge_P> Pair, *Pair_P;
typedef std::vector<Pair_P> M_Edge, * M_Edge_P;

M_Edge_P rearrange_edges(Edge * edges,
			 int no_nodes,
			 int no_edges) {
  M_Edge_P edge_map = new M_Edge();
  int no_nodes_in_edges = 0;
  int * positions = vertex_pos(edges, no_nodes, no_edges, no_nodes_in_edges);
  for (int i = 0; i < no_nodes_in_edges; i ++) {
    V_Edge_P entr = new V_Edge();
    Pair_P  p_i   = new Pair(-1, entr);
    (*edge_map).push_back(p_i);
    (*edge_map)[i]->first  = -1; /* Unassigned */
  }
  for(int i = 0; i < no_edges; i ++) {
    int p   = edges[i].p;
    int q   = edges[i].q;
    int pos = positions[p];
    int maj = p;
    if (pos > positions[q]) {
      pos = positions[q];
      maj = q;
    }
    (*edge_map)[pos]->second->push_back(std::make_pair(i, edges + i));
    (*edge_map)[pos]->first = maj;
  }
  return edge_map;
}

Edge * convert_to_Edge(int ** edge_mat, int no_edges) {
  Edge * edges = new Edge[no_edges];
  for(int i = 0; i < no_edges; i ++) {
    edges[i].set(edge_mat[i][0], edge_mat[i][1]);
  }
  return edges;
}

double ** get_adj(double * flattened, int no_nodes) {
  double ** adj = new double*[no_nodes];
  for (int i = 0; i < no_nodes; i ++) {
    adj[i]      = new double[no_nodes];
    for (int j = 0; j < no_nodes; j ++) {
      adj[i][j] = flattened[i * no_nodes + j];
    }
  }
  return adj;
}

int ** get_ed(int * flattened, int no_edges) {
  int ** ed  = new int*[no_edges];
  for (int i = 0; i < no_edges; i ++) {
    ed[i]      = new int[2];
    for (int j = 0; j < 2; j ++) {
      ed[i][j] = flattened[i * 2 + j];
    }
  }
  return ed; 
}

template <class T>
void kill_mat(T ** mat, int rows) {
  for (int i = 0; i < rows; i ++) {
    delete [] mat[i];
  }
  delete [] mat;
}


/* The function to be used */
void predict_edge(void * adj_matrix_v,
		  void * edge_mat_v,
		  int no_nodes,
		  int no_edges,
		  void * weights_v) {
  double ** adj_matrix = get_adj((double *) adj_matrix_v, no_nodes);
  int **    edge_mat   = get_ed((int *) edge_mat_v, no_nodes);
  double * weights     = (double *) weights_v;

  /*
------------------DEBUG---------------------
  for(int i = 0; i < no_nodes; i ++) {
    for(int j = 0; j < no_nodes; j ++) {
      printf("%f\t ", adj_matrix[i][j]);
    }
    if (i < 3) {
      printf("\t%d  %d\n", edge_mat[i][0], edge_mat[i][1]);
    }
    else {
      printf("\n");
    }
  }
  printf("%f", adj_matrix[0][1]);
  fflush(stdout);
-------------------------------------------
  */
  
  Edge * edges      = convert_to_Edge(edge_mat, no_edges);

  M_Edge_P e_map    = rearrange_edges(edges, no_nodes, no_edges);
  Graph * g         = Graph_init(adj_matrix, no_nodes);
  int no_emap       = e_map->size();
  double ** r_edges = new double*[no_edges];
  
  for(int i = 0; i < no_emap; i ++) {
    int source     = (*e_map)[i]->first;
    double *wt     = Graph_predict_source_weights_dijkstra(g,
							   source);
    V_Edge_P e_src = (*e_map)[i]->second;
    for(int j = 0; j < e_src->size(); j ++) {
      int q_src    = (*e_src)[j].second->q;
      if(source == (*e_src)[j].second->q) {
        q_src      = (*e_src)[j].second->p;
      }
      weights[(*e_src)[j].first] = wt[q_src];
    }
  }
  Graph_destroy(g);
  kill_mat(adj_matrix, no_edges);
  kill_mat(edge_mat, no_edges);
  delete [] edges;
  return;
}

/*
------------------- DEBUG_MAIN ----------------------
int main() {
  double** adj   = new double*[5];
  adj[0] = new double[5]{0  , 1  , 0  , 0  , 0};
  adj[1] = new double[5]{1  , 0  , 0  , 0.5, 0};
  adj[2] = new double[5]{0  , 0  , 0  , 1  , 0};
  adj[3] = new double[5]{0  , 0  , 1  , 0  , 0};
  adj[4] = new double[5]{0.5, 0  , 0  , 0  , 1};

  int** edges    = new int*[3];
  edges[0]  = new int[2]{1, 2};
  edges[1]  = new int[2]{4, 2};
  edges[2]  = new int[2]{4, 3};

  double * pred_ = new double[3]{-1, -1, -1};
  predict_edge((const void *) adj, (const void *) edges, 5, 3, (void *) pred_);
  for (int i = 0; i < 3; i ++) {
    std::cout << pred_[i] << "\n";
  }
}
----------------------------------------------------
*/
