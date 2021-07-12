#include <math.h>
#include <stdio.h>

void compute_l3_score(const void * p_neighbors,
		      const void * q_neighbors,
		      int size_pn,
		      int size_qn,
		      const void * edge_matrix,
		      const void * n_degrees,
		      int no_nodes,
		      void * out) {
  const int * p_neigh   = (int*) p_neighbors;
  const int * q_neigh   = (int*) q_neighbors;
  const double * edge_m = (double*) edge_matrix;
  const int * n_deg     = (int*) n_degrees;
  double * op           = (double *) out;
  double score = 0;
  for (int i = 0 ; i < size_pn; i += 1) {
    for (int j = 0 ; j < size_qn; j += 1) {
      int e1_id      = p_neigh[i];
      int e2_id      = q_neigh[j];
      double a_e1_e2 = edge_m[e1_id * no_nodes + e2_id];
      if (a_e1_e2 > 0)
      	score += 1 / sqrt(n_deg[e1_id] * n_deg[e2_id]);
    }
  }
  op[0] = score;
  return;
}

void compute_cw_score(int p,
		      int q,
		      const void * p_neighbors,
		      const void * q_neighbors,
		      int size_pn,
		      int size_qn,
		      const void * edge_matrix,
		      int no_nodes,
		      void * out) {
  const int * p_neigh   = (int*) p_neighbors;
  const int * q_neigh   = (int*) q_neighbors;
  const double * edge_m = (double*) edge_matrix;
  double * op           = (double *) out;
  double score          = 0;
  /* swap p_neigh and q_neigh if size_pn > size_qn */
  if (size_pn > size_qn) {
    for (int i = 0; i < size_qn; i += 1) {
      /* (p, q_neigh[i]) should have an edge,
	 if p and q has a neighbor in common */
      if (edge_m[p * no_nodes + q_neigh[i]] > 0) {
	score += (edge_m[q * no_nodes + q_neigh[i]]
		  + edge_m[p * no_nodes + q_neigh[i]]);
      }
    }
  }
  else {
    for (int i = 0; i < size_pn; i += 1) {
      /* (q, p_neigh[i]) should have an edge,
	 if p and q has a neighbor in common */
      if (edge_m[q * no_nodes + p_neigh[i]] > 0) {
	score += (edge_m[q * no_nodes + p_neigh[i]]
		  + edge_m[p * no_nodes + p_neigh[i]]);
      }
    }
  } 
  op[0] = score;
}
