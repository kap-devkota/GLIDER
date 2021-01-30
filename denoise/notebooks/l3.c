#include <math.h>

void compute_l3_score(const void * p_neighbors, const void * q_neighbors, int size_pn, int size_qn, const void * edge_matrix, const void * n_degrees, int no_nodes, const void * out) {
  const int * p_neigh   = (int*) p_neighbors;
  const int * q_neigh   = (int*) q_neighbors;
  const double * edge_m = (double*) edge_matrix;
  const int * n_deg     = (int*) n_degrees;
  const double * op     = (double*) out; 
  double score = 0;
  for (int i = 0 ; i < size_pn; i += 1) {
    for (int j = 0 ; j < size_qn; j += 1) {
      int e1_id      = p_neigh[i];
      int e2_id      = q_neigh[j];
      double a_e1_e2 = edge_m[e1_id * no_nodes + e2_id];
      if (a_e1_e2 > 0) a_e1_e2 = 1;
      score         += a_e1_e2 / sqrt(n_deg[e1_id] * n_deg[e2_id]);
    }
  }
  printf()
  op[0] = score;
}
