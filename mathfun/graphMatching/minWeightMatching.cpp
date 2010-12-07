#include <mex.h>

#include <blossom5/PerfectMatching.h>

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  struct PerfectMatching::Options options;
  int i, e ,node_num, edge_num;
  int* edges;
  int* weights;

  PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
  for (e=0; e<edge_num; e++) pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
  pm->options = options;
  pm->Solve();

  // Check the matching
  int res = CheckPerfectMatchingOptimality(node_num, edge_num, edges, weights, pm);
  printf("check optimality: res=%d (%s)\n", res, (res==0) ? "ok" : ((res==1) ? "error" : "fatal error"));

  double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
  printf("cost = %.1f\n", cost);

  delete pm;  
}
