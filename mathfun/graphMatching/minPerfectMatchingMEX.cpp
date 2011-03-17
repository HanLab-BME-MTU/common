#include <iostream>

#include <mex.h>

#include <blossom5/PerfectMatching.h>

// Command line to compile on Mac OS X:
// mex -DPM_TIMER_NONE -I../../../extern/ -o perfectMatchingMEX ../../../extern/blossom5/{PMinterface,PMmain,PMexpand,PMshrink,PMinit,PMduals,misc,MinCost/MinCost}.cpp perfectMatchingMEX.cpp

// Command line to compile on Linux:
// mex -I../../../extern/ -o perfectMatchingMEX ../../../extern/blossom5/{PMinterface,PMmain,PMexpand,PMshrink,PMinit,PMduals,misc,MinCost/MinCost}.cpp perfectMatchingMEX.cpp -lrt


// plhs[0] = M

// prhs[0] = numVertices
// prhs[1] = edges
// prhs[2] = weights

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  // Read input arguments
	
  if (nrhs != 3) {
    mexErrMsgTxt("3 input arguments required.");
  }
	
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }
	
  int node_num = (int) *mxGetPr(prhs[0]);
	
  int ndim = mxGetNumberOfDimensions(prhs[2]);
	
  if (ndim != 2) {
    mexErrMsgTxt("Invalid number of dimensions.");
  }
	
  const mwSize* size = mxGetDimensions(prhs[1]);
	
  int edge_num = size[0];

  if (size[1] != 2) {
    mexErrMsgTxt("Invalid number of columns for the 2nd argument.");
  }

  size = mxGetDimensions(prhs[2]);
	
  if (size[0] != edge_num) {
    mexErrMsgTxt("Invalid number of elements for the 3rd argument.");
  }
	
  double* p = mxGetPr(prhs[1]);
  void* q = mxGetPr(prhs[2]);
	
  PerfectMatching *pm = new PerfectMatching(node_num, edge_num);

  for (int e = 0; e < edge_num; ++e) {
    pm->AddEdge((int) p[e] - 1,
		(int) p[edge_num + e] - 1,
		((int*) q)[e]);
  }
	
  pm->Solve();

  // Write output
	
  plhs[0] = mxCreateDoubleMatrix(edge_num, 1, mxREAL);
	
  p = mxGetPr(plhs[0]);

  for (int e = 0; e < edge_num; ++e) {
    p[e] = pm->GetSolution(e);
  }
	
  delete pm;
}
