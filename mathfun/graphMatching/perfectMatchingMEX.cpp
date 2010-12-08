#include <mex.h>

#include <blossom5/PerfectMatching.h>

// Command line to compile on Mac OS X:
// mex -DPM_TIMER_NONE -I../../../extern/ -o perfectMacthingMEX ../../../extern/blossom5/{PMinterface,PMmain,PMexpand,PMshrink,PMinit,PMduals,misc,MinCost/MinCost}.cpp perfectMacthingMEX.cpp

// Command line to compile on Linux:
// mex -I../../../extern/ -o perfectMacthingMEX ../../../extern/blossom5/{PMinterface,PMmain,PMexpand,PMshrink,PMinit,PMduals,misc,MinCost/MinCost}.cpp perfectMacthingMEX.cpp -lrt


// plhs[0] = M             (U V)
// plhs[1] = cost
// plhs[2] = res

// prhs[0] = numVertices
// prhs[1] = edges         (U V)
// prhs[2] = weights

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
	struct PerfectMatching::Options options;
	int i, e ,node_num, edge_num;
	int* edges;
	double* weights;
	
	// Read input arguments
	
	if (nrhs != 3) {
		mexErrMsgTxt("3 input arguments required.");
	}
	
	if (nlhs > 3) {
		mexErrMsgTxt("Too many output arguments.");
	}
	
	node_num = (int) *mxGetPr(prhs[0]);
	
	int ndim = mxGetNumberOfDimensions(prhs[2]);
	
	if (ndim != 2) {
		mexErrMsgTxt("Invalid number of dimensions.");
	}
	
	const mwSize* size = mxGetDimensions(prhs[1]);
	
	edge_num = size[0];

	if (size[1] != 2) {
		mexErrMsgTxt("Invalid number of columns for the 2nd argument.");
	}

	size = mxGetDimensions(prhs[2]);
	
	if (size[0] != edge_num) {
		mexErrMsgTxt("Invalid number of elements for the 3rd argument.");
	}
		
	double* p = mxGetPr(prhs[1]);
	double* q = mxGetPr(prhs[2]);
	
	edges = new int[2*edge_num];
	weights = new double[edge_num];
	
	for (e = 0; e < edge_num; ++e) {
		edges[2 * e] = p[e] - 1;
		edges[2 * e + 1] = p[edge_num + e] - 1;
		weights[e] = q[e];
	}
	
	PerfectMatching *pm = new PerfectMatching(node_num, edge_num);

	for (e=0; e<edge_num; e++)
		pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
	
	pm->options = options;
	pm->Solve();

	// Compute the minimum cost
	double cost = 0; //ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);

	// Check the matching
	int res = 0; // CheckPerfectMatchingOptimality(node_num, edge_num, edges, weights, pm);

	// Write output
	
	plhs[0] = mxCreateDoubleMatrix(edge_num, 1, mxREAL);
	
	p = mxGetPr(plhs[0]);

	for (e = 0; e < edge_num; ++e) {
		p[e] = pm->GetSolution(e);
	}
	
	delete pm;
	
	delete [] edges;
	delete [] weights;	
}
