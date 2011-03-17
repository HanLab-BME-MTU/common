#include <iostream>

#include <mex.h>

#include <lemon/list_graph.h>

// Command line to compile on Mac OS X:
// mex maxWeightedMatching.cpp

// Command line to compile on Linux:
// mex maxWeightedMatching.cpp

using namespace lemon;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ListDigraph g;
	
	ListDigraph::Node u = g.addNode();
	ListDigraph::Node v = g.addNode();
	ListDigraph::Arc  a = g.addArc(u, v);
	
	cout << "Hello World! This is LEMON library here." << endl;
	cout << "We have a directed graph with " << countNodes(g) << " nodes "
	<< "and " << countArcs(g) << " arc." << endl;
}
