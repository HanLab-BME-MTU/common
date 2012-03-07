#include <mex.h>
#include <math.h>
#include <stdio.h>

#include "sisl.h"

// Compilation: mex -I. -IC:\Users\PB93\Downloads\sisl-4.5.0\include -LC:\Users\PB93\Desktop\SISLbuild\Release -lsisl truncateBezierSISL.cpp
// Linux: mex -I. -I/home/pb93/Downloads/sisl-4.5.0/include -L/home/pb93/Downloads/sislBuild -lsisl truncateBezierSISL.cpp

void truncate(double *cP, int cPdim, int nCP, double tStart, double tEnd);
void transposeArray(double *a, int sizeX, int sizeY);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// =====================
	// Check inputs
	// =====================

	// Check number of input arguments
	if (nrhs < 3) mexErrMsgTxt("Not enough input arguments.");
	if (nrhs > 4) mexErrMsgTxt("Too many input arguments.");

	// Check number of output arguments
	if (nlhs > 2) mexErrMsgTxt("Too many output arguments.");

	// Check the control points
	const mwSize *inputSize = mxGetDimensions(prhs[0]);
	int nCP = inputSize[0]; // Number of control points
	int cPdim = inputSize[1];
	if (cPdim < 1 || nCP < 2) mexErrMsgTxt("At least 2 1D control points are required!");

	// Check the parameter arguments
	inputSize = mxGetDimensions(prhs[1]);
	if (inputSize[0] != 1 || inputSize[1] != 1) mexErrMsgTxt("Second parameter must be scalar!");
	inputSize = mxGetDimensions(prhs[2]);
	if (inputSize[0] != 1 || inputSize[1] != 1) mexErrMsgTxt("Third parameter must be scalar!");

	// =====================
	// Inputs/Outputs
	// =====================

	// Variables
	double *cP, *cPout; // Control points
	double tStart, tEnd; // Start/End parameter

	double *t, *tout; // Optional parameter
	int sizeTx, sizeTy; // Size of the optional parameter

	// Associate inputs
	double *cPtmp = mxGetPr(prhs[0]);
	tStart = *mxGetPr(prhs[1]);
	tEnd = *mxGetPr(prhs[2]);

	// Transpose control points
	cP = new double[cPdim*nCP];
	for (int i = 0; i < cPdim*nCP; i++) cP[i] = cPtmp[i];
	transposeArray(cP, nCP, cPdim);

	// Associate outputs
	plhs[0] = mxCreateDoubleMatrix(nCP, cPdim, mxREAL);
	cPout = mxGetPr(plhs[0]);

	// =====================
	// Truncate the curve
	// =====================

	truncate(cP, cPdim, nCP, tStart, tEnd);
	transposeArray(cP, cPdim, nCP);

	// Copy control points
	for (int i = 0; i < cPdim*nCP; i++) cPout[i] = cP[i];

	// =====================
	// Process optional arguments
	// =====================

	if (nlhs == 2 && nrhs == 4) {

		// Associate inputs
		const mwSize *inputSize = mxGetDimensions(prhs[3]);
		sizeTx = inputSize[0];
		sizeTy = inputSize[1];
		t = mxGetPr(prhs[3]);

		// Associate outputs
		plhs[1] = mxCreateDoubleMatrix(sizeTx, sizeTy, mxREAL);
		tout = mxGetPr(plhs[1]);

		// Compute the new curve parameters
		double tDiff = tEnd-tStart;
		for (int i = 0; i < sizeTx * sizeTy; i++) {
			tout[i] = (t[i] - tStart)/tDiff;
		}
	}

	delete [] cP;
}

void truncate(double *cP, int cPdim, int nCP, double tStart, double tEnd) {
	// Convert the Bézier curve into a B-spline
	const int number = nCP;
	const int order = nCP;
	const int nKnots = number+order;

	// Define the knots of the B-spline
	double *knots = new double[nKnots];
	for (int i=0;i<nKnots/2;i++) knots[i] = 0;
	for (int i=nKnots/2;i<nKnots;i++) knots[i] = 1;

	// Inputs
	SISLCurve *curve = newCurve(number, // number of control points
		order,							// order of spline curve (degree + 1)
		knots,							// pointer to knot vector (parametrization)
		cP,								// pointer to coefficient vector (control points)
		3,								// kind = Polynomial Bezier curve
		cPdim,							// dimension
		0);								// no copying of information, 'borrow' arrays

	// Outputs
	int jstat = 0;	// Status flag
	SISLCurve *newcurve = newCurve(number,	// number of control points
		order,								// order of spline curve (degree + 1)
		knots,								// pointer to knot vector (parametrization)
		cP,									// pointer to coefficient vector (control points)
		3,									// kind = Polynomial Bezier curve
		cPdim,								// dimension
		0);									// no copying of information, 'borrow' arrays

	// Truncate the curve (See SISL manual)
	s1712(curve, tStart, tEnd, &newcurve, &jstat);
	if (jstat < 0) mexErrMsgTxt("SISL library error in function s1712!");

	// Copy the control points
	for (int i = 0; i < cPdim*nCP; i++) cP[i] = newcurve->ecoef[i];

	// Clean up
	freeCurve(curve);
	freeCurve(newcurve);
	delete [] knots;
}

void transposeArray(double *a, int sizeX, int sizeY) {
	double *aTrans = new double[sizeX * sizeY];
	for (int x = 0; x < sizeX; x++) {
		for (int y = 0; y < sizeY; y++) {
			aTrans[x * sizeY + y] = a[y * sizeX + x];
		}
	}
	for (int i = 0; i < sizeX * sizeY; i++) a[i] = aTrans[i];
	delete [] aTrans;
}