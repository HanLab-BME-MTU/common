#include <mex.h>
#include <math.h>
#include <stdio.h>

#include "sisl.h"

// Compilation: mex -I. -IC:\Users\PB93\Downloads\sisl-4.5.0\include -LC:\Users\PB93\Desktop\SISLbuild\Release -lsisl lengthBezierSISL.cpp

void truncate(double *cP, int cPdim, int nCP, double tStart, double tEnd);
void transposeArray(double *a, int sizeX, int sizeY);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// =====================
	// Check inputs
	// =====================

	// Check the number of output arguments
	if (nlhs > 1) mexErrMsgTxt("Only one output argument supported!");

	// Check the number of inputs
	if (nrhs != 1 && nrhs != 3) mexErrMsgTxt("Only 1 or 3 input argument(s) supported!");

	// Check the control points
	const mwSize *inputSize = mxGetDimensions(prhs[0]);
	int nCP = inputSize[0]; // Number of control points
	int cPdim = inputSize[1]; // Control points dimension
	if (cPdim < 1 || nCP < 1) mexErrMsgTxt("At least 1 1D control point is required!");

	// =====================
	// Inputs/Outputs
	// =====================

	// Variables
	double *cP; // Control points
	double *length; // Output

	double tStart, tEnd; // Optional start/end parameter

	// Associate inputs
	cP = mxGetPr(prhs[0]);

	// Transpose control points
	transposeArray(cP, cPdim, nCP);

	// Associate outputs
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	length = mxGetPr(plhs[0]);

	if (nCP == 1) *length = 0; return; 

	// =====================
	// Process optional arguments
	// =====================

	if (nrhs == 3) {

		// =====================
		// Check optinal inputs
		// =====================

		inputSize = mxGetDimensions(prhs[1]);
		if (inputSize[0] != 1 || inputSize[1] != 1) mexErrMsgTxt("Second parameter must be scalar!");

		inputSize = mxGetDimensions(prhs[2]);
		if (inputSize[0] != 1 || inputSize[1] != 1) mexErrMsgTxt("Third parameter must be scalar!");

		// =====================
		// Optional inputs/outputs
		// =====================

		double tStart, tEnd; // Start/End parameter
		tStart = *mxGetPr(prhs[1]);
		tEnd = *mxGetPr(prhs[2]);

		// =====================
		// Truncate the curve
		// =====================

		truncate(cP, cPdim, nCP, tStart, tEnd);
	}

	// =====================
	// Compute the length
	// =====================

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
	double epsge = 1.0e-5; // Geometric precision

	// Outputs
	int jstat = 0;	// Status flag

	// Compute the length of the curve (See SISL manual)
	s1240(curve, epsge, length, &jstat);
	if (jstat < 0) mexErrMsgTxt("SISL library error in function s1240!");

	freeCurve(curve);
	delete [] knots;

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