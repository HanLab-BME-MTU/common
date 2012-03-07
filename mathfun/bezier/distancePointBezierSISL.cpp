#include <mex.h>
#include <math.h>
#include <stdio.h>

#include "sisl.h"

// Windows: mex -I. -IC:\Users\PB93\Downloads\sisl-4.5.0\include -LC:\Users\PB93\Desktop\SISLbuild\Release -lsisl distancePointBezierSISL.cpp
// Linux: mex -I. -I/home/pb93/Downloads/sisl-4.5.0/include -L/home/pb93/Downloads/sislBuild -lsisl distancePointBezierSISL.cpp

void truncate(double *cP, int cPdim, int nCP, double tStart, double tEnd);
void transposeArray(double *a, int sizeX, int sizeY);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// =====================
	// Check inputs
	// =====================

	// Check number of input arguments
	if (nrhs != 2 && nrhs != 4) mexErrMsgTxt("2 or 4 input arguments required.");

	// Check number of output arguments
	if (nlhs > 2) mexErrMsgTxt("Too many output arguments.");

	// Check the dimension of the control points
	const mwSize *inputSize = mxGetDimensions(prhs[0]);
	int nCP = inputSize[0]; // Number of control points
	int cPdim = inputSize[1];
	if (cPdim < 1 || nCP < 1) mexErrMsgTxt("At least 1 1D control point is needed!");

	// Check the dimension of the free point
	inputSize = mxGetDimensions(prhs[1]);
	int nPoints = inputSize[0];
	int pointsDimension = inputSize[1];
	if (nPoints != 1) mexErrMsgTxt("Distance computation of only one point supported!");
	if (pointsDimension != cPdim) mexErrMsgTxt("The dimension of the point and the control points should be the same!");

	// =====================
	// Inputs/Outputs
	// =====================

	// Variables
	double *cP; // Control points
	double *point; // Point
	double *distanceMin; // Minimal distance between the point and the Bezier curve
	double *tDistMin; // Value of the parameter corresponding to the minimal distance

	// Associate inputs
	double *cPtmp = mxGetPr(prhs[0]);
	point = mxGetPr(prhs[1]);

	// Transpose control points
	cP = new double[cPdim*nCP];
	for (int i = 0; i < cPdim*nCP; i++) cP[i] = cPtmp[i];
	transposeArray(cP, nCP, cPdim);

	// Associate outputs
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	distanceMin = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	tDistMin = mxGetPr(plhs[1]);

	// Optional inputs
	if (nrhs == 4) {
		double tStart, tEnd;
		tStart = *mxGetPr(prhs[2]);
		tEnd = *mxGetPr(prhs[3]);
		truncate(cP, cPdim, nCP, tStart, tEnd);
	}

	// =====================
	// Compute the distance
	// =====================

	if (nCP == 1) { // Single control point
		double dist = 0;
		for (int i = 0; i < cPdim; i++) {
			dist = dist + pow(cP[i] - point[i],2);
		}
		*distanceMin = sqrt(dist);
		*tDistMin = 0;
	} else { // A curve

		// Convert the Bézier curve into a B-spline
		const int number = nCP;
		const int order = nCP;

		int nKnots = number+order;

		// Define the knots of the B-spline
		double *knots = new double[nKnots];
		for (int i=0;i<nKnots/2;i++) {
			knots[i] = 0;
		}
		for (int i=nKnots/2;i<nKnots;i++) {
			knots[i] = 1;
		}

		// Inputs
		SISLCurve *curve = newCurve(number, // number of control points
			order,							// order of spline curve (degree + 1)
			knots,							// pointer to knot vector (parametrization)
			cP,								// pointer to coefficient vector (control points)
			3,								// kind = Polynomial Bezier curve
			cPdim,							// dimension
			0);								// no copying of information, 'borrow' arrays
		double epsge = 1.0e-5;				// geometric tolerance
		double epsco = 0;					// Not used

		// Outputs
		int numintpt; // 
		double *intpar;
		int numintcu;
		SISLIntcurve **intcurve;
		int jstat = 0; // Status flag

		// Compute the curve parameter of the closest point (See the SISL manual)
		s1953(curve,point,cPdim,epsco,epsge,&numintpt,&intpar,&numintcu,&intcurve,&jstat);
		if (jstat < 0) mexErrMsgTxt("SISL library error in function s1953!");

		// Inputs
		int leftknot = nKnots/2 - 1;
		double *curvePoint = new double[cPdim];

		// Evaluate the curve (See the SISL manual)
		s1227(curve, 0, *intpar, &leftknot, curvePoint, &jstat);
		if (jstat < 0) mexErrMsgTxt("SISL library error in function s1227!");

		// Compute the distance to the closest point
		*tDistMin = *intpar;
		double dist = 0;
		for (int i = 0; i<cPdim; i++) {
			dist = dist + pow(curvePoint[i]-point[i],2);
		}
		
		*distanceMin = sqrt(dist);

		// Clean up
		freeCurve(curve);
		delete [] knots;
		delete [] curvePoint;
		delete [] cP;
	}
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