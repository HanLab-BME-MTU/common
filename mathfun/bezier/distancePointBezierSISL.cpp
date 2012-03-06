#include <mex.h>
#include <math.h>
#include <stdio.h>

#include "sisl.h"

// Compilation: mex -I. -IC:\Users\PB93\Downloads\sisl-4.5.0\include -LC:\Users\PB93\Desktop\SISLbuild\Release -lsisl distancePointBezierSISL.cpp

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// Check number of input arguments
	if (nrhs != 2) {
		mexErrMsgTxt("2 input arguments required.");
	}

	// Check number of output arguments
	if (nlhs > 2) {
		mexErrMsgTxt("Too many output arguments.");
	}

	// Check the dimension of the control points
	const mwSize *inputSize = mxGetDimensions(prhs[0]);
	int nCP = inputSize[0]; // Number of control points
	int controlPointsDimension = inputSize[1];

	if (controlPointsDimension < 1 || nCP < 1) {
		mexErrMsgTxt("At least 1 1D control point is needed!");
	}

	// Check the dimension of the free point
	inputSize = mxGetDimensions(prhs[1]);
	int nPoints = inputSize[0];
	int pointsDimension = inputSize[1];

	if (nPoints != 1) {
		mexErrMsgTxt("Distance computation of only one point supported!");
	}

	if (pointsDimension != controlPointsDimension) {
		mexErrMsgTxt("The dimension of the point and the control points should be the same!");
	}

	// Variables
	double *cP; // Control points
	double *point; // Point
	double *distanceMin; // Minimal distance between the point and the Bezier curve
	double *tDistMin; // Value of the parameter corresponding to the minimal distance

	// Associate inputs
	cP = mxGetPr(prhs[0]);
	double *cP_transposed = new double[nCP*controlPointsDimension];

	// Transpose control points
	for (int j = 0; j < nCP; j++) {
		for (int i = 0; i < controlPointsDimension; i++) {
			cP_transposed[j*controlPointsDimension+i] = cP[i*nCP+j];
		}
	}

	point = mxGetPr(prhs[1]);

	// Associate outputs
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	distanceMin = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	tDistMin = mxGetPr(plhs[1]);

	// =====================
	// Distance between two points
	// =====================

	if (nCP == 1) {
		double dist = 0;
		for (int i = 0; i<controlPointsDimension; i++) {
			dist = dist + pow(cP[i]-point[i],2);
		}
		*distanceMin = sqrt(dist);
		*tDistMin = 0;

	// =====================
	// Distance between the points and the curve
	// =====================

	} else {
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
			order,  // order of spline curve (degree + 1)
			knots,  // pointer to knot vector (parametrization)
			cP_transposed,   // pointer to coefficient vector (control points)
			1,      // kind = polynomial B-spline curve
			controlPointsDimension,      // dimension
			0);     // no copying of information, 'borrow' arrays
		double epsge = 1.0e-5; // geometric tolerance
		double epsco = 0; // Not used

		// Outputs
		int numintpt; // 
		double *intpar;
		int numintcu;
		SISLIntcurve **intcurve;
		int jstat = 0; // Status flag

		// Compute the curve parameter of the closest point (See the SISL manual)
		s1953(curve,point,controlPointsDimension,epsco,epsge,&numintpt,&intpar,&numintcu,&intcurve,&jstat);

		// Inputs
		int leftknot = nKnots/2 - 1;
		double *curvePoint = new double[controlPointsDimension];

		// Evaluate the curve (See the SISL manual)
		s1227(curve,0,*intpar,&leftknot,curvePoint,&jstat);

		// Compute the distance to the closest point
		*tDistMin = *intpar;
		double dist = 0;
		for (int i = 0; i<controlPointsDimension; i++) {
			dist = dist + pow(curvePoint[i]-point[i],2);
		}
		*distanceMin = sqrt(dist);

		// Clean up
		freeCurve(curve);
		delete [] cP_transposed;
		delete [] knots;
		delete [] curvePoint;
	}

}


