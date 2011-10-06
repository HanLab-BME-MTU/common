#include <mex.h>

#include <Wm5BezierCurve3.h>
#include <Wm5Vector3.h>

using namespace Wm5;

// Get Wild Magic: from http://www.geometrictools.com/Downloads/Downloads.html

// Windows: Compile WildMagic in 64 bit RELEASE mode with Visual Studio!
// Windows: mex -I. -IC:\Users\PB93\Downloads\WildMagic5p5\GeometricTools\WildMagic5\SDK\Include -LC:\Users\PB93\Downloads\WildMagic5p5\GeometricTools\WildMagic5\SDK\Library\Release\ -lWm5Core90 -lWm5Mathematics90 lengthBezier.cpp

// Linux: 1. Remove lines in WildMagic5/makefile.wm5 to compile only LibCore and LibMathematics
// Linux: 2. Modify "CFLAGS := -c -D__LINUX__" in WildMagic5/LibCore/makeprj.wm5 WildMagic5/LibMathematics/makeprj.wm5 to "CFLAGS := -c -D__LINUX__ -fPIC"
// Linux: 3. Compile in RELEASE mode: make CFG=Release -f makefile.wm5
// Linux: mex -I. -I/home/pb93/Downloads/GeometricTools/WildMagic5/SDK/Include -L/home/pb93/Downloads/GeometricTools/WildMagic5/SDK/Library/Release -lWm5Core -lWm5Mathematics lengthBezier.cpp

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Check number of input arguments
    if (nrhs != 1 && nrhs !=3) {
        mexErrMsgTxt("1 or 3 input arguments required.");
    }
    
    // Check number of output arguments
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    // Check the dimension of the control points
    const mwSize *inputSize = mxGetDimensions(prhs[0]);
    int nControlPoints = inputSize[0];
    int controlPointDimension = inputSize[1];
    
    if (controlPointDimension != 3 || nControlPoints < 2) {
        mexErrMsgTxt("At least two 3D control points needed!");
    }
    
    // Variables
    double *cP; // Control points
    double *length; // Length of the Bezier curve
    double t0, t1; // Parametrization interval
    
    // Associate inputs
    cP = mxGetPr(prhs[0]);
    
    if (nrhs == 3) {
        t0 = *mxGetPr(prhs[1]);
        t1 = *mxGetPr(prhs[2]);
    } else {
        t0 = 0.0;
        t1 = 1.0;
    }
    
    // Associate outputs
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    length = mxGetPr(plhs[0]);
    
    // Put the control points into a vector array
    Vector3d *cPVector = new Vector3d[nControlPoints]; // Control points in a Vector3 array
    for (int i=0;i<nControlPoints;i++) {
        cPVector[i] = Vector3d(cP[0*nControlPoints+i], cP[1*nControlPoints+i], cP[2*nControlPoints+i]);
    }
    
    // Create the Bezier curve, compute its length and set it as output
    BezierCurve3d bezierCurve = BezierCurve3d(nControlPoints-1, cPVector);
    *length = abs(bezierCurve.GetLength(t0, t1));
    
    // BezierCurve3d accepts responsibility for deleting the input arrays
    
    printf("Length internal: %f\n",*length);
    
}

