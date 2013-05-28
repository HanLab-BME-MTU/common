/*
 * See vectorialPSF.m for usage instructions.
 * Notations follow those in [F. Aguet, PhD Thesis, 2009].
 *
 * Compilation:
 * Mac/Linux: mex vectorialPSF.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /MT" vectorialPSF.cpp
 * (c) Francois Aguet 2006-2010 (last modified Nov 1, 2010)
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex>
#include "mex.h"
#include "matrix.h"

#define NARGIN 6
#define PI 3.14159265358979311599796346854

using namespace std;


// Constants for polynomial Bessel function approximation from [Abramowitz (p. 369)]
const double j0c[7] = {1, -2.2499997, 1.2656208, -0.3163866, 0.0444479, -0.0039444, 0.0002100};
const double t0c[7] = {-.78539816, -.04166397, -.00003954, 0.00262573, -.00054125, -.00029333, .00013558};
const double f0c[7] = {.79788456, -0.00000077, -.00552740, -.00009512, 0.00137237, -0.00072805, 0.00014476};
const double j1c[7] = {0.5, -0.56249985, 0.21093573, -0.03954289, 0.00443319, -0.00031761, 0.00001109};
const double f1c[7] = {0.79788456, 0.00000156, 0.01659667, 0.00017105, -0.00249511, 0.00113653, -0.00020033};
const double t1c[7] = {-2.35619449, 0.12499612, 0.00005650, -0.00637897, 0.00074348, 0.00079824, -0.00029166};

typedef struct {
    double ti0;
    double ni0;
    double ni0_2; // Ni_0^2
    double ni;
    double ni_2;
    double tg0;
    double tg;
    double ng0;
    double ng0_2;
    double ng;
    double ng_2;
    double ns;
    double ns_2;
    double lambda;
    double k0;
    double M;
    double NA;
    double NA_2;
    double alpha;
    double pixelSize;
    int sf;
    int mode;
} parameters;


// Bessel functions J0(x) and J1(x)
// Uses the polynomial approximations on p. 369-70 of Abramowitz & Stegun (1972).
// The error in J0 is supposed to be less than or equal to 5 x 10^-8.
__inline double J0(double x) {
	double r;

	if (x < 0.0) x *= -1.0;

	if (x <= 3.0){
		double y = x*x/9.0;
		r = j0c[0] + y*(j0c[1] + y*(j0c[2] + y*(j0c[3] + y*(j0c[4] + y*(j0c[5] + y*j0c[6])))));
	} else {
		double y = 3.0/x;
		double theta0 = x + t0c[0] + y*(t0c[1] + y*(t0c[2] + y*(t0c[3] + y*(t0c[4] + y*(t0c[5] + y*t0c[6])))));
		double f0 = f0c[0] + y*(f0c[1] + y*(f0c[2] + y*(f0c[3] + y*(f0c[4] + y*(f0c[5] + y*f0c[6])))));
		r = sqrt(1.0/x)*f0*cos(theta0);
	}
	return r;
}

__inline double J1(double x) {
	double r;
    double sign = 1.0;
	if (x < 0.0) {
        x *= -1.0;
        sign *= -1.0;
    }
	if (x <= 3.0){
		double y = x*x/9.0;
		r = x*(j1c[0] + y*(j1c[1] + y*(j1c[2] + y*(j1c[3] + y*(j1c[4] + y*(j1c[5] + y*j1c[6]))))));
	} else {
		double y = 3.0/x;
		double theta1 = x + t1c[0] + y*(t1c[1] + y*(t1c[2] + y*(t1c[3] + y*(t1c[4] + y*(t1c[5] + y*t1c[6])))));
		double f1 = f1c[0] + y*(f1c[1] + y*(f1c[2] + y*(f1c[3] + y*(f1c[4] + y*(f1c[5] + y*f1c[6])))));
		r = sqrt(1.0/x)*f1*cos(theta1);
	}
	return sign*r;
}


__inline void L_theta(complex<double>* L, double theta, parameters p, double ci, double z, double z_p) {
    double ni2sin2theta = p.ni_2*sin(theta)*sin(theta);
    complex<double> sroot = sqrt(complex<double>(p.ns_2 - ni2sin2theta));
    complex<double> groot = sqrt(complex<double>(p.ng_2 - ni2sin2theta));
    complex<double> g0root = sqrt(complex<double>(p.ng0_2 - ni2sin2theta));
    complex<double> i0root = sqrt(complex<double>(p.ni0_2 - ni2sin2theta));
    L[0] = p.ni*(ci - z)*cos(theta) + z_p*sroot + p.tg*groot - p.tg0*g0root - p.ti0*i0root;
    L[1] = p.ni*sin(theta) * (z-ci + p.ni*cos(theta)*(p.tg0/g0root + p.ti0/i0root - p.tg/groot - z_p/sroot));
}


// Intensity PSF for an isotropically emitting point source (average of all dipole orientations)
void psf(double *pixels, double x_p, double y_p, double z_p, double *z, int runits, int nz, parameters p) {
	
	int k;
	double r;
    int n;
        
     
    int nx = 2*runits-1;
    if (!p.mode) {
        nx *= p.sf; 
    }
    int ny = nx;
    
	// Integration parameters
	double constJ;
	int nSamples;
	double step;

	double theta, sintheta, costheta, sqrtcostheta, ni2sin2theta;
	complex<double> bessel_0, bessel_1, bessel_2, expW;
	complex<double> ngroot, nsroot;
    complex<double> ts1ts2, tp1tp2;
    complex<double> sum_I0, sum_I1, sum_I2;
    complex<double> i(0.0, 1.0);
    
    double xystep = p.pixelSize/p.M;
    
    // constant component of OPD
    double ci = z_p*(1.0 - p.ni/p.ns) + p.ni*(p.tg0/p.ng0 + p.ti0/p.ni0 - p.tg/p.ng);

    double xp_n = x_p*p.sf/xystep;
    double yp_n = y_p*p.sf/xystep;
    int rn = 1 + (int)sqrt(xp_n*xp_n + yp_n*yp_n);
    
    const int xyMax = ((2*runits-1)*p.sf-1)/2; //(nx-1)/2; // must be fine scale
    const int rMax = ceil(sqrt(2.0)*xyMax) + rn + 1; // +1 for interpolation, dx, dy
    
    
    // allocate dynamic structures
    double** integral;
    integral = new double*[nz];
    for (k=0;k<nz;k++) {
        integral[k] = new double[rMax];
    }
    
    int x, y, index, ri;
    double iconst;
    double ud = 3.0*p.sf;
    
    double w_exp;
    
    // initialize arrays, nx: coarse scale!
    int npixels = nx*ny*nz;
    for (index=0;index<npixels;index++) {
        pixels[index] = 0.0;
	}
    
    complex<double> L_th[2];
    double theta0, cst;
    
	for (k=0;k<nz;k++) {  
        
        theta0 = p.alpha;
        L_theta(L_th, theta0, p, ci, z[k], z_p);
        w_exp = abs(L_th[1]); // missing p.k0!
        
        cst = 0.975;
        while (cst >= 0.9) {
            L_theta(L_th, cst*theta0, p, ci, z[k], z_p);
            if (abs(L_th[1]) > w_exp) {
                w_exp = abs(L_th[1]);
            }
            cst -= 0.025;
        }
        w_exp *= p.k0;
        
		for (ri=0; ri<rMax; ri++) {

            r = xystep/p.sf*(double)(ri);
            constJ = p.k0*r*p.ni; // = w_J;

            if (w_exp > constJ) {
                nSamples = 4 * (int)(1.0 + p.alpha*w_exp/PI);
            } else {
                nSamples = 4 * (int)(1.0 + p.alpha*constJ/PI);
            }
            if (nSamples < 20) {
                nSamples = 20;
            }

            step =  p.alpha/(double)nSamples;
            iconst = step/ud;

            // Simpson's rule
            sum_I0 = 0.0;
            sum_I1 = 0.0;
            sum_I2 = 0.0;

            for (n=1; n<nSamples/2; n++) {
                theta = 2.0*n*step;
                sintheta = sin(theta);
                costheta = cos(theta);
                sqrtcostheta = sqrt(costheta);
                ni2sin2theta = p.ni_2*sintheta*sintheta;
                nsroot = sqrt(complex<double>(p.ns_2 - ni2sin2theta));
                ngroot = sqrt(complex<double>(p.ng_2 - ni2sin2theta));

                ts1ts2 = 4.0*p.ni*costheta*ngroot;
                tp1tp2 = ts1ts2;
                tp1tp2 /= (p.ng*costheta + p.ni/p.ng*ngroot) * (p.ns/p.ng*ngroot + p.ng/p.ns*nsroot);
                ts1ts2 /= (p.ni*costheta + ngroot) * (ngroot + nsroot);
                
                bessel_0 = 2.0*J0(constJ*sintheta)*sintheta*sqrtcostheta; // 2.0 factor : Simpson's rule
                bessel_1 = 2.0*J1(constJ*sintheta)*sintheta*sqrtcostheta;
                if (constJ != 0.0) {
                    bessel_2 = 2.0*bessel_1/(constJ*sintheta) - bessel_0;
                } else {
                    bessel_2 = 0.0;
                }
                bessel_0 *= (ts1ts2 + tp1tp2/p.ns*nsroot);
                bessel_1 *= (tp1tp2*p.ni/p.ns*sintheta);
                bessel_2 *= (ts1ts2 - tp1tp2/p.ns*nsroot);

                expW = exp(i*p.k0*((ci-z[k])*p.ni*costheta + z_p*nsroot + p.tg*ngroot - p.tg0*sqrt(complex<double>(p.ng0_2-ni2sin2theta)) - p.ti0*sqrt(complex<double>(p.ni0_2-ni2sin2theta))));
                sum_I0 += expW*bessel_0;
                sum_I1 += expW*bessel_1;
                sum_I2 += expW*bessel_2;
            }
            for (n=1; n<=nSamples/2; n++) {
                theta = (2.0*n-1)*step;
                sintheta = sin(theta);
                costheta = cos(theta);
                sqrtcostheta = sqrt(costheta);
                ni2sin2theta = p.ni_2*sintheta*sintheta;
                nsroot = sqrt(complex<double>(p.ns_2 - ni2sin2theta));
                ngroot = sqrt(complex<double>(p.ng_2 - ni2sin2theta));

                ts1ts2 = 4.0*p.ni*costheta*ngroot;
                tp1tp2 = ts1ts2;
                tp1tp2 /= (p.ng*costheta + p.ni/p.ng*ngroot) * (p.ns/p.ng*ngroot + p.ng/p.ns*nsroot);
                ts1ts2 /= (p.ni*costheta + ngroot) * (ngroot + nsroot);
                
                bessel_0 = 4.0*J0(constJ*sintheta)*sintheta*sqrtcostheta; // 4.0 factor : Simpson's rule
                bessel_1 = 4.0*J1(constJ*sintheta)*sintheta*sqrtcostheta;
                if (constJ != 0.0) {
                    bessel_2 = 2.0*bessel_1/(constJ*sintheta) - bessel_0;
                } else {
                    bessel_2 = 0.0;
                }
                bessel_0 *= (ts1ts2 + tp1tp2/p.ns*nsroot);
                bessel_1 *= (tp1tp2*p.ni/p.ns*sintheta);
                bessel_2 *= (ts1ts2 - tp1tp2/p.ns*nsroot);

                expW = exp(i*p.k0*((ci-z[k])*p.ni*costheta + z_p*nsroot + p.tg*ngroot - p.tg0*sqrt(complex<double>(p.ng0_2-ni2sin2theta)) - p.ti0*sqrt(complex<double>(p.ni0_2-ni2sin2theta))));
                sum_I0 += expW*bessel_0;
                sum_I1 += expW*bessel_1;
                sum_I2 += expW*bessel_2;
            }
            // theta = alpha;
            sintheta = sin(p.alpha);
            costheta = cos(p.alpha);
            sqrtcostheta = sqrt(costheta);
            nsroot = sqrt(complex<double>(p.ns_2 - p.NA_2));
            ngroot = sqrt(complex<double>(p.ng_2 - p.NA_2));

            ts1ts2 = 4.0*p.ni*costheta*ngroot;
            tp1tp2 = ts1ts2;
            tp1tp2 /= (p.ng*costheta + p.ni/p.ng*ngroot) * (p.ns/p.ng*ngroot + p.ng/p.ns*nsroot);
            ts1ts2 /= (p.ni*costheta + ngroot) * (ngroot + nsroot);
            
            bessel_0 = J0(constJ*sintheta)*sintheta*sqrtcostheta;
            bessel_1 = J1(constJ*sintheta)*sintheta*sqrtcostheta;
            if (constJ != 0.0) {
                bessel_2 = 2.0*bessel_1/(constJ*sintheta) - bessel_0;
            } else {
                bessel_2 = 0.0;
            }
            bessel_0 *= (ts1ts2 + tp1tp2/p.ns*nsroot);
            bessel_1 *= (tp1tp2*p.ni/p.ns*sintheta);
            bessel_2 *= (ts1ts2 - tp1tp2/p.ns*nsroot);

            expW = exp(i*p.k0*((ci-z[k])*sqrt(complex<double>(p.ni_2-p.NA_2)) + z_p*nsroot + p.tg*ngroot - p.tg0*sqrt(complex<double>(p.ng0_2-p.NA_2)) - p.ti0*sqrt(complex<double>(p.ni0_2-p.NA_2))));
            sum_I0 += expW*bessel_0;
            sum_I1 += expW*bessel_1;
            sum_I2 += expW*bessel_2;

            sum_I0 = abs(sum_I0);
            sum_I1 = abs(sum_I1);
            sum_I2 = abs(sum_I2);

            integral[k][ri] = 8.0*PI/3.0*real(sum_I0*sum_I0 + 2.0*sum_I1*sum_I1 + sum_I2*sum_I2) * iconst*iconst;
        }
    } // z loop

    // Interpolate (linear)   
    int r0;
    double dr, rx, xi, yi;
    index = 0;
    if (p.mode == 1) {
        for (k=0;k<nz;k++) {
            for (y=-xyMax; y<=xyMax; y++) {
                for (x=-xyMax; x<=xyMax; x++) {
                    xi = (double)x - xp_n;
                    yi = (double)y - yp_n;
                    rx = sqrt(xi*xi+yi*yi);
                    r0 = (int)rx;
                    if (r0+1 < rMax) {
                        dr = rx - r0;
                        index = (x+xyMax)/p.sf + ((y+xyMax)/p.sf)*nx + k*nx*ny;
                        pixels[index] += dr*integral[k][r0+1] + (1.0-dr)*integral[k][r0];
                    } // else '0'
                }
            }
        }
    } else {
        for (k=0;k<nz;k++) {
            for (y=-xyMax; y<=xyMax; y++) {
                for (x=-xyMax; x<=xyMax; x++) {
                    xi = (double)x - xp_n;
                    yi = (double)y - yp_n;
                    rx = sqrt(xi*xi+yi*yi);
                    r0 = (int)rx;
                    if (r0+1 < rMax) {
                        dr = rx - r0;
                        pixels[index] += dr*integral[k][r0+1] + (1.0-dr)*integral[k][r0];
                    } // else '0'
                    index++;
                }
            }
        }
    }
    // free dynamic structures
    for (k=0;k<nz;k++) {
        delete [] integral[k];
    }    
    delete [] integral;
} // psf




	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double x_p, y_p, z_p;
    int runits;
    
	// Check for proper number of arguments
	if (nrhs!=NARGIN) mexErrMsgTxt("There must be 6 input arguments: xp, yp, zp, z, ru, p.");
    if ( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[4])) mexErrMsgTxt("'xp', 'yp', 'zp' must be of type Double; 'ru' must be an integer.");
    if ( !mxIsDouble(prhs[3]) ) mexErrMsgTxt("Input 'z' must be a Double scalar or vector.");
    if ( !mxIsStruct(prhs[5]) ) mexErrMsgTxt("Input 'p' must be a parameter structure.");
    
    int pIndex = 5; // index of parameter structure
    double *z;
    
    z_p = mxGetScalar(prhs[2]);
    int nz = mxGetNumberOfElements(prhs[3]);
    z = mxGetPr(prhs[3]);
    runits = (int)mxGetScalar(prhs[4]);

    x_p = mxGetScalar(prhs[0]);
    y_p = mxGetScalar(prhs[1]);


    int np = mxGetNumberOfFields(prhs[pIndex]);
    parameters p;
    if (np != 13 && np != 15) mexErrMsgTxt("Incorrect parameter vector");
    for (int i=0;i<13;++i) {
        if (!mxIsDouble(mxGetFieldByNumber(prhs[pIndex],0,i))) {
            mexErrMsgTxt("All fields of the parameter structure 'p' must be Double scalars.");
        }
    }
    
    p.ti0 = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,0));
    p.ni0 = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,1));
    p.ni0_2 = p.ni0*p.ni0;
    p.ni  = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,2));
    p.ni_2 = p.ni*p.ni;
    p.tg0 = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,3));
    p.tg  = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,4));
    p.ng0 = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,5));
    p.ng0_2 = p.ng0*p.ng0;
    p.ng  = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,6));
    p.ng_2 = p.ng*p.ng;
    p.ns  = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,7));
    p.ns_2 = p.ns*p.ns;
    p.lambda    = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,8));
    p.k0 = 2*PI/p.lambda;
    p.M         = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,9));
    p.NA        = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,10));
    p.NA_2 = p.NA*p.NA;
    p.alpha     = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,11));
    p.pixelSize = mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,12));
    if (p.ni <= p.NA) mexErrMsgTxt("Refractive index must be strictly larger than NA.");
    if (np == 15) {
        p.sf = (int)mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,13));
        p.mode = (int)mxGetScalar(mxGetFieldByNumber(prhs[pIndex],0,14));
    } else {
        p.sf = 1;
        p.mode = 1;
    }
    
    
    int nx = 2*runits-1;
    if (!p.mode) {
        nx *= p.sf; 
    }
    int ny = nx;
    
    int ndim = 3;
    const int dims[3] = {ny, nx, nz}; // This order provides more efficiency due to Matlab addressing
   
    
    double *psfArray;
    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    psfArray = mxGetPr(plhs[0]);
    psf(psfArray, x_p, y_p, z_p, z, runits, nz, p);
    
}
