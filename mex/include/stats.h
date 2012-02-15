/* Collection of functions for computing statistical tests
 *
 * (c) Francois Aguet, 2011 (last modified Mar 20, 2011)
 * */

#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
    #define fmax max
    #define fmin min
    #define erf gsl_sf_erf
    #include <minmax.h>
    #include <gsl/gsl_sf_erf.h>
#endif



/* Code from:
 * Marsaglia, G., W.W. Tsang, and J. Wang (2003), "Evaluating Kolmogorov's
 * Distribution", Journal of Statistical Software, vol. 8, issue 18.
 * This is the implementation used in kstest.m in Matlab.
 */
void mMultiply(double *A, double *B, double *C, int m)
{ 
    int i,j,k;
    double s;
    for(i=0;i<m;i++)
        for(j=0; j<m; j++) {
        s=0.0; 
        for(k=0;k<m;k++) 
            s+=A[i*m+k]*B[k*m+j];
        C[i*m+j]=s;
        }
}

void mPower(double *A, int eA, double *V, int *eV, int m, int n)
{ 
    double *B;
    int eB,i;
    if(n==1) {
        for (i=0;i<m*m;i++)
            V[i]=A[i];
        *eV=eA; 
        return;
    } 
    mPower(A, eA, V, eV, m, n/2);
    B = (double*)malloc((m*m)*sizeof(double));
    mMultiply(V, V, B, m);
    eB = 2*(*eV);
    if (n%2==0) {
        for (i=0;i<m*m;i++)
            V[i]=B[i];
        *eV=eB;
    } else {
        mMultiply(A,B,V,m);
        *eV=eA+eB;
    } 
    if (V[(m/2)*m+(m/2)]>1e140) {
        for(i=0;i<m*m;i++)
            V[i] = V[i]*1e-140;
        *eV+=140;
    }
    free(B);
}

double K(int n, double d) {
    int k, m, i, j, g, eH, eQ;
    double h, s, *H, *Q;
    /* OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL */
    s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt((double)n)+1.409/n)*s);
    k = (int)(n*d)+1;
    m=2*k-1;
    h=k-n*d;
    H = (double*)malloc((m*m)*sizeof(double));
    Q = (double*)malloc((m*m)*sizeof(double));
    for (i=0;i<m;i++)
        for(j=0;j<m;j++)
            if(i-j+1<0)
                H[i*m+j]=0;
            else
                H[i*m+j]=1;
    for(i=0;i<m;i++) {
        H[i*m] -= pow(h, i+1);
        H[(m-1)*m+i] -= pow(h, (m-i));
    }
    H[(m-1)*m] += (2*h-1>0 ? pow(2*h-1, m) : 0);
    for (i=0;i<m;i++)
        for (j=0;j<m;j++)
            if (i-j+1>0)
                for (g=1;g<=i-j+1;g++)
                    H[i*m+j]/=g;
    eH=0;
    mPower(H, eH, Q, &eQ, m, n);
    s=Q[(k-1)*m+k-1];
    for(i=1;i<=n;i++) {
        s = s*i/n;
        if (s<1e-140) {
            s*=1e140;
            eQ-=140;
        }
    }
    s*=pow(10., eQ);
    free(H);
    free(Q);
    return s;
}



/* Adapted from "Numerical Recipes in C" 3rd ed. pp. 334-335 & 736-738 */
double pks(double z) {
    if (z < 0.0)
        mexErrMsgTxt("z value for KS dist. must be positive.");
    if (z == 0.0)
        return 0.0;
    double z2 = z*z;
    if (z < 1.18) {
        double y = exp(-1.23370055013616983/z2); /* exp(-pi^2/(8*z^2)) */
        return 2.506628274631 / z * (y + pow(y,9) + pow(y,25) + pow(y,49)); /* sqrt(2*pi) */
    } else {
        double x = exp(-2*z2);
        return 1.0 - 2.0*(x - pow(x,4) - pow(x,9));
    }
}

/* Q_ks = 1 - P_ks */
double qks(double z) {
    if (z < 0.0)
        mexErrMsgTxt("z value for KS dist. must be positive.");
    if (z == 0.0)
        return 1.0;
    double z2 = z*z;
    if (z < 1.18) {
        double y = exp(-1.23370055013616983/z2); /* exp(-pi^2/(8*z^2)) */
        return 1.0 - 2.506628274631 / z * (y + pow(y,9) + pow(y,25) + pow(y,49)); /* sqrt(2*pi) */
    } else {
        double x = exp(-2*z2);
        return 2.0*(x - pow(x,4) - pow(x,9));
    }
}


int compDouble(const void *x, const void *y) {

    double dx = *(double *)x;
    double dy = *(double *)y;

    if (dx < dy) {
        return -1;
    } else if (dx > dy) {
        return +1;
    }
    return 0;
}


/* adapted from Numerical Recipes in C */
double ksone(double *data, int N, double mu, double sigma) {
    
    double *sdata = (double*)malloc(sizeof(double)*N);
    memcpy(sdata, data, N*sizeof(double));
    
    qsort(sdata, N, sizeof(double), compDouble);
    
    double normCDF;
    int j;
    double D = 0.0;
    double dt, fn, fo = 0.0;
    
    for (j=0; j<N; ++j) {
        fn = (j+1.0)/N; /* CDF */
        normCDF = 0.5 * (1.0 + erf((sdata[j]-mu)/(sqrt(2.0)*sigma)));
        
        dt = fmax(fabs(fo - normCDF), fabs(fn - normCDF));
        if (dt > D)
            D = dt;
        fo = fn;
    }
    free(sdata);
    /* return 1.0 - K(N,D); this is a much more accurate solution, but also much slower */
    double en = sqrt((double)N);
    return qks((en + 0.12 + 0.11/en)*D);
}


unsigned char adtest(double *data, int N, int adCase, double mu, double sigma, double alpha) {
    
    if (adCase<1 || adCase>4) {
        mexErrMsgTxt("'adCase' must be an integer between 1 and 4.");
    }
    
    double *sdata = (double*)malloc(sizeof(double)*N);
    memcpy(sdata, data, N*sizeof(double));
    qsort(sdata, N, sizeof(double), compDouble);
   
    double A2 = 0.0;
    double r2 = sqrt(2.0);
   
    int i;
   
    // Normal CDF
    double *z = (double*)malloc(sizeof(double)*N);
    for (i=0;i<N;++i) {
        z[i] = 0.5 * (1 + erf((sdata[i]-mu)/(r2*sigma)));
    }
        
    // A-D test statistic
    for (i=0;i<N;++i) {
        A2 += (2.0*i+1.0)*(log(z[i]) + log(1.0-z[N-1-i]));
    }
    A2 = -N - A2/N;
    
    if (adCase==4) {
        A2 *= 1.0+4.0/N-25.0/(N*N);
    }
    
    free(z);
    free(sdata);
    
    // Look-up table for critical values
    static const double alphaVec[4] = {0.01, 0.025, 0.05,  0.1};
    static const double ctable[16] = {3.857, 3.070, 2.492, 1.933, // case 1
                                      1.573, 1.304, 1.105, 0.908, // case 2
                                      3.690, 2.904, 2.323, 1.760, // case 3
                                      1.092, 0.918, 0.787, 0.656};// case 4
    
    // get critical value from lookup table
    int ai = 5;
    for (i=0;i<4;++i) {
        if (alpha==alphaVec[i]) {
            ai = i;
        }
    }
    if (ai==5) {
        mexErrMsgTxt("Admissible 'alpha' values: 0.01, 0.025, 0.05, 0.1");
    }
    return A2>ctable[ai + 4*(adCase-1)];
}
