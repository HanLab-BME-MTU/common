/* [prmVect covarianceMatrix residuals Jacobian] = fitGaussian2Dmex(prmVect, initValues, mode);
 *
 * (c) Francois Aguet & Sylvain Berlemont, 2011 (last modified Jan 22, 2011)
 *
 * Compile with: mex -I/usr/local/include -lgsl -lgslcblas fitGaussian2Dmex.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h> // tolower()

#include <regex.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "mex.h"
#include "matrix.h"

typedef struct argStruct {
    double xi, yi, A, g, sigma2, sigma3;
} argStruct_t;

typedef int(*pfunc_t)(gsl_matrix*, int, int, argStruct_t*);

typedef struct dataStruct {
    int nx, np;
    double *pixels;
    double *gx, *gy;
    int *estIdx;
    double *x_init;
    double prmVect[5];
    pfunc_t *dfunc;
    gsl_vector *residuals;
    gsl_matrix *J;
} dataStruct_t;



static int df_dx(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    gsl_matrix_set(J, i, k, A/s2*xi*g);
    return 0;
}

static int df_dy(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    gsl_matrix_set(J, i, k, A/s2*yi*g);
    return 0;
}

static int df_dA(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    gsl_matrix_set(J, i, k, argStruct->g);
    return 0;
}

static int df_ds(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s3 = argStruct->sigma3;
    gsl_matrix_set(J, i, k, (xi*xi + yi*yi)*A/s3*g);
    return 0;
}

static int df_dc(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    gsl_matrix_set(J, i, k, 1);
    return 0;
}



int gaussian_f(const gsl_vector *x, void *params, gsl_vector *f) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int b = nx/2, i, k;
    
    double *pixels = dataStruct->pixels;
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    
    // update prmVect with new estimates
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double A = dataStruct->prmVect[2];
    double sigma = dataStruct->prmVect[3];
    double c = dataStruct->prmVect[4];

    double xi, yi;
    double d = 2.0*sigma*sigma;
    for (i=0; i<nx; ++i) {
        k = i-b;
        xi = k-xp;
        yi = k-yp;
        gx[i] = exp(-xi*xi/d);
        gy[i] = exp(-yi*yi/d);
    }
    
    int nPixels = nx*nx;
    div_t divRes;
    for (i=0; i<nPixels; ++i) {
        divRes = div(i, nx);
        gsl_vector_set(f, i, A*gx[divRes.quot]*gy[divRes.rem]+c - pixels[i]);
    }
    return GSL_SUCCESS;
}



int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int b = nx/2, i, k;
    
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    
    // update prmVect with new estimates
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double A = dataStruct->prmVect[2];
    double sigma = dataStruct->prmVect[3];
        
    double xi, yi;
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    argStruct.A = A;
    
    for (i=0; i<nx; ++i) {
        k = i-b;
        xi = k-xp;
        yi = k-yp;
        gx[i] = exp(-xi*xi/d);
        gy[i] = exp(-yi*yi/d);
    }
    
    int nPixels = nx*nx;
    div_t divRes;
    for (i=0; i<nPixels; ++i) {
        divRes = div(i, nx);
        argStruct.xi = divRes.quot-b - xp;
        argStruct.yi = divRes.rem-b - yp;
        argStruct.g = gx[divRes.quot]*gy[divRes.rem];
        
        for (k=0; k<dataStruct->np; ++k)
            dataStruct->dfunc[k](J, i, k, &argStruct);
    }
    return GSL_SUCCESS;
}



int gaussian_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {

    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int b = nx/2, i, k;
    
    double *pixels = dataStruct->pixels;
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    
    // update prmVect with new estimates
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double A = dataStruct->prmVect[2];
    double sigma = dataStruct->prmVect[3];
    double c = dataStruct->prmVect[4];
        
    double xi, yi;
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    argStruct.A = A;
    
    for (i=0; i<nx; ++i) {
        k = i-b;
        xi = k-xp;
        yi = k-yp;
        gx[i] = exp(-xi*xi/d);
        gy[i] = exp(-yi*yi/d);
    }
    
    int nPixels = nx*nx;
    div_t divRes;
    for (i=0; i<nPixels; ++i) {
        divRes = div(i, nx);
        gsl_vector_set(f, i, A*gx[divRes.quot]*gy[divRes.rem]+c - pixels[i]);
        argStruct.xi = divRes.quot-b - xp;
        argStruct.yi = divRes.rem-b - yp;
        argStruct.g = gx[divRes.quot]*gy[divRes.rem];
        
        for (k=0; k<dataStruct->np; ++k)
            dataStruct->dfunc[k](J, i, k, &argStruct);
    }
    return GSL_SUCCESS;
}



int MLalgo(struct dataStruct *data) {
    
    // declare solvers
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    
    // number of parameters to optimize
    const size_t p = data->np;
    
    gsl_vector_view x = gsl_vector_view_array(data->x_init, p);
    
    const gsl_rng_type *type;
    gsl_rng *r;
    
    gsl_rng_env_setup();
    
    type = gsl_rng_default;
    r = gsl_rng_alloc(type);
    
    gsl_multifit_function_fdf f;
    f.f = &gaussian_f;
    f.df = &gaussian_df;
    f.fdf = &gaussian_fdf;
    size_t n = data->nx;
    n *= n;
    f.n = n;
    f.p = p;
    f.params = data;
    
    
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, n, p);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);
    
    int status;
    int iter = 0;
    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);
        if (status)
            break;
        
        status = gsl_multifit_test_delta(s->dx, s->x, 1e-8, 1e-8);
    }
    while (status == GSL_CONTINUE && iter < 500);
    
    int i;
    for (i=0; i<data->np; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(s->x, i);
    }
    
    // copy model
    data->residuals = gsl_vector_alloc(n);
    gsl_vector_memcpy(data->residuals, s->f);
    
    // copy Jacobian
    data->J = gsl_matrix_alloc(n, data->np);
    gsl_matrix_memcpy(data->J, s->J);
    
    gsl_multifit_fdfsolver_free(s);
    return 0;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
       
    /* inputs:
     * image array
     * prmVect
     * mode
     *
     */
    
    
    // check inputs
    if (nrhs < 3) mexErrMsgTxt("Inputs should be: data, prmVect, mode.");
    if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Data input must be double array.");
    if (mxGetNumberOfElements(prhs[1])!=5 || !mxIsDouble(prhs[1])) mexErrMsgTxt("Incorrect parameter vector format.");
    if (!mxIsChar(prhs[2])) mexErrMsgTxt("Mode needs to be a string.");
    
    size_t nx = mxGetN(prhs[0]);
    int nx2 = nx*nx;
    
    // read mode input
    int np = mxGetNumberOfElements(prhs[2]);
    char *mode;
    mode = (char*)malloc(sizeof(char)*(np+1));
    mxGetString(prhs[2], mode, np+1);
    
    int i;
    for (i=0; i<strlen(mode); ++i) {
        mode[i] = tolower(mode[i]);
    }
    
    char *refMode = "xyasc";
    np = 0; // number of parameters to fit
    for (i=0; i<5; ++i) {
        if (strchr(mode, refMode[i])!=NULL) { np++; }
    }
    
    // allocate
    dataStruct_t data;
    data.nx = nx;
    data.np = np;
    data.pixels = mxGetPr(prhs[0]);
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*nx);
    //bzero(data.estIdx, 5);
    data.estIdx = (int*)malloc(sizeof(int)*np);
    memcpy(data.prmVect, mxGetPr(prhs[1]), 5*sizeof(double));    
    data.dfunc = (pfunc_t*) malloc(sizeof(pfunc_t) * np);
    
    np = 0;
    if (strchr(mode, 'x')!=NULL) {data.estIdx[np] = 0; data.dfunc[np++] = df_dx;}
    if (strchr(mode, 'y')!=NULL) {data.estIdx[np] = 1; data.dfunc[np++] = df_dy;}
    if (strchr(mode, 'a')!=NULL) {data.estIdx[np] = 2; data.dfunc[np++] = df_dA;}
    if (strchr(mode, 's')!=NULL) {data.estIdx[np] = 3; data.dfunc[np++] = df_ds;}
    if (strchr(mode, 'c')!=NULL) {data.estIdx[np] = 4; data.dfunc[np++] = df_dc;}
    
    data.x_init = (double*)malloc(sizeof(double)*np);
    for (i=0; i<np; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }   
    
    MLalgo(&data);
    np = data.np;
    
    // parameters
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(1, 5, mxREAL);
        memcpy(mxGetPr(plhs[0]), data.prmVect, 5*sizeof(double));
    }
    
    // covariance matrix
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(np, np, mxREAL);
        gsl_matrix *covar = gsl_matrix_alloc(np, np);
        gsl_multifit_covar(data.J, 0.0, covar);
        // cov. matrix is symmetric, no need to transpose
        memcpy(mxGetPr(plhs[1]), covar->data, np*np*sizeof(double));
        free(covar);
    }
    
    // residuals
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(nx, nx, mxREAL);
        memcpy(mxGetPr(plhs[2]), data.residuals->data, nx2*sizeof(double));
    }
    
    // Jacobian
    if (nlhs > 3) {
        // convert row-major double* data.J->data to column-major double*
        plhs[3] = mxCreateDoubleMatrix(nx2, np, mxREAL);
        double *J = mxGetPr(plhs[3]);
        
        int p;
        for (p=0; p<np; ++p) {
            for (i=0; i<nx2; ++i) {
                J[i+p*nx2] = (data.J)->data[p + i*np];
            }
        }
    }
    
    free(data.gx);
    free(data.gy);
    free(data.dfunc);
    free(data.x_init);
    free(data.residuals);
    free(data.J);
    free(data.estIdx);
    free(mode);
}