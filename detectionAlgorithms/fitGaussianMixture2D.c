/* [prmVect prmStd covarianceMatrix residuals Jacobian] = fitGaussian2Dmex(prmVect, initValues, mode);
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

typedef int(*pfunc_t)(double*, int, argStruct_t*);

typedef struct dataStruct {
    int nx, np, ng, nDF, step;
    double *pixels;
    double *buffer;
    double *gx, *gy;
    int *estIdx;
    unsigned char binMode[5];
    int *idx;
    int nValid; // number of non-NaN pixels
    double *x_init;
    double *prmVect;
    pfunc_t *dfunc;
    gsl_vector *residuals;
    gsl_matrix *J;
    double *Jbuffer;
} dataStruct_t;


static int df_dx(double *J, int i, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    J[i] = A/s2*xi*g;
    return 0;
}

static int df_dy(double *J, int i, argStruct_t *argStruct) {
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    J[i] = A/s2*yi*g;
    return 0;
}

static int df_dA(double *J, int i, argStruct_t *argStruct) {
    J[i] = argStruct->g;
    return 0;
}

static int df_ds(double *J, int i, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s3 = argStruct->sigma3;
    J[i] += (xi*xi + yi*yi)*A/s3*g;
    return 0;
}

static int df_dc(double *J, int i, argStruct_t *argStruct) {
    J[i] = 1;
    return 0;
}



int gaussian_f(const gsl_vector *x, void *params, gsl_vector *f) {
    
    dataStruct_t *data = (dataStruct_t *)params;
    int nx = data->nx;
    int b = nx/2, i, k;
    int np = data->np;
    
    double *buffer = data->buffer;
    double *gx = data->gx;
    double *gy = data->gy;
    
    // update prmVect with new estimates
    for (i=0; i<data->nDF; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp, yp, A, xi, yi;
    double sigma = data->prmVect[np-2];
    double d = 2.0*sigma*sigma;
    double c = data->prmVect[np-1];
    
    int g, idx;
    div_t divRes;
    
    for (i=0; i<data->nValid; ++i) {
        idx = data->idx[i];
        divRes = div(idx, nx);
        buffer[i] = c - data->pixels[idx];
    }
    
    for (g=0; g<3*data->ng; g+=3) {
        xp = data->prmVect[g];
        yp = data->prmVect[g+1];
        A  = data->prmVect[g+2];
        
        for (i=0; i<nx; ++i) {
            k = i-b;
            xi = k-xp;
            yi = k-yp;
            gx[i] = exp(-xi*xi/d);
            gy[i] = exp(-yi*yi/d);
        }
        
        for (i=0; i<data->nValid; ++i) {
            idx = data->idx[i];
            divRes = div(idx, nx);
            buffer[i] += A*gx[divRes.quot]*gy[divRes.rem];
        }
    }
    
    for (i=0; i<data->nValid; ++i) {
        gsl_vector_set(f, i, buffer[i]);
    }
    
    return GSL_SUCCESS;
}



int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J) {
    
    dataStruct_t *data = (dataStruct_t *)params;
    // initialize Jacobian
    int nJ = data->nDF * data->nValid;
    memset(data->Jbuffer, 2.54, nJ*sizeof(double));
    
    int nx = data->nx;
    int b = nx/2, i, k;
    
    double *gx = data->gx;
    double *gy = data->gy;
    
    // update prmVect with new estimates
    for (i=0; i<data->nDF; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp, yp, xi, yi;
    double sigma = data->prmVect[data->np-2];
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    
    
    int g, idx;
    div_t divRes;
    
    for (g=0; g<3*data->ng; g+=3) {
        xp = data->prmVect[g];
        yp = data->prmVect[g+1];
        argStruct.A = data->prmVect[g+2];
        
        for (i=0; i<nx; ++i) {
            k = i-b;
            xi = k-xp;
            yi = k-yp;
            gx[i] = exp(-xi*xi/d);
            gy[i] = exp(-yi*yi/d);
        }
        
        for (i=0; i<data->nValid; ++i) {
            idx = data->idx[i];
            divRes = div(idx, nx);
            argStruct.xi = divRes.quot-b - xp;
            argStruct.yi = divRes.rem-b - yp;
            argStruct.g = gx[divRes.quot]*gy[divRes.rem];
            
            for (k=0; k<data->step; ++k) {
                data->dfunc[g+k](data->Jbuffer, i+(g+k)*data->nValid, &argStruct);
            }
            for (k=data->ng*data->step; k<data->nDF; ++k) {
                data->dfunc[k](data->Jbuffer, i+k*data->nValid, &argStruct);
            }
        }
    }
    
    // copy Jacobian
    for (i=0; i<nJ; ++i) {
        divRes = div(i, data->nValid);
        gsl_matrix_set(J, divRes.rem, divRes.quot, data->Jbuffer[i]);
    }
    return GSL_SUCCESS;
}



int gaussian_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
    
    dataStruct_t *data = (dataStruct_t *)params;
    // initialize Jacobian
    int nJ = data->nDF * data->nValid;
    memset(data->Jbuffer, 2.54, nJ*sizeof(double));
    
    int nx = data->nx;
    int b = nx/2, i, k;
    
    double *buffer = data->buffer;
    double *gx = data->gx;
    double *gy = data->gy;
    
    // update prmVect with new estimates
    for (i=0; i<data->nDF; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp, yp, xi, yi, A;
    double sigma = data->prmVect[data->np-2];
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    double c = data->prmVect[data->np-1];
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    
    int g, idx;
    div_t divRes;
    
    for (i=0; i<data->nValid; ++i) {
        idx = data->idx[i];
        divRes = div(idx, nx);
        buffer[i] = c - data->pixels[idx];
    }
    
    for (g=0; g<3*data->ng; g+=3) {
        xp = data->prmVect[g];
        yp = data->prmVect[g+1];
        A  = data->prmVect[g+2];
        argStruct.A = A;
        
        for (i=0; i<nx; ++i) {
            k = i-b;
            xi = k-xp;
            yi = k-yp;
            gx[i] = exp(-xi*xi/d);
            gy[i] = exp(-yi*yi/d);
        }
        
        for (i=0; i<data->nValid; ++i) {
            idx = data->idx[i];
            divRes = div(idx, nx);
            
            buffer[i] += A*gx[divRes.quot]*gy[divRes.rem];
            
            argStruct.xi = divRes.quot-b - xp;
            argStruct.yi = divRes.rem-b - yp;
            argStruct.g = gx[divRes.quot]*gy[divRes.rem];
            
            for (k=0; k<data->step; ++k) {
                data->dfunc[g+k](data->Jbuffer, i+(g+k)*data->nValid, &argStruct);
            }
            for (k=data->ng*data->step; k<data->nDF; ++k) {
                data->dfunc[k](data->Jbuffer, i+k*data->nValid, &argStruct);
            }
        }
    }
    
    for (i=0; i<data->nValid; ++i) {
        gsl_vector_set(f, i, buffer[i]);
    }
    
    // copy Jacobian
    for (i=0; i<nJ; ++i) {
        divRes = div(i, data->nValid);
        gsl_matrix_set(J, divRes.rem, divRes.quot, data->Jbuffer[i]);
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
    
    gsl_rng_env_setup();
    
    type = gsl_rng_default;
    
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
    for (i=0; i<data->nDF; ++i) {
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



int main(void) {
    
    
    
    int nx = 15;
    int nx2 = nx*nx;
    double* px;
    px = (double*)malloc(sizeof(double)*nx2);
    
    
    
    int i;
    // fill with noise
    for (i=0; i<nx2; ++i) {
        px[i] = rand();
    }
    
    int np = 5;
    
    dataStruct_t data;
    data.nx = nx;
    data.np = 5;
    data.pixels = px;
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*nx);
    
    data.estIdx = (int*)malloc(sizeof(int)*np);
    //memcpy(data.prmVect, mxGetPr(prhs[1]), 5*sizeof(double));
    data.dfunc = (pfunc_t*) malloc(sizeof(pfunc_t) * np);
    
    // read mask/pixels
    data.nValid = nx2;
    data.idx = (int*)malloc(sizeof(int)*data.nValid);
    int k = 0;
    for (i=0; i<nx2; ++i) {
        if (!mxIsNaN(data.pixels[i])) {
            data.idx[k++] = i;
        }
    }
    
    data.prmVect[0] = 0;
    data.prmVect[1] = 0;
    data.prmVect[2] = 5;
    data.prmVect[3] = 1;
    data.prmVect[4] = 0;
    
    np = 0;
    data.estIdx[np] = 0; data.dfunc[np++] = df_dx;
    data.estIdx[np] = 1; data.dfunc[np++] = df_dy;
    data.estIdx[np] = 2; data.dfunc[np++] = df_dA;
    data.estIdx[np] = 3; data.dfunc[np++] = df_ds;
    data.estIdx[np] = 4; data.dfunc[np++] = df_dc;
    
    data.x_init = (double*)malloc(sizeof(double)*np);
    for (i=0; i<np; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }
    
    MLalgo(&data);
    
    gsl_matrix_free(data.J);
    gsl_vector_free(data.residuals);
    
    free(data.x_init);
    free(data.idx);
    free(data.dfunc);
    free(data.estIdx);
    free(data.gy);
    free(data.gx);
    
    free(px);
    return 0;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* inputs:
     * image array
     * prmVect
     * mode
     */
    
    
    // check inputs
    if (nrhs < 3) mexErrMsgTxt("Inputs should be: data, prmVect, mode.");
    if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Data input must be double array.");
    int np = mxGetNumberOfElements(prhs[1]);
    div_t divRes = div(np-2, 3);
    if (divRes.rem != 0) mexErrMsgTxt("Invalid parameter vector length.");
    
    if (!mxIsDouble(prhs[1])) mexErrMsgTxt("Parameter vector must be double array.");
    if (!mxIsChar(prhs[2])) mexErrMsgTxt("Mode needs to be a string.");
    
    
    
    size_t nx = mxGetN(prhs[0]);
    int nx2 = nx*nx;
    int ng = divRes.quot;
    int i;
    
    // read mode input
    int nm = mxGetNumberOfElements(prhs[2]);
    char *mode;
    mode = (char*)malloc(sizeof(char)*(nm+1));
    mxGetString(prhs[2], mode, nm+1);
    for (i=0; i<strlen(mode); ++i) {
        mode[i] = tolower(mode[i]);
    }
    
    dataStruct_t data;
    
    // detect parameters to optimize
    char *refMode = "xyasc";
    for (i=0; i<5; ++i) {
        data.binMode[i] = 1;
        (strchr(mode, refMode[i])!=NULL) ? (data.binMode[i]=1) : (data.binMode[i]=0);
    }
    int step = data.binMode[0] + data.binMode[1] + data.binMode[2];
    int nModes = step + data.binMode[3] + data.binMode[4];
    int nDF = nModes + step*(ng-1);
    data.nDF = nDF;
    data.step = step;
    
    
    // allocate
    data.nx = nx;
    data.np = np;
    data.ng = ng;
    data.pixels = mxGetPr(prhs[0]);
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*nx);
    data.estIdx = (int*)malloc(sizeof(int)*nDF);
    data.prmVect = mxGetPr(prhs[1]);
    data.dfunc = (pfunc_t*)malloc(sizeof(pfunc_t)*nDF);
    
    
    // read mask/pixels
    data.nValid = nx2;
    for (i=0; i<nx2; ++i) {
        data.nValid -= (int)mxIsNaN(data.pixels[i]);
    }
    data.buffer = (double*)malloc(sizeof(double)*data.nValid);
    data.idx = (int*)malloc(sizeof(int)*data.nValid);
    int k = 0;
    for (i=0; i<nx2; ++i) {
        if (!mxIsNaN(data.pixels[i])) {
            data.idx[k++] = i;
        }
    }
    
    np = 0;
    
    for (i=0; i<ng; ++i) {
        if (strchr(mode, 'x')!=NULL) {
            data.estIdx[np] = 3*i;
            data.dfunc[np++] = df_dx;
        }
        if (strchr(mode, 'y')!=NULL) {
            data.estIdx[np] = 1+3*i;
            data.dfunc[np++] = df_dy;
        }
        if (strchr(mode, 'a')!=NULL) {
            data.estIdx[np] = 2+3*i;
            data.dfunc[np++] = df_dA;
        }
    }
    int pos = step*ng;
    if (strchr(mode, 's')!=NULL) {data.estIdx[pos++] = 3*ng; data.dfunc[np++] = df_ds;}
    if (strchr(mode, 'c')!=NULL) {data.estIdx[pos] = 3*ng+1; data.dfunc[np++] = df_dc;}
    
    np = data.np;
    
    data.x_init = (double*)malloc(sizeof(double)*nDF);
    for (i=0; i<nDF; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }
    data.Jbuffer = (double*)malloc(sizeof(double)*nDF*data.nValid);
    
    MLalgo(&data);
    
    // parameters
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(1, np, mxREAL);
        memcpy(mxGetPr(plhs[0]), data.prmVect, np*sizeof(double));
    }
    
    // standard dev. of parameters & covariance matrix
    if (nlhs > 1) {
        
        gsl_matrix *covar = gsl_matrix_alloc(np, np);
        gsl_multifit_covar(data.J, 0.0, covar);
        double sigma_e = 0.0, e;
        for (i=0; i<data.nValid; ++i) {
            e = data.residuals->data[data.idx[i]];
            sigma_e += e*e;
        }
        sigma_e /= data.nValid - data.np - 1;
        plhs[1] = mxCreateDoubleMatrix(1, data.np, mxREAL);
        double *prmStd = mxGetPr(plhs[1]);
        for (i=0; i<data.np; ++i) {
            prmStd[i] = sqrt(sigma_e*gsl_matrix_get(covar, i, i));
        }
        if (nlhs > 2) {
            plhs[2] = mxCreateDoubleMatrix(np, np, mxREAL);
            // cov. matrix is symmetric, no need to transpose
            memcpy(mxGetPr(plhs[2]), covar->data, np*np*sizeof(double));
        }
        free(covar);
        
    }
    
    // residuals
    if (nlhs > 3) {
        plhs[3] = mxCreateDoubleMatrix(nx, nx, mxREAL);
        memcpy(mxGetPr(plhs[3]), data.residuals->data, nx2*sizeof(double));
    }
    
    // Jacobian
    if (nlhs > 4) {
        // convert row-major double* data.J->data to column-major double*
        plhs[4] = mxCreateDoubleMatrix(nx2, np, mxREAL);
        double *J = mxGetPr(plhs[4]);
        
        int p;
        for (p=0; p<np; ++p) {
            for (i=0; i<nx2; ++i) {
                J[i+p*nx2] = (data.J)->data[p + i*np];
            }
        }
    }
    
    gsl_matrix_free(data.J);
    gsl_vector_free(data.residuals);
    free(data.Jbuffer);
    free(data.x_init);
    free(data.idx);
    free(data.buffer);
    free(data.dfunc);
    free(data.estIdx);
    free(data.gy);
    free(data.gx);
    free(mode);
}