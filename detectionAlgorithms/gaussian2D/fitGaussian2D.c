/*
 * fitGaussian2D.c  —  MEX interface using GSL LM solver
 *
 * MATLAB signature:
 *   [prmVect, prmStd, C, res, J] = fitGaussian2D(data, prmVect, mode, options)
 *
 * Inputs
 *   data    : 2D image (NaNs are masked)
 *   prmVect : [xp yp A s c] initial guess (xp,yp origin = image center)
 *   mode    : string subset of 'xyasc' selecting which params to fit
 *   options : [maxIter eAbs eRel] (optional; defaults: 200, 1e-6, 1e-6)
 *
 * Outputs
 *   prmVect : 1×5 vector with fitted values in selected positions
 *   prmStd  : 1×np std devs for fitted parameters only (np = length(mode))
 *   C       : np×np covariance for fitted parameters only
 *   res     : struct with fields .data (m×1), .pval (NaN), .mean, .std, .RSS
 *   J       : m×np Jacobian at solution (m = #valid pixels (non-NaN))
 *
 * Compile (example macOS + Homebrew GSL):
clear mex
mex('-v','-g', ...
    'CFLAGS=$CFLAGS -O0 -fno-omit-frame-pointer -std=c11', ...
    'fitGaussian2D.c', ...
    '-I/opt/homebrew/opt/gsl/include', ...
    '-L/opt/homebrew/opt/gsl/lib', ...
    '-lgsl','-lgslcblas','-lm');
 */

#include "mex.h"
#include "matrix.h"

#include <math.h>
#include <string.h>
#include <ctype.h>  /* tolower */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
/* ADD: the following two headers provide gsl_vector_size, gsl_matrix_size1/2, etc. */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* ------------------------ data container ------------------------ */

typedef struct {
    /* image geometry */
    size_t w, h;         /* width (cols), height (rows) */
    size_t m;            /* # valid pixels */

    /* coordinate arrays (relative to center) and observations */
    double *X;           /* length m: x (meshgrid convention) */
    double *Y;           /* length m: y */
    double *I;           /* length m: observed intensity */

    /* selection of parameters (subset of [xp yp A s c] = [0..4]) */
    size_t np;           /* # fitted params = length(estIdx) */
    int    estIdx[5];    /* indices into prmVect, e.g., {0,1,2} */

    /* solver controls */
    size_t maxIter;
    double eAbs;
    double eRel;

    /* full parameter vector: [xp, yp, A, s, c] */
    double prmVect[5];

} FitData;

/* --------------------- model & derivatives ---------------------- */
/* Gaussian model: g = c + A * exp( -((x-xp)^2 + (y-yp)^2) / (2*s^2) ) */

static inline double gauss_eval(double xp, double yp, double A, double s, double c,
                                double x, double y)
{
    const double dx = x - xp;
    const double dy = y - yp;
    const double s2 = s * s;
    const double e  = exp(-(dx*dx + dy*dy) / (2.0 * s2));
    return c + A * e;
}

/* f: residuals = observed - model */
static int gaussian_f(const gsl_vector *x, void *params, gsl_vector *fvec)
{
    FitData *D = (FitData*)params;

    /* expand x (fitted subset) back into full prmVect */
    double pv[5];
    for (int i = 0; i < 5; ++i) pv[i] = D->prmVect[i];
    for (size_t k = 0; k < D->np; ++k) {
        pv[D->estIdx[k]] = gsl_vector_get(x, k);
    }
    const double xp = pv[0], yp = pv[1], A = pv[2], s = pv[3], c = pv[4];

    const size_t m = D->m;
    for (size_t i = 0; i < m; ++i) {
        const double xpi = D->X[i];
        const double ypi = D->Y[i];
        const double yobs = D->I[i];
        const double g = gauss_eval(xp, yp, A, s, c, xpi, ypi);
        gsl_vector_set(fvec, i, yobs - g);
    }
    return GSL_SUCCESS;
}

/* Model: f(x,y) = c + A * exp( - ((X-x0)^2 + (Y-y0)^2) / (2*s^2) )
   Params (full order): [x0, y0, A, s, c]
   Jacobian columns (full):
     ∂f/∂x0 =  A * E * ( (X-x0) / s^2 )
     ∂f/∂y0 =  A * E * ( (Y-y0) / s^2 )
     ∂f/∂A  =  E
     ∂f/∂s  =  A * E * ( r2 / s^3 )
     ∂f/∂c  =  1
*/
/* FULL param order we target everywhere: [x0, y0, A, s, c]
   Model: f(x,y) = c + A * exp(-((X-x0)^2 + (Y-y0)^2) / (2*s^2))
   Residual convention used elsewhere: r = yobs - f(x,y)  (so J = ∂r/∂θ = -∂f/∂θ)
*/
static int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
    FitData *D = (FitData*)params;

    /* --- 1) Rebuild full parameter vector from reduced x + fixed values --- */
    double pfull[5]; /* [x0, y0, A, s, c] */
    /* start from current full vector (contains fixed values) */
    pfull[0] = D->prmVect[0];
    pfull[1] = D->prmVect[1];
    pfull[2] = D->prmVect[2];
    pfull[3] = D->prmVect[3];
    pfull[4] = D->prmVect[4];

    /* overwrite estimated ones from reduced vector */
    for (size_t j = 0; j < D->np; ++j) {
        int full = D->estIdx[j]; /* 0..4 */
        pfull[full] = gsl_vector_get(x, j);
    }

    double x0 = pfull[0];
    double y0 = pfull[1];
    double A  = pfull[2];
    double s  = pfull[3];
    /* c = pfull[4];  // not needed for derivatives below */

    /* guard s */
    const double s_eps = 1e-9;
    if (!(s > s_eps)) s = s_eps;
    const double s2 = s*s;
    const double inv_s2 = 1.0 / s2;
    const double inv_s3 = 1.0 / (s2*s);

    /* --- 2) Geometry: compute (X,Y) coordinates from linear index k --- */
    const size_t W = D->w;  /* number of columns */
    const size_t H = D->h;  /* number of rows    */
    const double cx = 0.5 * (double)(W - 1); /* center column index */
    const double cy = 0.5 * (double)(H - 1); /* center row index    */

    /* J must be (m x np), where m = number of residuals (usually W*H or #valid) */
    /* If you pack only valid pixels elsewhere, loop over those indices here.
       Otherwise, fill all pixels. This version fills ALL W*H pixels. */

    size_t row = 0;
    for (size_t r = 0; r < H; ++r) {
        const double Y = (double)r - cy;
        for (size_t c = 0; c < W; ++c, ++row) {
            const double X = (double)c - cx;

            /* exp term */
            const double dx = X - x0;
            const double dy = Y - y0;
            const double r2 = dx*dx + dy*dy;
            const double E  = exp(-0.5 * r2 * inv_s2);

            /* fill each estimated column */
            for (size_t j = 0; j < D->np; ++j) {
                double dfdtheta = 0.0;
                switch (D->estIdx[j]) {
                    case 0: /* ∂f/∂x0 */ dfdtheta =  A * E * ( dx * inv_s2 );        break;
                    case 1: /* ∂f/∂y0 */ dfdtheta =  A * E * ( dy * inv_s2 );        break;
                    case 2: /* ∂f/∂A  */ dfdtheta =  E;                               break;
                    case 3: /* ∂f/∂s  */ dfdtheta =  A * E * ( r2 * inv_s3 );         break;
                    case 4: /* ∂f/∂c  */ dfdtheta =  1.0;                             break;
                    default: dfdtheta = 0.0;                                          break;
                }
                /* residual r = yobs - f  ->  J = ∂r/∂θ = -∂f/∂θ */
                gsl_matrix_set(J, row, j, -dfdtheta);
            }
        }
    }

    /* Optional sanity check: warn if any Jacobian column is exactly zero */
    /*
    for (size_t j = 0; j < D->np; ++j) {
        double norm2 = 0.0;
        for (size_t i = 0; i < row; ++i) {
            double v = gsl_matrix_get(J,i,j);
            norm2 += v*v;
        }
        if (norm2 == 0.0) {
            mexWarnMsgIdAndTxt("fitGaussian2D:ZeroJacCol",
                "Jacobian column %zu (full param %d) is zero.", j, (int)D->estIdx[j]);
        }
    }
    */
    return GSL_SUCCESS;
}
/* Model: f(x,y) = c + A * exp( - ((X-x0)^2 + (Y-y0)^2) / (2*s^2) )
   Params (full order): [x0, y0, A, s, c]
   Jacobian columns (full):
     ∂f/∂x0 =  A * E * ( (X-x0) / s^2 )
     ∂f/∂y0 =  A * E * ( (Y-y0) / s^2 )
     ∂f/∂A  =  E
     ∂f/∂s  =  A * E * ( r2 / s^3 )
     ∂f/∂c  =  1
*/
/* J: Jacobian of residuals wrt fitted params (m×np) */
/*static int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
    FitData *D = (FitData*)params;

    double pv[5];
    for (int i = 0; i < 5; ++i) pv[i] = D->prmVect[i];
    for (size_t k = 0; k < D->np; ++k) {
        pv[D->estIdx[k]] = gsl_vector_get(x, k);
    }
    const double xp = pv[0], yp = pv[1], A = pv[2], s = pv[3], c = pv[4];
    (void)c; /* c used in eval, but derivs below don't need it */

/*    const size_t m = D->m;
    const double s2 = s*s;
    const double s3 = s2*s;

    for (size_t i = 0; i < m; ++i) {
        const double xpi = D->X[i];
        const double ypi = D->Y[i];
        const double dx  = xpi - xp;
        const double dy  = ypi - yp;
        const double e   = exp(-(dx*dx + dy*dy) / (2.0 * s2));

        /* residual r = y - g; J = dr/dtheta = -dg/dtheta */
        /*for (size_t j = 0; j < D->np; ++j) {
            double d = 0.0;
            switch (D->estIdx[j]) {
                case 0: /* xp */
          /*          d = - (A * e) * (dx / s2);
                    break;
                case 1: /* yp */
            /*        d = - (A * e) * (dy / s2);
                    break;
                case 2: /* A  */
            /*        d = - e;
                    break;
                case 3: /* s  */
            /*        d = - (A * e) * ( (dx*dx + dy*dy) / s3 );
                    break;
                case 4: /* c  */
           /*         d = - 1.0;
                    break;
                default:
                    d = 0.0;
                    break;
            }
            gsl_matrix_set(J, i, j, d);
        }
    }
    return GSL_SUCCESS;
}
*/

static int gaussian_fdf(const gsl_vector *x, void *params, gsl_vector *fvec, gsl_matrix *J)
{
    gaussian_f(x, params, fvec);
    gaussian_df(x, params, J);
    return GSL_SUCCESS;
}

/* ------------------------ utilities ----------------------------- */

static void parse_mode(const mxArray *mxMode, FitData *D)
{
    if (!mxIsChar(mxMode))
        mexErrMsgIdAndTxt("fitGaussian2D:modeType", "mode must be a char array.");

    char buf[16];
    mxGetString(mxMode, buf, sizeof(buf));
    /* normalize to lowercase, keep only xyasc in order encountered */
    int used[5] = {0,0,0,0,0};
    D->np = 0;
    for (size_t i = 0; i < strlen(buf); ++i) {
        char c = (char)tolower((int)buf[i]);
        int idx = -1;
        if (c=='x') idx = 0;
        else if (c=='y') idx = 1;
        else if (c=='a') idx = 2;
        else if (c=='s') idx = 3;
        else if (c=='c') idx = 4;
        if (idx>=0 && !used[idx]) {
            D->estIdx[D->np++] = idx;
            used[idx] = 1;
        }
    }
    if (D->np == 0)
        mexErrMsgIdAndTxt("fitGaussian2D:modeEmpty", "mode must contain at least one of 'x','y','a','s','c'.");
}

static void parse_options(const mxArray *mxOpt, FitData *D)
{
    /* defaults */
    D->maxIter = 200;
    D->eAbs    = 1e-6;
    D->eRel    = 1e-6;

    if (mxOpt == NULL) return;
    if (!mxIsDouble(mxOpt) || mxIsComplex(mxOpt) || mxGetNumberOfElements(mxOpt) < 1)
        mexErrMsgIdAndTxt("fitGaussian2D:optionsType", "options must be a real double vector [maxIter eAbs eRel].");

    double *v = mxGetPr(mxOpt);
    size_t n = mxGetNumberOfElements(mxOpt);
    if (n >= 1) D->maxIter = (size_t)( (v[0] > 0) ? v[0] : 200 );
    if (n >= 2) D->eAbs    = (v[1] > 0) ? v[1] : 1e-6;
    if (n >= 3) D->eRel    = (v[2] > 0) ? v[2] : 1e-6;
}

static void extract_valid_pixels(const mxArray *mxImg, FitData *D)
{
    if (!mxIsDouble(mxImg) || mxIsComplex(mxImg) || mxGetNumberOfDimensions(mxImg) != 2)
        mexErrMsgIdAndTxt("fitGaussian2D:dataType", "data must be a real 2D double matrix.");

    const mwSize *dims = mxGetDimensions(mxImg);
    D->h = (size_t)dims[0];
    D->w = (size_t)dims[1];
    const double *Iall = mxGetPr(mxImg);

    /* count valid (non-NaN) */
    size_t m = 0;
    for (size_t j = 0; j < D->w; ++j) {
        for (size_t i = 0; i < D->h; ++i) {
            double val = Iall[i + j*(size_t)D->h];
            if (!mxIsNaN(val)) ++m;
        }
    }
    D->m = m;

    D->X = (double*)mxCalloc(m, sizeof(double));
    D->Y = (double*)mxCalloc(m, sizeof(double));
    D->I = (double*)mxCalloc(m, sizeof(double));

    const double cx = (double)(D->w - 1) / 2.0;
    const double cy = (double)(D->h - 1) / 2.0;

    size_t k = 0;
    for (size_t j = 0; j < D->w; ++j) {
        for (size_t i = 0; i < D->h; ++i) {
            double val = Iall[i + j*(size_t)D->h];
            if (!mxIsNaN(val)) {
                /* meshgrid convention: x varies with columns j, y with rows i */
                D->X[k] = (double)j - cx;
                D->Y[k] = (double)i - cy;
                D->I[k] = val;
                ++k;
            }
        }
    }
}

static void copy_outputs_and_cleanup(
    const FitData *D,
    gsl_multifit_fdfsolver *s,
    int status,
    mxArray *plhs[])
{
    /* ---- 0) prmVect (1x5) full vector with fitted subset inserted ---- */
    plhs[0] = mxCreateDoubleMatrix(1, 5, mxREAL);
    double *out0 = mxGetPr(plhs[0]);

    /* current solution vector x (length np) */
    const gsl_vector *xvec = gsl_multifit_fdfsolver_position(s);

    /* start with prior pv and insert fitted subset */
    double pv[5];
    for (int i = 0; i < 5; ++i) pv[i] = D->prmVect[i];
    for (size_t k = 0; k < D->np; ++k) {
        pv[D->estIdx[k]] = gsl_vector_get(xvec, k);
    }
    /* enforce s >= 0 */
    if (pv[3] < 0.0) pv[3] = fabs(pv[3]);

    for (int i = 0; i < 5; ++i) out0[i] = pv[i];

    /* ---------------- sizes and covariance for fitted params ---------------- */

    const size_t np = xvec->size;

    const gsl_vector *F = gsl_multifit_fdfsolver_residual(s);
    const size_t m = F->size;

    gsl_matrix *Jtmp = gsl_matrix_alloc(m, np);
    gsl_multifit_fdfsolver_jac((gsl_multifit_fdfsolver *)s, Jtmp);

    double nrm = gsl_blas_dnrm2(F); /* ||F|| */
    double chisq = nrm * nrm;
    double dof = (m > np) ? (double)(m - np) : 1.0;
    double s2  = chisq / dof;

    gsl_matrix *cov_np = gsl_matrix_alloc(np, np);
    gsl_multifit_covar(Jtmp, s2, cov_np);
    gsl_matrix_free(Jtmp);

    /* ---- 1) prmStd (1xnp) ---- */
    plhs[1] = mxCreateDoubleMatrix(1, np, mxREAL);
    double *out1 = mxGetPr(plhs[1]);
    for (size_t j = 0; j < np; ++j) {
        double v = gsl_matrix_get(cov_np, j, j);
        out1[j] = (v > 0.0) ? sqrt(v) : mxGetNaN();
    }

    /* ---- 2) C (np x np) ---- */
    plhs[2] = mxCreateDoubleMatrix((mwSize)np, (mwSize)np, mxREAL);
    double *out2 = mxGetPr(plhs[2]);
    for (size_t r = 0; r < np; ++r) {
        for (size_t c = 0; c < np; ++c) {
            out2[r + c*np] = gsl_matrix_get(cov_np, r, c);
        }
    }

    /* ---- 3) res struct ---- */
    const char *fnames[] = {"data","pval","mean","std","RSS"};
    plhs[3] = mxCreateStructMatrix(1, 1, 5, fnames);

    /* res.data (m x 1) = residuals (observed - model) */
    mxArray *mxRes = mxCreateDoubleMatrix((mwSize)m, 1, mxREAL);
    double *rout = mxGetPr(mxRes);
    for (size_t i = 0; i < m; ++i) rout[i] = gsl_vector_get(F, i);
    mxSetField(plhs[3], 0, "data", mxRes);

    /* pval: not computed here; set NaN */
    mxArray *mxP = mxCreateDoubleScalar(mxGetNaN());
    mxSetField(plhs[3], 0, "pval", mxP);

    /* mean, std, RSS of residuals */
    double sum = 0.0, sum2 = 0.0;
    for (size_t i = 0; i < m; ++i) {
        double r = rout[i];
        sum  += r;
        sum2 += r*r;
    }
    double mean = (m>0) ? sum / (double)m : 0.0;
    double var  = (m>1) ? (sum2 - (sum*sum)/(double)m) / (double)(m-1) : 0.0;
    double sdev = (var > 0.0) ? sqrt(var) : 0.0;

    mxSetField(plhs[3], 0, "mean", mxCreateDoubleScalar(mean));
    mxSetField(plhs[3], 0, "std",  mxCreateDoubleScalar(sdev));
    mxSetField(plhs[3], 0, "RSS",  mxCreateDoubleScalar(sum2));

    /* ---- 4) J (m x np) ---- */
    plhs[4] = mxCreateDoubleMatrix((mwSize)m, (mwSize)np, mxREAL);
    double *Jout = mxGetPr(plhs[4]);
    gsl_matrix *Jcur = gsl_matrix_alloc(m, np);
    gsl_multifit_fdfsolver_jac((gsl_multifit_fdfsolver *)s, Jcur);
    for (size_t r = 0; r < m; ++r) {
        for (size_t c = 0; c < np; ++c) {
            Jout[r + c*m] = gsl_matrix_get(Jcur, r, c);
        }
    }
    gsl_matrix_free(Jcur);

    /* warn if convergence not reached (but return best-so-far) 
    if (status != GSL_SUCCESS) {
        mexWarnMsgIdAndTxt("fitGaussian2D:NoConvergence",
                           "GSL did not reach requested tolerance: %s. Returning best-so-far.",
                           gsl_strerror(status));
    }*/

    /* free covariance copy */
    gsl_matrix_free(cov_np);
}

/* -------------------------- MEX ENTRY --------------------------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3 || nrhs > 4)
        mexErrMsgIdAndTxt("fitGaussian2D:nrhs",
                          "Usage: [prmVect,prmStd,C,res,J] = fitGaussian2D(data, prmVect, mode, [options])");

    /* 1) image */
    const mxArray *mxImg = prhs[0];

    /* 2) prmVect */
    if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 5)
        mexErrMsgIdAndTxt("fitGaussian2D:prmVect",
                          "prmVect must be a 1x5 real double vector [xp yp A s c].");
    const double *pvin = mxGetPr(prhs[1]);

    /* 3) mode */
    const mxArray *mxMode = prhs[2];

    /* 4) options (optional) */
    const mxArray *mxOpt  = (nrhs >= 4) ? prhs[3] : NULL;

    /* set up FitData */
    FitData D;
    extract_valid_pixels(mxImg, &D);
    parse_mode(mxMode, &D);
    parse_options(mxOpt, &D);

    for (int i = 0; i < 5; ++i) D.prmVect[i] = pvin[i];
    if (D.m == 0)
        mexErrMsgIdAndTxt("fitGaussian2D:noData", "No valid (non-NaN) pixels in window.");

    /* GSL setup */
    gsl_set_error_handler_off(); /* avoid abort() on numerical issues */

    gsl_multifit_function_fdf fdf;
    fdf.f      = gaussian_f;
    fdf.df     = gaussian_df;
    fdf.fdf    = gaussian_fdf;
    fdf.n      = D.m;              /* # residuals */
    fdf.p      = D.np;             /* # params */
    fdf.params = &D;

    /* initial x = selected subset from prmVect */
    gsl_vector *x = gsl_vector_alloc(D.np);
    for (size_t k = 0; k < D.np; ++k) {
        gsl_vector_set(x, k, D.prmVect[D.estIdx[k]]);
    }

    /* choose solver */
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc(T, fdf.n, fdf.p);

    int status = gsl_multifit_fdfsolver_set(s, &fdf, x);
    if (status != GSL_SUCCESS) {
        gsl_multifit_fdfsolver_free(s);
        gsl_vector_free(x);
        mxFree(D.X); mxFree(D.Y); mxFree(D.I);
        mexErrMsgIdAndTxt("fitGaussian2D:initFail",
                          "Failed to initialize GSL solver: %s", gsl_strerror(status));
    }

    /* iterate */
    size_t iter = 0;
    
    /* We'll compute "dx" ourselves: keep previous x and subtract */
    gsl_vector *x_prev = gsl_vector_alloc(D.np);
    gsl_vector *dx_tmp = gsl_vector_alloc(D.np);
    
    /* initialize x_prev = initial x (what we passed to set()) */
    {
        const gsl_vector *xpos0 = gsl_multifit_fdfsolver_position(s);
        gsl_vector_memcpy(x_prev, xpos0);
    }
    
    do {
        status = gsl_multifit_fdfsolver_iterate(s);
        if (status && status != GSL_ENOPROG && status != GSL_EBADFUNC) {
            /* unexpected error: break, return best-so-far */
            break;
        }
    
        /* current position */
        const gsl_vector *xpos = gsl_multifit_fdfsolver_position(s);
    
        /* dx_tmp = xpos - x_prev */
        gsl_vector_memcpy(dx_tmp, xpos);
        gsl_blas_daxpy(-1.0, x_prev, dx_tmp);
    
        /* delta test using our dx and current position */
        int dstat = gsl_multifit_test_delta(dx_tmp, xpos, D.eAbs, D.eRel);
        if (dstat == GSL_SUCCESS) { status = GSL_SUCCESS; break; }
    
        /* update x_prev = xpos */
        gsl_vector_memcpy(x_prev, xpos);
    
        ++iter;
    } while (iter < D.maxIter);
    
    /* free temporaries used for delta test */
    gsl_vector_free(x_prev);
    gsl_vector_free(dx_tmp);
    
    /* copy best-so-far back into D.prmVect (for completeness; MATLAB output built below) */
    const gsl_vector *xbest = gsl_multifit_fdfsolver_position(s);
    for (size_t k = 0; k < D.np; ++k) {
        D.prmVect[D.estIdx[k]] = gsl_vector_get(xbest, k);
    }
    if (D.prmVect[3] < 0.0) D.prmVect[3] = fabs(D.prmVect[3]);

    /* build outputs (also computes covariance, residual stats, and J) */
    if (nlhs < 5)
        mexWarnMsgIdAndTxt("fitGaussian2D:nlhs",
                           "Function returns 5 outputs; caller requested %d. Truncation will occur.", nlhs);

    /* Ensure plhs has 5 slots to fill safely; MATLAB ignores extras */
    mxArray *tmp[5] = {NULL,NULL,NULL,NULL,NULL};
    copy_outputs_and_cleanup(&D, s, status, tmp);
    for (int i = 0; i < ((nlhs<5)?nlhs:5); ++i) plhs[i] = tmp[i];

    /* cleanup */
    gsl_multifit_fdfsolver_free(s);
    gsl_vector_free(x);
    mxFree(D.X); mxFree(D.Y); mxFree(D.I);
}