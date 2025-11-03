/* fitGaussian2D.c  — MEX wrapper for 2D Gaussian fit with GSL LM (isotropic sigma)
 *
 * p = [A, x0, y0, sigma, b], model:
 *   y(x,y) = A * exp(-((x-x0)^2 + (y-y0)^2) / (2*sigma^2)) + b
 *
 * [params, resnorm, status, info] = fitGaussian2D(I, init, [mask], [opts])
 *   I:     MxN double
 *   init:  double[5] initial guess [A, x0, y0, sigma, b]
 *   mask:  optional logical MxN, same size as I (default: all true)
 *   opts:  struct with fields:
 *            maxIter (default 200)
 *            eAbs    (default 1e-8)
 *            eRel    (default 1e-8)
 *
 * Outputs:
 *   params : double[5] best-fit parameters
 *   resnorm: double, sum of squared residuals
 *   status : int (GSL status code; 0 = GSL_SUCCESS)
 *   info   : struct with fields: iter, gradnorm, m, n, message
 *
 * Notes:
 *  - Uses only GSL 2.x public API; does not access solver internals.
 *  - Robust stopping via test_delta (on x) and test_gradient (on J^T f).
 *  - gsl_set_error_handler_off() prevents abort; errors become return codes.
 *
 * Compile (macOS + Homebrew GSL, MATLAB R2024b; adjust paths as needed)
   clear mex
      mex -v ...
      -I/opt/homebrew/opt/gsl/include ...
      -L/opt/homebrew/opt/gsl/lib ...
      CFLAGS='$CFLAGS -std=c11' ...
      LDFLAGS='$LDFLAGS -Wl,-rpath,/opt/homebrew/opt/gsl/lib' ...
        fitGaussian2D.c -lgsl -lgslcblas -lm
 */
#include "mex.h"
#include "matrix.h"

#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* ------------------ Helpers to avoid touching solver internals ------------------ */
#ifndef SOLVER_M
#  define SOLVER_M(s) ((size_t)((s)->f ? (s)->f->size : 0))   /* # residuals */
#endif
#ifndef SOLVER_N
#  define SOLVER_N(s) ((size_t)((s)->x ? (s)->x->size : 0))   /* # params */
#endif

static inline const gsl_vector* solver_F(const gsl_multifit_fdfsolver *s) {
    return gsl_multifit_fdfsolver_residual(s); /* const pointer */
}

static inline int solver_fill_J(const gsl_multifit_fdfsolver *s, gsl_matrix *Jout) {
    return gsl_multifit_fdfsolver_jac(s, Jout);
}
/* ------------------------------------------------------------------------------- */

/* Parameter indices */
enum { P_A = 0, P_X0 = 1, P_Y0 = 2, P_SIG = 3, P_B = 4, P_COUNT = 5 };

/* Data passed to callbacks */
typedef struct {
    const double *I;       /* image data, column-major (N columns, M rows) */
    mwSize M, N;           /* image size */
    size_t m;              /* number of residuals (valid pixels) */
    int *vidx;             /* linear indices (0..M*N-1) of valid pixels */
    double eAbs, eRel;     /* tolerances */
    size_t maxIter;
} FitArgs;

/* Gaussian model at coordinate (x,y) */
static inline double gauss2d_iso(double A, double x0, double y0, double sig, double b,
                                 double x, double y)
{
    const double dx = x - x0;
    const double dy = y - y0;
    const double s2 = sig * sig;
    return A * exp(-(dx*dx + dy*dy) / (2.0 * s2)) + b;
}

/* f: residuals r_k = I(i_k, j_k) - model(x_k, y_k) */
static int gaussian_f(const gsl_vector *x, void *params, gsl_vector *f)
{
    FitArgs *A = (FitArgs*)params;
    const double Aamp  = gsl_vector_get(x, P_A);
    const double x0    = gsl_vector_get(x, P_X0);
    const double y0    = gsl_vector_get(x, P_Y0);
    const double sigma = fabs(gsl_vector_get(x, P_SIG)) + 1e-12; /* keep > 0 */
    const double b     = gsl_vector_get(x, P_B);

    const size_t m = A->m;
    const mwSize M = A->M;
    const mwSize N = A->N;

    for (size_t k = 0; k < m; ++k) {
        int L = A->vidx[k];            /* linear index in column-major */
        int i = L % M;                 /* row    (0..M-1) => y */
        int j = L / M;                 /* column (0..N-1) => x */
        double yobs = A->I[L];
        double yhat = gauss2d_iso(Aamp, x0, y0, sigma, b, (double)j, (double)i);
        gsl_vector_set(f, k, yobs - yhat);
    }
    return GSL_SUCCESS;
}

/* df: Jacobian J(k, p) = d r_k / d p */
static int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
    FitArgs *A = (FitArgs*)params;
    const double Aamp  = gsl_vector_get(x, P_A);
    const double x0    = gsl_vector_get(x, P_X0);
    const double y0    = gsl_vector_get(x, P_Y0);
    const double sigma = fabs(gsl_vector_get(x, P_SIG)) + 1e-12; /* > 0 */
    const double b     = gsl_vector_get(x, P_B);
    (void)b; /* not needed in derivatives except db */

    const double s2 = sigma * sigma;

    const size_t m = A->m;
    const mwSize M = A->M;

    for (size_t k = 0; k < m; ++k) {
        int L  = A->vidx[k];
        int i  = L % M;     /* y */
        int j  = L / M;     /* x */
        double dx = (double)j - x0;
        double dy = (double)i - y0;
        double E  = exp(-(dx*dx + dy*dy) / (2.0 * s2));  /* exp part */

        /* residual r = I - (A*E + b)
         * dr/dA    = -E
         * dr/dx0   = -A * E * (dx / s2)
         * dr/dy0   = -A * E * (dy / s2)
         * dr/dsig  = -A * E * ((dx^2 + dy^2)/ (sigma^3))
         * dr/db    = -1
         */
        double dA   = -E;
        double dx0p = -Aamp * E * ( dx / s2 );
        double dy0p = -Aamp * E * ( dy / s2 );
        double dsig = -Aamp * E * ( (dx*dx + dy*dy) / (sigma*sigma*sigma) );
        double db   = -1.0;

        gsl_matrix_set(J, k, P_A,   dA);
        gsl_matrix_set(J, k, P_X0,  dx0p);
        gsl_matrix_set(J, k, P_Y0,  dy0p);
        gsl_matrix_set(J, k, P_SIG, dsig);
        gsl_matrix_set(J, k, P_B,   db);
    }
    return GSL_SUCCESS;
}

/* fdf: both f and df */
static int gaussian_fdf(const gsl_vector *x, void *params,
                        gsl_vector *f, gsl_matrix *J)
{
    int s = gaussian_f(x, params, f);
    if (s != GSL_SUCCESS) return s;
    return gaussian_df(x, params, J);
}

/* Build mask index list (vidx). If mask==NULL, use all pixels */
static void build_valid_index(const mxArray *mask,
                              mwSize M, mwSize N,
                              int **vidx_out, size_t *m_out)
{
    size_t cap = (size_t)M * (size_t)N;
    int *vidx = (int*)mxCalloc(cap, sizeof(int));
    size_t m  = 0;

    if (mask && mxIsLogical(mask) &&
        mxGetM(mask)==M && mxGetN(mask)==N)
    {
        mxLogical *mk = mxGetLogicals(mask);
        for (mwSize j=0; j<N; ++j) {
            for (mwSize i=0; i<M; ++i) {
                size_t L = (size_t)i + (size_t)M*(size_t)j;
                if (mk[L]) vidx[m++] = (int)L;
            }
        }
    } else {
        /* all valid */
        for (mwSize j=0; j<N; ++j) {
            for (mwSize i=0; i<M; ++i) {
                size_t L = (size_t)i + (size_t)M*(size_t)j;
                vidx[m++] = (int)L;
            }
        }
    }
    *vidx_out = vidx;
    *m_out    = m;
}

/* Create MATLAB info struct */
static mxArray* make_info_struct(size_t iter, double gradnorm, size_t m, size_t n, int status)
{
    const char *fn[] = {"iter","gradnorm","m","n","message"};
    mxArray *S = mxCreateStructMatrix(1,1,5,fn);
    mxSetField(S, 0, "iter",     mxCreateDoubleScalar((double)iter));
    mxSetField(S, 0, "gradnorm", mxCreateDoubleScalar(gradnorm));
    mxSetField(S, 0, "m",        mxCreateDoubleScalar((double)m));
    mxSetField(S, 0, "n",        mxCreateDoubleScalar((double)n));
    mxSetField(S, 0, "message",  mxCreateString(gsl_strerror(status)));
    return S;
}

/* --------------------------------- MEX entry ---------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2) {
        mexErrMsgIdAndTxt("fitGaussian2D:Args",
            "Usage: [p,resnorm,status,info] = fitGaussian2D(I, init, [mask], [opts])");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0])!=2) {
        mexErrMsgIdAndTxt("fitGaussian2D:Image","I must be real double MxN.");
    }
    const mxArray *Iarr = prhs[0];
    const double *I = mxGetPr(Iarr);
    const mwSize M = mxGetM(Iarr);
    const mwSize N = mxGetN(Iarr);

    const mxArray *initArr = prhs[1];
    if (!mxIsDouble(initArr) || mxGetNumberOfElements(initArr) != P_COUNT) {
        mexErrMsgIdAndTxt("fitGaussian2D:Init","init must be double[5] = [A x0 y0 sigma b].");
    }
    const double *init = mxGetPr(initArr);

    const mxArray *mask = NULL;
    const mxArray *opts = NULL;
    if (nrhs >= 3) mask = prhs[2];
    if (nrhs >= 4) opts = prhs[3];

    /* Defaults */
    size_t maxIter = 200;
    double eAbs = 1e-8;
    double eRel = 1e-8;

    if (opts && mxIsStruct(opts)) {
        mxArray *v;
        if ((v = mxGetField(opts,0,"maxIter")) && mxIsDouble(v) && mxGetNumberOfElements(v)==1) {
            double t = mxGetScalar(v);
            if (t > 0) maxIter = (size_t)(t + 0.5);
        }
        if ((v = mxGetField(opts,0,"eAbs")) && mxIsDouble(v) && mxGetNumberOfElements(v)==1) {
            double t = mxGetScalar(v);
            if (t > 0) eAbs = t;
        }
        if ((v = mxGetField(opts,0,"eRel")) && mxIsDouble(v) && mxGetNumberOfElements(v)==1) {
            double t = mxGetScalar(v);
            if (t > 0) eRel = t;
        }
    }

    /* Build valid indices */
    int *vidx = NULL; size_t m = 0;
    build_valid_index(mask, M, N, &vidx, &m);
    if (m < P_COUNT) {
        mxFree(vidx);
        mexErrMsgIdAndTxt("fitGaussian2D:Mask","Not enough valid pixels for fitting.");
    }

    /* Prepare args */
    FitArgs A;
    A.I = I;
    A.M = M;
    A.N = N;
    A.m = m;
    A.vidx = vidx;
    A.eAbs = eAbs;
    A.eRel = eRel;
    A.maxIter = maxIter;

    /* GSL config */
    gsl_set_error_handler_off();

    gsl_multifit_function_fdf fcfg;
    fcfg.f      = gaussian_f;
    fcfg.df     = gaussian_df;
    fcfg.fdf    = gaussian_fdf;
    fcfg.n      = (unsigned int)m;      /* # residuals */
    fcfg.p      = (unsigned int)P_COUNT;/* # params    */
    fcfg.params = (void*)&A;

    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmder;
    gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc(T, m, P_COUNT);
    if (!s) {
        mxFree(vidx);
        mexErrMsgIdAndTxt("fitGaussian2D:GSL","Failed to allocate solver.");
    }

    gsl_vector *x = gsl_vector_alloc(P_COUNT);
    for (int i=0;i<P_COUNT;++i) gsl_vector_set(x, i, init[i]);
    /* basic clamping on sigma initial guess */
    if (gsl_vector_get(x, P_SIG) <= 0) gsl_vector_set(x, P_SIG, fabs(gsl_vector_get(x, P_SIG))+1.0);

    int gerr = gsl_multifit_fdfsolver_set(s, &fcfg, x);
    if (gerr != GSL_SUCCESS) {
        gsl_vector_free(x);
        gsl_multifit_fdfsolver_free(s);
        mxFree(vidx);
        mexErrMsgIdAndTxt("fitGaussian2D:GSLSet", "gsl_multifit_fdfsolver_set: %s", gsl_strerror(gerr));
    }

    /* Iterate with robust tests: Δx and gradient */
    size_t iter = 0;
    int status = GSL_CONTINUE;
    int status_grad = GSL_CONTINUE;
    double gradnorm_out = mxGetNaN();

    while ((status == GSL_CONTINUE || status_grad == GSL_CONTINUE) && iter < maxIter) {
        ++iter;

        int step = gsl_multifit_fdfsolver_iterate(s);
        if (step && step != GSL_ENOPROG && step != GSL_EBADFUNC) {
            status = step; /* unexpected failure */
            break;
        }

        status = gsl_multifit_test_delta(s->dx, s->x, eAbs, eRel);

        /* gradient test: ||J^T f|| */
        {
            const size_t mm = SOLVER_M(s);
            const size_t nn = SOLVER_N(s);
            gsl_matrix *Jtmp = gsl_matrix_alloc(mm, nn);
            gsl_vector *Ftmp = gsl_vector_alloc(mm);
            gsl_vector *grad = gsl_vector_alloc(nn);

            solver_fill_J(s, Jtmp);
            gsl_vector_memcpy(Ftmp, solver_F(s));
            gsl_multifit_gradient(Jtmp, Ftmp, grad);

            double gnorm = 0.0;
            gsl_blas_dnrm2(grad, &gnorm);
            gradnorm_out = gnorm;

            status_grad = gsl_multifit_test_gradient(grad, eAbs);

            gsl_vector_free(grad);
            gsl_vector_free(Ftmp);
            gsl_matrix_free(Jtmp);
        }
    }

    /* Prepare outputs */
    plhs[0] = mxCreateDoubleMatrix(1, P_COUNT, mxREAL);
    double *pout = mxGetPr(plhs[0]);
    for (int i=0;i<P_COUNT;++i) {
        pout[i] = gsl_vector_get(s->x, i);
    }
    if (pout[P_SIG] <= 0) pout[P_SIG] = fabs(pout[P_SIG]);

    /* resnorm = ||f||^2 */
    double resnorm = 0.0;
    {
        const gsl_vector *F = solver_F(s);
        double nrm = 0.0;
        gsl_blas_dnrm2(F, &nrm);
        resnorm = nrm*nrm;
    }
    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleScalar(resnorm);
    }
    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleScalar((double)status);
    }
    if (nlhs >= 4) {
        plhs[3] = make_info_struct(iter, gradnorm_out, SOLVER_M(s), SOLVER_N(s), status);
    }

    /* Warn but do not error on non-convergence */
    if (status != GSL_SUCCESS && status != GSL_CONTINUE) {
        mexWarnMsgIdAndTxt("fitGaussian2D:NoConvergence",
            "GSL did not reach tolerance: %s. Returning best-so-far.",
            gsl_strerror(status));
    }

    /* Cleanup */
    gsl_vector_free(x);
    gsl_multifit_fdfsolver_free(s);
    mxFree(vidx);
}