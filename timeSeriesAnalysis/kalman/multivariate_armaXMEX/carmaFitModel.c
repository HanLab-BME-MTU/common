/*
 * File: carmaFitModel.c
 * ---------------------
 * The mexFunction is designed to be called in armaxFitKalman.m, if
 * option 'nl' (Numerical Recipes local minimizer) is chosen.
 *
 * Hunter Elliott, Pei-hsin Hsu, September 2008
 */

#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "prob.c"
#include "./Utilities/Carma_Utilities/carmaUtils.h"
#include "./Utilities/Matrix_Operations/matOps.h"
#include "./Minimizer/amoeba.h"
#include "./Minimizer/nrutil.h"
#include "carmaFitModel.h"



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{  
    
  double *TOPOp;
  double *maPARAMSp;

  double edgeLen, fTol;
  double *paramV, *vertexLikelihoods;
  double **p;

  double *TOPOfit, *TOPOfitp, *maPARAMSfit, *maPARAMSfitp;

  int i, j, k, nRows, nCols, nDims;
  int nMovies, nNodes, nParams, nEvals, nLags, nPoints, arOrderMax, maOrderMax;

  struct probCDT prob;

  mxArray *pm, *pt; /* pointer to mxArray, pointer to trajOut */

  const mwSize *trajOutDim, *trajDim, *topoDim, *topoBinDim;

  mwSize topoD[3];



  /* Verify input is a structure */

  if (!mxIsStruct(prhs[0]))
    mxErrMsgTxt("First input must be a structure.");  



  /* Retrieve topology matrix containing AR and X parameters */

  pm = mxGetField(prhs[0], 0, "TOPOp");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve TOPOp from input.");

  TOPOp = mxGetPr(pm);

  nDims = (int) mxGetNumberOfDimensions(pm);

  topoDim = mxGetDimensions(pm);

  prob.arOrderMax = arOrderMax = topoDim[0] - 1;
  nLags = topoDim[0];
  prob.nNodes = nNodes = topoDim[1];

  /* mexPrintf("\nprob.nNodes = %d\n", prob.nNodes); */

  if ((nNodes > 1) && (nDims != 3))
    mxErrMsgTxt("If nNodes > 1, input matrix must be 3D.");

  if (nDims == 3)
    if (topoDim[2] != nNodes)
      mxErrMsgTxt("Input topology matrix has wrong dimensions.");



  /* Retrieve MA parameters */
  
  pm = mxGetField(prhs[0], 0, "maPARAMSp");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve maPARAMSp from input.");

  maPARAMSp = mxGetPr(pm);

  nRows = (int) mxGetM(pm);
  nCols = (int) mxGetN(pm);

  prob.maOrderMax = maOrderMax = nRows;

  if (nCols != nNodes)
    mxErrMsgTxt("Second dimension of maPARAMSp should be equal to nNodes.");



  /* Retrieve binary matrix of MA parameters */

  pm = mxGetField(prhs[0], 0, "maBIN");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve maBIN from input.");

  prob.maBIN = (int *) mxGetPr(pm);

  nRows = (int) mxGetM(pm);
  if (nRows != maOrderMax)
    mxErrMsgTxt("maBIN must be same size as maPARAMS.");

  nCols = (int) mxGetN(pm);
  if (nCols != nNodes)
    mxErrMsgTxt("2nd dimension of maBIN should be equal to nNodes.");



  /* Retrieve an array of time series */

  pt = mxGetField(prhs[0], 0, "trajOut");
  if (pt == NULL) mxErrMsgTxt("Cannot retrieve trajOut from input.");

  trajOutDim = mxGetDimensions(pt);
  /* mexPrintf("\ntrajOutD = %d x %d\n", trajOutDim[0], trajOutDim[1]); */

  prob.nMovies = nMovies = trajOutDim[1];

  prob.data = (struct movie *) malloc(nMovies * sizeof(struct movie));
  if (prob.data == NULL) mxErrMsgTxt("No memory available for prob.data");

  for (i = 0; i < nMovies; i++){  
 
    pm = mxGetField(pt, i, "observations");
    if (pm == NULL) mxErrMsgTxt("Cannot retrieve observations.");

    prob.data[i].traj = mxGetPr(pm);

    nDims = mxGetNumberOfDimensions(pm);
		
    if (nDims != 3)
      mxErrMsgTxt("Input time series matrix should be 3D. Assign observational error to 0 if none exists.");

    trajDim = mxGetDimensions(pm);
    prob.data[i].trajLen = trajDim[0];

    if ((trajDim[1] != nNodes) || (trajDim[2] != 2))
      mxErrMsgTxt("Size of input time series is incorrect.");

    prob.data[i].nMissing = 0;      /* temp */
  }



  /* Retrieve binary matrix of topology */

  pm = mxGetField(prhs[0], 0, "topoBIN");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve topoBIN from first input.");

  prob.topoBIN = (int *) mxGetPr(pm);

  nDims = (int) mxGetNumberOfDimensions(pm);

  if ((nNodes > 1) && (nDims != 3))
    mxErrMsgTxt("if nNodes > 1, topoBIN should be 3D.");



  /* Binary topology should be equal to topology */

  topoBinDim = mxGetDimensions(pm);

  if ((topoBinDim[0] != nLags) || (topoBinDim[1] != nNodes))
    mxErrMsgTxt("topoBIN should be equal to TOPO dimensionally.");

  if ((nDims > 2) && (topoBinDim[2] != nNodes))
    mxErrMsgTxt("topoBIN should be equal to TOPO dimensionally.");

  /* Retrieve white noise variance if input */

  pm = mxGetField(prhs[0],0,"wnVariance");
  if (pm == NULL){ 
      prob.wnVariance = 1;
  } else {
      prob.wnVariance = *mxGetPr(pm);
  }

 

  /* Assemble the parameters into a vector for passing to minimizer */

  vectorFromParams(TOPOp, maPARAMSp, &prob, &paramV);

  nParams = prob.nParams;


  p = matrix(1, nParams + 1, 1, nParams);

  /* Allocate a vector for function values at simplex vertices */

  vertexLikelihoods = vector(1, nParams + 1);

  edgeLen = 3.0E-1;

  createSimplex(paramV, nParams, edgeLen, p);



  /* Evaluate the likelihood at each vertex */

  for (i = 1; i <= nParams + 1; i++) {
    vertexLikelihoods[i] = carmaNegLnLikelihood(p[i], &prob);
    /*   mexPrintf("\nvertex %d: likelihood = %g\n", i, vertexLikelihoods[i]); */
  }


  /* 
   * Minimizer: amoeba
   * -----------------
   * (1) If minimization is successful, the third output plhs[2] is 1. The variable "proceed"
   *     of armaxCoefKalman.m should be assigned to this value.
   * (2) Originally, plhs[2] is Bayesian Information Criterion (BIC). The calculation of BIC
   *     is commented out to fit this module to armaxCoefKalman.m.
   */

  nEvals = 0;
  fTol = 1.0E-8;

  if (amoeba(p, vertexLikelihoods, nParams, fTol, carmaNegLnLikelihood, &nEvals, &prob)){
    printf("\ncarmaFitModel: Minimizer succeeds.\n");
    plhs[2] = mxCreateDoubleScalar(1);
  } else{
    printf("\ncarmaFitModel: Minimizer fails.");
    printf("\nMaximum number of evaluations is exceeded.\n");
    plhs[2] = mxCreateDoubleScalar(0);
  }



  /* Set up first two outputs */

  /* Dimensions of the TOPO matrix */

  topoD[0] = nLags;
  topoD[1] = topoD[2] = nNodes;


  /* Output topology matrix */

  plhs[0] = mxCreateNumericArray(3, topoD, mxDOUBLE_CLASS, mxREAL);

  /* Output MA parameters */

  plhs[1] = mxCreateDoubleMatrix(maOrderMax, nNodes, mxREAL);



  TOPOfitp = mxGetPr(plhs[0]);
  maPARAMSfitp = mxGetPr(plhs[1]);

  /* TOPOfitp = (double *) malloc(sizeof(double) * nLags * nNodes * nNodes); */
  /* maPARAMSfitp = (double *) malloc(sizeof(double) * nNodes * maOrderMax); */



  /* Extract partial parameters from simplex */

  paramsFromVector(&TOPOfitp, nNodes, arOrderMax,
		   &maPARAMSfitp, maOrderMax,
		   prob.topoBIN, prob.maBIN,
		   p[1] + 1, nParams);

  

  /* Convert partial TOPO parameter to full for return */

  /* 
  for (i = 0; i < nNodes; i++){
    for (j = 0; j < nNodes; j++){

	if (i == j) {

	  *(TOPOfit + i * nLags + j * nLags * nNodes) = 0;
	    
	  levinsonDurbinExpoAR( (TOPOfitp + i * nLags + j * nLags * nNodes + 1), arOrderMax,
				(TOPOfit + i * nLags + j * nLags * nNodes + 1) );

	} else {

	  for (k = 0; k < nLags; k++)
	    *(TOPOfit + i * nLags + j * nLags * nNodes + k) = 
	      *(TOPOfitp + i * nLags + j * nLags * nNodes + k);
	}
  
    }
  }
  */

 
  /* Convert partial MA parameters to full for return */

  /*
  for (i = 0; i < nNodes; i++)
    levinsonDurbinExpoMA(maPARAMSfitp + i * maOrderMax, maOrderMax,
			 maPARAMSfit + i * maOrderMax);
  */


  /*
   * plhs[2]
   * -------
   * plhs[2] contains the Bayesian Information Criterion (BIC).
   * Ref: Jaqaman et al. Biophys. J. 91: 2312-2325 (2006), equation (5)
   */
  /*
  nPoints = 0;

  for (i = 0; i < nMovies; i++)
    nPoints += prob.data[i].trajLen - prob.data[i].nMissing;

  plhs[2] = mxCreateDoubleScalar(vertexLikelihoods[1] + log(nPoints) * nParams);
  */


  free(prob.data);
  free(paramV);
  free_vector(vertexLikelihoods, 1, nParams + 1);
  free_matrix(p, 1, nParams + 1, 1, nParams);
}
