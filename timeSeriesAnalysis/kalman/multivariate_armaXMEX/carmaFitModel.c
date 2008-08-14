/*
 * File: carmaFitModel.c
 * ---------------------
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

  /* Initialize variables and pointers */  
    
  double *TOPO, *TOPOp;
  double *maPARAMS, *maPARAMSp;

  double edgeLen, fTol;
  double *paramV, *vertexLikelihoods;
  double **p;

  double *TOPOfit, *TOPOfitp, *maPARAMSfit, *maPARAMSfitp;

  int i, j, k, nRows, nCols, nDims;
  int nNodes, nParams, nEvals, nLags, arOrderMax, maOrderMax;

  struct probCDT prob;

  mxArray *pm;      /* generic mxArray pointer */

  const mwSize *trajDim, *topoDim, *topoBinDim;
  mwSize topoD[3];



  prob.numMissing = 0;



  /* Retrieve model parameters and time series from input. */


  /*Verify input is a structure */
  if (!mxIsStruct(prhs[0]))
    mxErrMsgTxt("First input must be a structure.");  
  
  /* Retrieve the topology matrix containing AR and X parameters and 
   * determine its size */
  pm = mxGetField(prhs[0], 0, "TOPO");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve TOPO from input.");

  TOPO = mxGetPr(pm);

  nDims = (int) mxGetNumberOfDimensions(pm);

  topoDim = mxGetDimensions(pm);

  prob.arOrderMax = arOrderMax = topoDim[0] - 1;
  nLags = topoDim[0];
  prob.nNodes = nNodes = topoDim[1];

  if ((nNodes > 1) && (nDims != 3))
    mxErrMsgTxt("If nNodes > 1, input matrix must be 3D.");

  if (nDims == 3)
    if (topoDim[2] != nNodes)
      mxErrMsgTxt("Input topology matrix has wrong dimensions.");

  /* Retrieve the MA parameters from input structure and determine 
   * and determine their size */
  
  pm = mxGetField(prhs[0], 0, "maPARAMS");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve maPARAMS from input.");

  maPARAMS = mxGetPr(pm);

  nRows = (int) mxGetM(pm);
  nCols = (int) mxGetN(pm);

  prob.maOrderMax = maOrderMax = nRows;

  if (nCols != nNodes)
    mxErrMsgTxt("Second dimension of maPARAMS shouldbe equal to equal nNodes.");



  /* Get the pointer for the binary moving average parameter matrix */

  pm = mxGetField(prhs[0], 0, "maBIN");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve maBIN from input.");

  prob.maBIN = (int *) mxGetPr(pm);

  nRows = (int) mxGetM(pm);
  if (nRows != maOrderMax)
    mxErrMsgTxt("maBIN must be same size as maPARAMS.");

  nCols = (int) mxGetN(pm);
  if (nCols != nNodes)
    mxErrMsgTxt("2nd dimension of maBIN should be equal to nNodes.");



  /* Get the pointer for input time series */

  pm = mxGetField(prhs[0], 0, "TRAJ");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve TRAJ from input.");

  prob.TRAJ = mxGetPr(pm);

  nDims = mxGetNumberOfDimensions(pm);
		
  if (nDims != 3)
    mxErrMsgTxt("Input time series matrix should be 3D. ", 
		"Put zeros as observational error if none exists.");



  /* Determine the time series length */

  trajDim = mxGetDimensions(pm);
  prob.trajLength = trajDim[0];

  if ((trajDim[1] != nNodes) || (trajDim[2] != 2))
    mxErrMsgTxt("Size of input time series is incorrect.");



  /* Get pointer for binary topology matrix */

  pm = mxGetField(prhs[0], 0, "topoBIN");
  if (pm == NULL) mxErrMsgTxt("Cannot retrieve topoBIN from first input.");

  prob.topoBIN = (int *) mxGetPr(pm);

  nDims = (int) mxGetNumberOfDimensions(pm);

  if ((nNodes > 1) && (nDims != 3))
    mxErrMsgTxt("if nNodes > 1, topoBIN should be 3D.");



  /* Binary topology should be equal to topology */

  topoBinDim = mxGetDimensions(pm);

  /* Check the size of time series */

  if ((topoBinDim[0] != nLags) || (topoBinDim[1] != nNodes))
    mxErrMsgTxt("topoBIN should be equal to TOPO dimensionally.");

  if ((nDims > 2) && (topoBinDim[2] != nNodes))
    mxErrMsgTxt("topoBIN should be equal to TOPO dimensionally.");

 

  /* Convert ARMA parameters to partial parameters for the minimizer */

  TOPOp = (double *) malloc(sizeof(double) * nLags * nNodes * nNodes);



  for (i = 0; i < nNodes; i++){
    for (j = 0; j < nNodes; j++){
   
      /* Along the diagonal are AR parameters */
    
      if (i == j) {
    
	/* No AR parameters at lag 0 */
    
	*(TOPOp + i * nLags + j * nLags * nNodes) = 0;

	inverseLevinsonDurbinExpoAR( (TOPO + i * nLags + j * nLags * nNodes + 1), arOrderMax,
				     (TOPOp + i * nLags + j * nLags * nNodes + 1) );

      } else {

	/* Copy directly the X-parameters */

	for (k = 0; k < nLags; k++)
	  *(TOPOp + i * nLags + j * nLags * nNodes + k) = *(TOPO + i * nLags + j * nLags * nNodes + k);

      }

    }
  }
 

  maPARAMSp = (double *) malloc(sizeof(double) * nNodes * maOrderMax);

  for (i = 0; i < nNodes; i++)
    inverseLevinsonDurbinExpoMA(maPARAMS + i * maOrderMax, maOrderMax,
				maPARAMSp + i * maOrderMax);



  /* Assemble the parameters into a vector for passing to minimizer */

  vectorFromParams(TOPOp, maPARAMSp, &prob, &paramV);

  nParams = prob.numParams;

  prob.wnVariance = 1;



  p = matrix(1, nParams + 1, 1, nParams);

  /* Allocate vector for function values at simplex vertices */

  vertexLikelihoods = vector(1, nParams + 1);

  edgeLen = 3.0E-1;

  createSimplex(paramV, nParams, edgeLen, p);



  /* Now evaluate the likelihood at each vertex */

  for (i = 1; i <= nParams + 1; i++)
    vertexLikelihoods[i] = carmaNegLnLikelihood(p[i], &prob);


  /* Call the minimizer for maximum likelihood estimation of model params */
  nEvals = 0;
  fTol = 1.0E-8;

  if (amoeba(p, vertexLikelihoods, nParams, fTol, carmaNegLnLikelihood, &nEvals, &prob)){
    printf("\ncarmaFitModel: Minimizer succeeds.\n");
    plhs[3] = mxCreateDoubleScalar(1);
  } else{
    printf("\ncarmaFitModel: Minimizer fails.\n");
    printf("\nMaximum number of evaluations is exceeded.\n");
    plhs[3] = mxCreateDoubleScalar(0);
  }

  /* Set up outputs for results of minimization */

  /* Dimensions of the TOPO matrix */

  topoD[0] = nLags;
  topoD[1] = topoD[2] = nNodes;


  /* Allocate output topology matrix */

  plhs[0] = mxCreateNumericArray(3, topoD, mxDOUBLE_CLASS, mxREAL);

  /* output MA parameters */

  plhs[1] = mxCreateDoubleMatrix(maOrderMax, nNodes, mxREAL);

  /* BIC */

  plhs[2] = mxCreateDoubleScalar(vertexLikelihoods[1] + nParams * log(prob.trajLength - prob.numMissing));



  TOPOfit = mxGetPr(plhs[0]);
  TOPOfitp = (double *) malloc(sizeof(double) * nLags * nNodes * nNodes);

  maPARAMSfit = mxGetPr(plhs[1]);
  maPARAMSfitp = (double *) malloc(sizeof(double) * nNodes * maOrderMax);



  /* Extract partial parameters from simplex */

  paramsFromVector(&TOPOfitp, nNodes, arOrderMax,
		   &maPARAMSfitp, maOrderMax,
		   prob.topoBIN, prob.maBIN,
		   p[1] + 1, nParams);

  /* Convert these to regular parameters for return */

  for (i = 0; i < nNodes; i++){
    for (j = 0; j < nNodes; j++){
	  
	/* Along the diagonal the parameters are AR */

	if (i == j) {
	    
	  /* No AR parameters at lag = 0 */

	  *(TOPOfit + i * nLags + j * nLags * nNodes) = 0;
	    
	  levinsonDurbinExpoAR( (TOPOfitp + i * nLags + j * nLags * nNodes + 1), arOrderMax,
				(TOPOfit + i * nLags + j * nLags * nNodes + 1) );

	} else {

	  /* Copy directly the X-parameters */

	  for (k = 0; k < nLags; k++)
	    *(TOPOfit + i * nLags + j * nLags * nNodes + k) = 
	      *(TOPOfitp + i * nLags + j * nLags * nNodes + k);
	}
  
    }
  } 
  /* Convert partial MA parameters to full for return */
  for (i = 0; i < nNodes; i++)
    levinsonDurbinExpoMA(maPARAMSfitp + i * maOrderMax, maOrderMax,
			 maPARAMSfit + i * maOrderMax);


  /* Free allocated memory */
  free(paramV);
  free_vector(vertexLikelihoods, 1, nParams + 1);
  free_matrix(p, 1, nParams + 1, 1, nParams);
}
