/*
 * File: carmaNegLnLikelihood
 * --------------------------
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "prob.h"
#include "prob.c"
#include "carmaFitModel.h"

#if !defined(MAX)

#define	MAX(A, B)	((A) > (B) ? (A) : (B))

#endif

double carmaNegLnLikelihood(double *paramNR, void *d)
{

  probADT prob = (probADT) d;
  double *paramV = paramNR + 1;  /* Convert Numerical Recipes array to C array */

  double *TOPO, *TOPOp;

  double *maPARAMS, *maPARAMSp;

  double *innovations, *innovationVars, *wnV;

  double likelihood, sum1, sum2;

  int i, j, k, n;

  int trajLen, nMissing;
  
  double *traj;

  int nMovies = prob->nMovies;
  int nNodes = prob->nNodes;
  int nParams = prob->nParams;
  int arOrderMax = prob->arOrderMax;
  int maOrderMax = prob->maOrderMax;
  int nLags = arOrderMax + 1;
  int *topoBIN = prob->topoBIN;
  int *maBIN = prob->maBIN;

  double wnVariance = prob->wnVariance;



  /* Allocate memory for full and partial parameters */

  TOPO = (double *) malloc(sizeof(double) * nNodes * nNodes * nLags);
  TOPOp = (double *) malloc(sizeof(double) * nNodes * nNodes * nLags);
  maPARAMS = (double *) malloc(sizeof(double) * nNodes * maOrderMax);
  maPARAMSp = (double *) malloc(sizeof(double) * nNodes * maOrderMax);



  likelihood = 0;
  sum1 = sum2 = 0;

  for (n = 0; n < nMovies; n++) {

    traj = prob->data[n].traj;
    trajLen = prob->data[n].trajLen;
    nMissing = prob->data[n].nMissing;

    /* Retrieve partial parameters from parameter vector */

    paramsFromVector(&TOPOp, nNodes, arOrderMax,
		     &maPARAMSp, maOrderMax,
		     topoBIN, maBIN,
		     paramV, nParams);


    /* Convert partial parameters to full parameters */

    for (i = 0; i < nNodes; i++){
      for (j = 0; j < nNodes; j++){

	/* Along the diagonal are AR parameters */

	if (i == j) {

	  /* No AR parameters at lag 0 */

	  *(TOPO + i * nLags + j * nNodes * nLags) = 0;
      
	  /* Convert partial AR to full parameters for innovations calculation */

	  levinsonDurbinExpoAR((TOPOp + i * nLags + j * nNodes * nLags + 1),
			       arOrderMax,
			       (TOPO + i * nLags + j * nNodes * nLags + 1));
	} else {

	  /* Copy directly the X-parameters */

	  for (k = 0; k < nLags; k++)
	    *(TOPO + i * nLags + j * nNodes * nLags + k) = 
		*(TOPOp + i * nLags + j * nNodes * nLags + k);

	}
	
      }
    }



    /* Convert partial MA parameters to full for innovations calculation */

    for (i = 0; i < nNodes; i++)
      levinsonDurbinExpoMA(maPARAMSp + i * maOrderMax,
			   maOrderMax,
			   maPARAMS + i * maOrderMax);


    /* Allocate memory for innovations, their variances and the white noise vector */

    innovations = (double *) malloc(sizeof(double) * trajLen);

    innovationVars = (double *) malloc(sizeof(double) * trajLen);

    wnV = (double *) malloc(sizeof(double) * trajLen);



    /* Calculate the innovations for each node, and sum the -2lnlikelihoods */

    for (i = 0; i < nNodes; i++){
      carmaCalcKalmanInnov(traj, trajLen,
			   TOPO, topoBIN,
			   arOrderMax, nNodes,
			   maPARAMS, maOrderMax,
			   wnVariance, i,
			   &innovations, &innovationVars, &wnV,
			   &sum1, &sum2, &nMissing);
    }

    prob->data[n].nMissing = MAX(prob->data[n].nMissing, nMissing);

    free(innovations);
    free(innovationVars);
    free(wnV);
  }



  free(TOPO);
  free(TOPOp);
  free(maPARAMS);
  free(maPARAMSp);

  likelihood = sum1 + (trajLen - nMissing) * log(sum2);
  
  return likelihood;
}
