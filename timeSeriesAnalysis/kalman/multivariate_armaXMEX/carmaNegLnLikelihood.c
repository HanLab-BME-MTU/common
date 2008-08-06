#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "prob.h"
#include "prob.c"
#include "carmaFitModel.h"

#if !defined(MAX)

#define	MAX(A, B)	((A) > (B) ? (A) : (B))

#endif

double carmaNegLnLikelihood(double *paramNR, void *data)
			  

{

  probADT prob = (probADT) data;
  double *paramV = paramNR + 1;  /* Convert Numerical Recipes array to C array */

  double *TOPO, *TOPOp;

  double *maPARAMS, *maPARAMSp;

  double *innovations, *innovationVars, *wnV;

  double likelihood, sum1, sum2;

  int i, j, k;

  int numMissing;

  int iNode;
  
  double *TRAJ = prob->TRAJ;
  int trajLength = prob->trajLength;
  int arOrderMax = prob->arOrderMax;
  int maOrderMax = prob->maOrderMax;
  int nLags = arOrderMax + 1;
  int *topoBIN = prob->topoBIN;
  int *maBIN = prob->maBIN;
  int nNodes = prob->nNodes;
  int nInputs = prob->nInputs;
  int numParams = prob->numParams;
  double wnVariance = prob->wnVariance;



  TOPOp = (double *) malloc(sizeof(double) * nNodes * nNodes * nLags);

  maPARAMSp = (double *) malloc(sizeof(double) * nNodes * maOrderMax);

  paramsFromVector(&TOPOp, nNodes, arOrderMax,
		   &maPARAMSp, maOrderMax,
		   topoBIN, maBIN,
		   paramV, numParams);


  /* Convert partial parameters to full parameters */

  TOPO = (double *) malloc(sizeof(double) * nNodes * nNodes * nLags);

  for (i = 0; i < nNodes; i++){
    for (j = 0; j < nNodes; j++){

	/* Along the diagonal are AR parameters */

	if (i == j) {

	  /* No AR parameters at lag 0 */

	  *(TOPO + i * nLags + j * nNodes * nLags) = 0;

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



  maPARAMS = (double *) malloc(sizeof(double) * nNodes * maOrderMax);

  for (i = 0; i < nNodes; i++)
    levinsonDurbinExpoMA(maPARAMSp + i * maOrderMax,
			 maOrderMax,
			 maPARAMS + i * maOrderMax);



  /* Allocate memory for innovations, their variances and the white noise vector */

  innovations = (double *) malloc(sizeof(double) * trajLength);

  innovationVars = (double *) malloc(sizeof(double) * trajLength);

  wnV = (double *) malloc(sizeof(double) * trajLength);

  numMissing = 0;

  likelihood = 0;

  for (iNode = nInputs; iNode < nNodes; iNode++){

    sum1 = sum2 = 0;

    carmaCalcKalmanInnov(TRAJ, trajLength,
			 TOPO, topoBIN,
			 arOrderMax, nNodes,
			 maPARAMS, maOrderMax,
			 wnVariance, iNode,
			 &innovations, &innovationVars, &wnV,
			 &sum1, &sum2, &numMissing);


    likelihood += sum1 + (trajLength - numMissing) * log(sum2);

    prob->numMissing = MAX(prob->numMissing, numMissing);
 
  }


  /* Free allocated memory */
  
  free(TOPO);
  free(TOPOp);
  free(maPARAMS);
  free(maPARAMSp);
  free(innovations);
  free(innovationVars);
  free(wnV);  
  
  /*printf("\nlikelihood = %g",likelihood);*/

  return likelihood;

}
