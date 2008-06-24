#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "prob.h"
#include "prob.c"
#include "carmaFitModel.h"

#if !defined(MAX)

#define	MAX(A, B)	((A) > (B) ? (A) : (B))

#endif

double carmaNegLnLikelihood(double *paramVnr, void *probP)
			  

{

  probADT prob = (probADT) probP;

  /* Convert the 1-index numerical recipes pointer to 0 index */
  double *paramV = paramVnr + 1;

  double likelihood = 0;
  double sum1;
  double sum2;
  int iNode = 2;

  /* Extract process orders etc. from problem structure for readability */
  
  double *TRAJ = prob->TRAJ;
  int trajLength = prob->trajLength;
  int arOrderMax = prob->arOrderMax;
  int maOrderMax = prob->maOrderMax;
  int *topoBIN = prob->topoBIN;
  int *maBIN = prob->maBIN;
  int nNodes = prob->nNodes;
  int nInputs = prob->nInputs;
  int numParams = prob->numParams;
  double wnVariance = prob->wnVariance;

  /* Use these orders to extract partial parameters from vector */

  double *TOPOp;
  TOPOp = (double *) malloc(sizeof(double) * nNodes * nNodes * (arOrderMax+1));
  double *maPARAMSp;
  maPARAMSp = (double *) malloc(sizeof(double) * nNodes * maOrderMax);

  paramsFromVector(&TOPOp,nNodes,arOrderMax,
			 &maPARAMSp,maOrderMax,
			 topoBIN,maBIN,
			 paramV,numParams);


  /* Convert partial ARMA parameters to ARMA parameters */
  int j;
  int k;
  int l;
  double *TOPO;
  TOPO = (double *) malloc(sizeof(double) * nNodes * nNodes * (arOrderMax+1));

  for (j = 0; j < nNodes; j++){
    for (k = 0; k < nNodes; k++){

	/* Along the diagonal the parameters are AR */
	if ( j == k) {

	  /* No AR parameters at lag = 0 */
	  *(TOPO + j*(arOrderMax+1) + k*nNodes*(arOrderMax+1)) = 0;

	  levinsonDurbinExpoAR( (TOPOp + j*(arOrderMax+1) + k*nNodes*(arOrderMax+1) + 1),
					arOrderMax,
					(TOPO + j*(arOrderMax+1) + k*nNodes*(arOrderMax+1)) + 1);
	  
	  
	} else {
	  /* Copy directly the X-parameters */
	  for (l = 0; l < (arOrderMax+1); l++){
	    *(TOPO + l + j*(arOrderMax+1) + k*nNodes*(arOrderMax+1)) = 
		*(TOPOp + l + j*(arOrderMax+1) + k*nNodes*(arOrderMax+1));
	  }
	}
	
    }
  } 

  /* Convert partial MA params to MA params */
  double *maPARAMS;
  maPARAMS = (double *) malloc(sizeof(double) * nNodes * maOrderMax);

  for (j = 0; j < nNodes; j++){
	
    levinsonDurbinExpoMA( maPARAMSp + j*maOrderMax,
				  maOrderMax,
				  maPARAMS + j*maOrderMax);	

	
  }


  /* Allocate memory for innovations, thier variances and the white noise vector */
  double *innovations, *innovationVars, *wnV;
  innovations = (double *) malloc(sizeof(double) * trajLength);
  innovationVars = (double *) malloc(sizeof(double) * trajLength);
  wnV = (double *) malloc(sizeof(double) * trajLength);    
  int numMissing = 0;

  for (iNode = nInputs; iNode < nNodes; iNode++){

    sum1 = 0;
    sum2 = 0;

    carmaCalcKalmanInnov(TRAJ,
				 trajLength,
				 TOPO,topoBIN,
				 arOrderMax,nNodes,
				 maPARAMS,maOrderMax,
				 wnVariance,iNode,
				 &innovations,
				 &innovationVars,&wnV,
				 &sum1,&sum2,
				 &numMissing);
    
    
    
    likelihood += sum1 + (trajLength-numMissing) * log(sum2);

    /* Store the number of missing observations in problem structure */
    prob->numMissing = MAX(prob->numMissing,numMissing);
    
            
    
  }/* End for iNode */


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
