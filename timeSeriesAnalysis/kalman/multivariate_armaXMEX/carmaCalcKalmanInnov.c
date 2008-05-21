#include <math.h>
#include <stdlib.h>
#include "carmaFitModel.h"
#include "./Utilities/Matrix_Operations/matOps.h"

#if !defined(MAX)

#define	MAX(A, B)	((A) > (B) ? (A) : (B))

#endif

static double **NewMat(int row, int col)
{
  int i;
  double **mat;

  mat = (double **) malloc(sizeof(double *) * row);
  for (i = 0; i < row; i++) mat[i]  = (double *) malloc(sizeof(double) * col);
  return mat;
}

static void FreeMat(double **mat, int row)
{
  int i;

  for (i = 0; i < row; i++) free(mat[i]);
  free(mat);
}



int carmaCalcKalmanInnov(
			  double *TRAJ,
			  int trajLength,
			  double *TOPO,
			  int *topoBIN,
			  int arOrderMax,
			  int nNodes,
			  double *maPARAMS,
			  int maOrderMax,                 
			  double wnVariance,
			  int iNode,
			  double **innovations,
			  double **innovationVars,
			  double **wnVector,
			  double *sum1,
			  double *sum2)
{

  /* Development Note: If in the future you are using models with long maximum  AR lags
     on some nodes but not others, it will increase efficiency to pass/calculate the
     ar Order for the specific node you are calculating innovations for  rather than the max */

  /* initialize counters */
  int j, k, l, m, n,p,q,r,tonsOCountersforEveryone;
  
  /* Calculate maximum time lag of model for sizing matrices */
  int maxOrder = MAX(arOrderMax,maOrderMax+1);

  /*Initialize Arrays for Modified parameters */
  double arParamMod[maxOrder];
  double maParamMod[maxOrder];

  /* Initialize state vectors, Covariance matrices etc. */    
  double stateVecT1_T[maxOrder];
  double stateVecT_T[maxOrder];
  double stateCovMatT_T[maxOrder][maxOrder];
  double stateCovMatT1_T[maxOrder][maxOrder];
  double observationVec[maxOrder];
  double G[maxOrder];
  double tmp = 0;


  for (k = 0; k < maxOrder; k++){
        
    stateVecT_T[k] = 0;        
    stateVecT1_T[k] = 0;
           
    /* Extract AR parameters for the specified node from TOPO */
    if (k < arOrderMax){
	/* Plus one because no ar param at lag = 0 */
	arParamMod[k] = *(TOPO + k + 1
				+ iNode*(arOrderMax+1) 
				+ iNode*nNodes*(arOrderMax+1));
    }
    else{arParamMod[k] = 0;}
        
    if (k < maOrderMax){
	maParamMod[k] = *(maPARAMS + k + iNode*maOrderMax);
    }
    else{maParamMod[k] = 0;}
        
    /*Initialize vector G */
    if (k < 1){
	G[k] = 1;
    }
    else {
	tmp = maParamMod[k-1];
	for (j = 1; j <= k; j++){
	  tmp += arParamMod[j-1]*G[k-j];
	}
	G[k] = tmp;
    }
        
    /* build "observation vector" */
        
    if (k == 0){
	observationVec[k] = 1;
    }
    else{
	observationVec[k] = 0;
    }

  }
  



  /* Get the initial covariance matrix */  
  
  double **tmp1;
  tmp1 = NewMat(maxOrder,maxOrder);
  double  tmp2[maxOrder][maxOrder];
  
  if (covKalmanInit(arParamMod,maParamMod,G,arOrderMax,maOrderMax,maxOrder,tmp1));
  
  else
    errHandle("Problem calling covKalmanInit!",__LINE__,__FILE__);

  /* Copy elements of covariance matrix so they are contiguous in memory... */
  
  for (m = 0; m < maxOrder; m++){
    for (n = 0;n < maxOrder; n++){
	tmp2[m][n] = tmp1[m][n];
    }
  } 
  

  FreeMat(tmp1,maxOrder);

  /*multiply initial covariance matrix by white noise variance */    
  matrixMultiplyConst(&tmp2[0][0],maxOrder,maxOrder,wnVariance,&stateCovMatT_T[0][0]);
  
  double matGGprime[maxOrder][maxOrder];
  double matGGprimeWN[maxOrder][maxOrder];
  double *GGprimeWN = &matGGprimeWN[0][0];    
  double *GGprime = &matGGprime[0][0];        
  matrixMultiply(G,maxOrder,1,G,1,maxOrder,GGprime);
  matrixMultiplyConst(GGprime,maxOrder,maxOrder,wnVariance, GGprimeWN);
   

  
  /* construct transition matrix F */    
    
  double F[maxOrder][maxOrder];
  double Fprime[maxOrder][maxOrder];
    
  for (k = 0; k < maxOrder; k++){
    for (l = 0; l < maxOrder; l++){
	if (l == k + 1){
	  F[k][l] = 1.0;
	  Fprime[l][k] = 1.0;
	}
	else if (k == maxOrder - 1){
	  F[k][l] = arParamMod[maxOrder-(l+1)];
	  Fprime[l][k] = arParamMod[maxOrder-(l+1)];
	}
	else {
	  F[k][l] = 0;
	  Fprime[l][k] = 0;
	}

    }
  }


  
  /* Determine number and location of connections */
  
  int nConnTo = 0;

  /* Over-allocate memory for connection indices */
  int *connFrom;
  connFrom = (int *) malloc((sizeof(double) * nNodes * nNodes));
  
  /* Go through binary adjacency matrices for each lag
     to determine which connections go to current node */

  for (m = 0; m < nNodes; m++){

    int   nParamsCurrConn = 0;

    for (r=0; r< arOrderMax+1; r++){
	
	if ( *(topoBIN + r + iNode*(arOrderMax+1) + m*nNodes*(arOrderMax+1))  & (m != iNode) ){
	  /*Increment counter if there's a parameter at this lag */
	  nParamsCurrConn++;
	}

    } /* end for r */

    /* If there was one or more parameters from this node,
	 store it's index */
    if (nParamsCurrConn){
	  connFrom[nConnTo] = m;
	  nConnTo++;
    }

  }/* end for m */

  /* Reallocate to proper length */
  connFrom = (int *) realloc(connFrom, sizeof(double)*nConnTo);


  /* Construct input coef matrix B */
  
  double B[nConnTo][maxOrder][arOrderMax+1];
    
  for (m = 0; m < nConnTo; m++){
    for (k = 0; k < maxOrder; k++){
	for (l = 0; l < arOrderMax+1 ; l++){
	  if (k == maxOrder-1){
	    B[m][k][l] = *(TOPO + arOrderMax - l + iNode*(arOrderMax+1) + connFrom[m]*nNodes*(arOrderMax+1));
	  }
	  else{
	    B[m][k][l] = 0;
	  }
	}
    }
  }

  /* Initialize variables for innovation calculation */

  double delta[maxOrder];
  double *inputVectors[nConnTo];
  int t2;
  double tmpMat[maxOrder][maxOrder];
  double tmpMat2[maxOrder][maxOrder];
  double tmpVec[maxOrder];
  int nanInInput = 0;
  int nanPresent = 0;
    
    
  /* Predict and filter timeseries */
  /* This is the meat! */                    
    
  /*Set the first few innovations/variances to NaN as they will be skipped below 
    depending on the xOrder */
  for (j = 0; j <= arOrderMax; j++){
    *(*innovations+j) = 0.0 / 0.0;
    *(*innovationVars+j) = 0.0/0.0;
  }            
    
    
  for (j = MAX(arOrderMax,0); j < trajLength; j++){

    /* get input timeseries */
        
        
    t2 = j - 1 + maxOrder - arOrderMax;

    
    if (arOrderMax >= 0){
	for (m = 0 ; m < nConnTo; m++){

            inputVectors[m] = TRAJ + (trajLength*connFrom[m]) + t2;

		/* check for NaNs in input vectors */        
		for (k = 0; k < arOrderMax+1; k++){                
		  if ( isnan( *(inputVectors[m] + k) ) ){
		    nanInInput = 1;
		  }                
		}

	}    




	        

    }

        
    if ( !nanInInput){
            
	/* If no NaN, predict the state at t+1 */

	/* apply transition matrix*/
	matrixMultiply(&F[0][0],maxOrder,maxOrder,&stateVecT_T[0],maxOrder,1,&stateVecT1_T[0]);

	if (arOrderMax >=0){
	  for (m = 0; m < nConnTo; m++){

	    /*multiply input matrix B by input trajectory */
	    matrixMultiply(&B[m][0][0],maxOrder,arOrderMax+1,inputVectors[m],arOrderMax+1,1,&tmpVec[0]);

	    /* add resulting vector to stateVec */
	    for (k = 0; k < maxOrder; k++){
		stateVecT1_T[k] += tmpVec[k];
	    }

	  }

	}

	/* Predict the covariance matrix of the state */

	/* calculate F * stateCovT_T */
	matrixMultiply(&F[0][0],maxOrder,maxOrder,&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0]);
	/* above times F' */
	matrixMultiply(&tmpMat[0][0],maxOrder,maxOrder,&Fprime[0][0],maxOrder,maxOrder,&tmpMat2[0][0]);
	/* add G*G' */
	matrixAdd(&tmpMat2[0][0],maxOrder,maxOrder,GGprimeWN,maxOrder,maxOrder,&stateCovMatT1_T[0][0]);
    }
    /*check if there is an observation at current timepoint */
        
        
    if ( isnan( *(TRAJ + iNode*trajLength + j) ) | nanInInput ){
	nanPresent = 1;
	nanInInput = 0;
    }
        
    if (nanPresent){
	for (k = 0; k < maxOrder; k++){

	  stateVecT_T[k] = stateVecT1_T[k];
                
	  for (m = 0; m < maxOrder; m++){
	    stateCovMatT_T[k][m] = stateCovMatT1_T[k][m];                                        
	  }                               
                
	}
            
	/* insure that covariance matrix is symmetric */
	matrixTranspose(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0]);
	matrixAdd(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0],maxOrder,maxOrder,&tmpMat2[0][0]);
	matrixMultiplyConst(&tmpMat2[0][0],maxOrder,maxOrder,.5,&stateCovMatT_T[0][0]);                                                            
            
	nanPresent = 0;
	/* Put NaN in place of innovations and variance */
	*(*innovations+j) = 0.0/0.0;
	*(*innovationVars+j) = 0.0/0.0;

    }        
    else{
	/* If all observations are present */

            
	/* Calculate innovations */
	*(*innovations+j) = *(TRAJ + trajLength*iNode + j) - stateVecT1_T[0];
	/* And variance in innovations */
	*(*innovationVars+j) = stateCovMatT1_T[0][0] + pow( *(TRAJ + nNodes*trajLength + iNode*trajLength + j ), 2); 

	*sum1 += log( *(*innovationVars+j) );
	*sum2 += pow( *(*innovations+j) , 2) / *(*innovationVars+j);

	/* calculate Delta */
	matrixMultiplyConst(&stateCovMatT1_T[0][0],maxOrder,1,( 1 / (*(*innovationVars+j)) ),&delta[0]);
	/* delta times the next innovation */
	matrixMultiplyConst(&delta[0],maxOrder,1,*(*innovations+j),&tmpVec[0]);        
	/* add this to previous state prediction */
	matrixAdd(&stateVecT1_T[0],maxOrder,1,&tmpVec[0],maxOrder,1,&stateVecT_T[0]);

	/* update covatiance matrix */
	matrixMultiply(&delta[0],maxOrder,1,&observationVec[0],1,maxOrder,&tmpMat[0][0]);
	matrixMultiply(&tmpMat[0][0],maxOrder,maxOrder,&stateCovMatT1_T[0][0],maxOrder,maxOrder,&tmpMat2[0][0]);
	matrixSubtract(&stateCovMatT1_T[0][0],maxOrder,maxOrder,&tmpMat2[0][0],maxOrder,maxOrder,&stateCovMatT_T[0][0]);

	/* insure that covariance matrix is symmetric */
	matrixTranspose(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0]);
	matrixAdd(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0],maxOrder,maxOrder,&tmpMat2[0][0]);
	matrixMultiplyConst(&tmpMat2[0][0],maxOrder,maxOrder,.5,&stateCovMatT_T[0][0]);	
  
    }


  }

  /* Free allocated memory */
  free(connFrom);

}
