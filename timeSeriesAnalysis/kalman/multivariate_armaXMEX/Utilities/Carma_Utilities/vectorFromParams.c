#include <stdlib.h>
#include "carmaUtils.h"
#include "../../prob.h"
#include "../../prob.c"

/* Converts CARMA model parameters into a vector for passing to a minimizer, */
/*                using the connections allowed by adjacency matrix tryCONN */

/****** Synopsis **********/
/* under construction  ;) */


int vectorFromParams(double *TOPO, double *maPARAMS, void *probP,
			   double **paramVector)
{

  /* Initialize Counters */
  int m, n, r;

  /* Retrieve ADT from  void pointer */
  probADT prob = (probADT) probP;


  /* Determine maximum possible number of parameters in the model */
  int maxNumParams = (prob->nNodes * prob->nNodes * (prob->arOrderMax+1)) + (prob->nNodes * prob->maOrderMax);


  /* allocate paramVector array to this size, and
    return the  pointer */
  *paramVector = (double *) malloc((sizeof(double) * maxNumParams));  



  /* Concatenate parameters into vector and count them */
  /* Set the first parameter to zero for numerical recipes */
  **paramVector = 0;
  int paramCtr = 1;
  int strt;

  /* First add the AR and X params (C params) */ 

  for (m = 0; m < prob->nNodes; m++){
    for (n = 0; n < prob->nNodes; n++){	  
	
	/* if it is a node instead of a connection, skip lag 0 */
	if (m == n)
	  strt = 1;
	else
	  strt = 0;
	
	for (r = strt; r <= prob->arOrderMax; r++){
	  
	  /* Get the params for node/connection, if a parameter is present
	     as indicated by topoBIN */	    
	  if( *(prob->topoBIN + r + (n*(prob->arOrderMax+1)) + (m*prob->nNodes*(prob->arOrderMax+1)) ) ){
	    
	    *(*paramVector + paramCtr) = *(TOPO + r + (n*(prob->arOrderMax+1)) + (m*prob->nNodes*(prob->arOrderMax+1)));
	    paramCtr++;
	    
	  }
	  
	}
	
    }/*end for n */
    
  }/*end for m */
  

  
  /* Now add the moving average parameters */

  for (m = 0; m < prob->nNodes; m++){
    
    /* Retrieve the Moving Average (MA) parameters at this node, 
	 if maBIN indicates a parameter present */
    for (n = 0; n < prob->maOrderMax; n++){
	
	if ( *(prob->maBIN + n + (m*prob->maOrderMax)) ){

	  *(*paramVector+paramCtr) = *(maPARAMS + n + (m*prob->maOrderMax));
	  paramCtr++;
	 
	}

    } /* end for n */
  }/* end for m */
  
  
  /* Reallocate parameter vector for actual number of params */
  *paramVector = (double *) realloc(*paramVector, sizeof(double) * paramCtr);
  

  
  /* Return number of parameters, which is one less than paramCtr due to 
     zero at beginning for numerical recipes */
  prob->numParams = --paramCtr;
  
  
} /* Close vectorFromParams */
