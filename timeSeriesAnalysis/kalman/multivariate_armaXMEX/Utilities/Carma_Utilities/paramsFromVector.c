#include <stdlib.h>
#include "carmaUtils.h"

/* Extracts CARMA model parameters from vector used by minimizer, */


/****** Synopsis **********/



paramsFromVector(double **TOPO, int nNodes, int arOrderMax, 
		     double **maPARAMS, int maOrderMax, 
		     int *topoBIN, int *maBIN, 
		     double *paramVector, int numParams)
{

/* Initialize Counters */
 int m, n, r;
 
 /* Retrieve Parameters from Vector */

 int strt;
 int paramCtr = numParams;
 paramCtr--;

 /* First, retrieve the MA parameters */

 for (m = nNodes-1; m >= 0; m--){
    
   /* Retrieve the Moving Average (MA) parameters at this node, 
	if maBIN indicates a parameter present */
   for (n = maOrderMax-1; n >= 0; n--){

     if ( *(maBIN + n + (m*maOrderMax)) ){

	 *(*maPARAMS + n + (m*maOrderMax)) = *(paramVector+paramCtr);
	 paramCtr--;
     } else {
	 *(*maPARAMS + n + (m*maOrderMax)) = 0.0;
     }

   }/* end for n */
 }/* end for m */

 /* Now retrieve the AR and X params (C params) */ 

 for (m = nNodes-1; m >= 0; m--){
    for (n = nNodes-1; n >=0; n--){


	/* Get the params for node/connection, excluding NaNs */
	  
	/* No AR params at lag 0 */
	if (m == n){
	  strt = 1;	
	  *(*TOPO + (n*(arOrderMax+1)) + (m*nNodes*(arOrderMax+1))) = 0;
	} else {
	  strt = 0;
	}

	for (r = arOrderMax; r >= strt; r--){
	  
	  /* If topoBIN indicates a parameter is present, get it */
	  if (*(topoBIN + r + (n*(arOrderMax+1)) + (m*nNodes*(arOrderMax+1)) ) ){

	    *(*TOPO + r + (n*(arOrderMax+1)) + (m*nNodes*(arOrderMax+1))) = *(paramVector+paramCtr);
	    paramCtr--;
	  } else{
	    *(*TOPO + r + (n*(arOrderMax+1)) + (m*nNodes*(arOrderMax+1))) = 0.0;
	  }
	}
    
    }/*end for n */
    
 }/*end for m */
		   

} /* Close paramsFromVector */
