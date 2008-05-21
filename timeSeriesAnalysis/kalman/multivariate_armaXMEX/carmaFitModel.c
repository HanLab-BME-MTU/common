#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "./Utilities/Carma_Utilities/carmaUtils.h"
#include "./Utilities/Matrix_Operations/matOps.h"
#include "prob.h"
#include "prob.c"
#include "./Minimizer/amoeba.h"
#include "./Minimizer/nrutil.h"
#include "carmaFitModel.h"



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{  
  
  /* Initialize timeseries & model parameter variables */
    
    double *TOPO;    
    double *maPARAMS;
    
    /* Initialize structure for model orders etc. */
    struct probCDT prob;

    probADT pProb = &prob;
    
    /* Retrieve model parameters and time series from
	 model structure. 
	 Fields are retrieved from structure with mxGetField.
	 Pointers for fields are retrieved with mxGetPr. */


    if (mxIsStruct(prhs[0])){    /* Check if first input is a structure */

	/* The mxGetField function throws a segmentation fault if provided a non-existent
	   field name, so you must verify that the field exists with mxGetFieldNumber */

	  if (mxGetFieldNumber(prhs[0],"arOrderMax") >= 0){
	    prob.arOrderMax = (int) *mxGetPr(mxGetField(prhs[0],0,"arOrderMax"));
	  } else {
	    errHandle("Cannot retrieve field 'arOrderMax' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"maOrderMax") >= 0){
	    prob.maOrderMax = (int) *mxGetPr(mxGetField(prhs[0],0,"maOrderMax"));
	  } else {
	    errHandle("Cannot retrieve field 'maOrderMax' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"maBIN") >= 0){
	    prob.maBIN = (int *) mxGetPr(mxGetField(prhs[0],0,"maBIN"));
	  } else {
	    errHandle("Cannot retrieve field 'maBIN' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"nNodes") >= 0){
	    prob.nNodes = (int) *mxGetPr(mxGetField(prhs[0],0,"nNodes"));
	  } else {
	    errHandle("Cannot retrieve field 'nNodes' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"TRAJ") >= 0){
	    prob.TRAJ = mxGetPr(mxGetField(prhs[0],0,"TRAJ"));
	  } else {
	    errHandle("Cannot retrieve field 'TRAJ' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"trajLength") >= 0){
	    prob.trajLength = (int) *mxGetPr(mxGetField(prhs[0],0,"trajLength"));
	  } else {
	    errHandle("Cannot retrieve field 'trajLength' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"TOPO") >= 0){
	    TOPO = mxGetPr(mxGetField(prhs[0],0,"TOPO"));
	  } else {
	    errHandle("Cannot retrieve field 'TOPO' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"maPARAMS") >= 0){
	    maPARAMS = mxGetPr(mxGetField(prhs[0],0,"maPARAMS"));
	  } else {
	    errHandle("Cannot retrieve field 'maPARAMS' from input model structure !",__LINE__,__FILE__);
	  }

	  if (mxGetFieldNumber(prhs[0],"topoBIN") >= 0){
	    prob.topoBIN = (int *) mxGetPr(mxGetField(prhs[0],0,"topoBIN"));
	  } else {
	    errHandle("Cannot retrieve field 'topoBIN' from input model structure !",__LINE__,__FILE__);
	  }
	  
    } else { /*if input is not a structure */
	errHandle("First input must be a structure containing the model!",__LINE__,__FILE__);
    }


    /* Convert ARMA parameters to partial parameters for call to
	 minimizer */

    
    int j;
    
    int k;
    int l;
    
    double *TOPOp;
    
    TOPOp = (double *) malloc(sizeof(double) * prob.nNodes * prob.nNodes * (prob.arOrderMax+1));

    for (j = 0; j < prob.nNodes; j++){
	for (k = 0; k < prob.nNodes; k++){
   
	  /* Along the diagonal the parameters are AR */
	
    
	if ( j == k) {
    
	    /* No AR parameters at lag = 0 */
    
	    *(TOPOp + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1)) = 0;

	    inverseLevinsonDurbinExpoAR( (TOPO + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1) + 1),
		                           prob.arOrderMax,
						   (TOPOp + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1)) + 1);
    
	    
      } else { 
	    /* Copy directly the X-parameters */

    
	    for (l = 0; l < (prob.arOrderMax+1); l++){
		*(TOPOp + l + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1)) = 
		  *(TOPO + l + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1));
	    }
	  }

	}
    }    
    
    
    double *maPARAMSp;
    maPARAMSp = (double *) malloc(sizeof(double) * prob.nNodes * prob.maOrderMax);
    for (j = 0; j < prob.nNodes; j++){
	
	inverseLevinsonDurbinExpoMA( maPARAMS + j*prob.maOrderMax,
					     prob.maOrderMax,
					     maPARAMSp + j*prob.maOrderMax);		
	
    }

    /* Assemble the parameters into a vector for passing to minimizer */
    
    double *paramV;

    vectorFromParams(TOPOp,maPARAMSp,pProb,
			   &paramV);
	 
	 
    prob.wnVariance = 1;
    double likelihood = 0;
    

    /*	likelihood = carmaNegLnLikelihood(paramV,pProb);*/

    /* Allocate memory for simplex */

    double **p;
    p = matrix(1,prob.numParams+1,1,prob.numParams);

    /* Create the initial simplex vertices based on initial guess */

    /* allocate vector for function values at simplex vertices */
    double *vertexLikelihoods = vector(1,prob.numParams+1);
    double simpRadius = .15;
    createSimplex(paramV,prob.numParams,simpRadius,p,vertexLikelihoods);

    /* Now evaluate the likelihood at each vertex */
    for (j = 1; j <= prob.numParams+1; j++){
	vertexLikelihoods[j] = carmaNegLnLikelihood(&p[j][1],pProb);	
    }
    int numEvals = 0;

    double fTol = .00000001;


    /* call the minimizer */
    if ( amoeba(p,vertexLikelihoods,prob.numParams,fTol,carmaNegLnLikelihood,&numEvals,pProb) ){
	printf("\ncarmaFitModel:Minimization successful!\n");
	plhs[3] = mxCreateDoubleScalar(1);
    } else{
	printf("\ncarmaFitModel:ERROR! Minimization failed! Maximum number of iterations reached!\n");
	plhs[3] = mxCreateDoubleScalar(0);
    }
    

    /* Dimensions of the TOPO matrix */
    mwSize topoDims[] = {prob.arOrderMax + 1,prob.nNodes, prob.nNodes};

    /* Allocate output matrices */
    /*TOPOLOGY matrix */
    plhs[0] = mxCreateNumericArray(3,topoDims,mxDOUBLE_CLASS,mxREAL);
    /* moving average parameters*/
    plhs[1] = mxCreateDoubleMatrix(prob.maOrderMax,prob.nNodes,mxREAL);
    /* and BIC */
    plhs[2] = mxCreateDoubleScalar(  vertexLikelihoods[1] + prob.numParams * log(prob.trajLength) );    

    double *TOPOfit = mxGetPr(plhs[0]);
    double *maPARAMSfit = mxGetPr(plhs[1]);
    double *TOPOfitp;
    TOPOfitp = (double *) malloc(sizeof(double) * prob.nNodes * prob.nNodes * (prob.arOrderMax+1));
    double *maPARAMSfitp;
    maPARAMSfitp = (double *) malloc(sizeof(double) * prob.nNodes * prob.maOrderMax);

    /* Extract partial parameters from simplex */
    paramsFromVector(&TOPOfitp,prob.nNodes,prob.arOrderMax,
			   &maPARAMSfitp,prob.maOrderMax,
			   prob.topoBIN,prob.maBIN,
			   &p[1][1],prob.numParams);

    /* convert these to 'regular' parameters for return */
    for (j = 0; j < prob.nNodes; j++){
	for (k = 0; k < prob.nNodes; k++){
	  
	  /* Along the diagonal the parameters are AR */
	  if ( j == k) {
	    
	    /* No AR parameters at lag = 0 */
	    *(TOPOfit + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1)) = 0;
	    
	    levinsonDurbinExpoAR( (TOPOfitp + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1) + 1),
					  prob.arOrderMax,
					  (TOPOfit + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1)) + 1);
	    
	    
	  } else {
	    /* Copy directly the X-parameters */
	    for (l = 0; l < (prob.arOrderMax+1); l++){
		*(TOPOfit + l + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1)) = 
		  *(TOPOfitp + l + j*(prob.arOrderMax+1) + k*prob.nNodes*(prob.arOrderMax+1));
	    }
	  }
	  
	}
    } 
    
    for (j = 0; j < prob.nNodes; j++){
	
	levinsonDurbinExpoMA( maPARAMSfitp + j*prob.maOrderMax,
				    prob.maOrderMax,
				    maPARAMSfit + j*prob.maOrderMax);	
	
	
    }

    




    /* Free allocated memory */
    free(paramV);
    free_vector(vertexLikelihoods,1,prob.numParams+1);
    free_matrix(p,1,prob.numParams+1,1,prob.numParams);

}

