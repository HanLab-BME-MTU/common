// MATLAB C-MEX 
//
// * createSparseDistanceMatrix.c *
//
// First version: Aaron Ponti - February 17th, 2003
//
// See createSparseDistanceMatrix.m for detailed help.
//

#include "mex.h"
#include <math.h>

/* If using a compiler that equates NaN to be zero, compile with
   
	 >> mex -DNAN_EQUALS_ZEROS createSparseDistanceMatrix.c

*/

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif

// Define constant
double percent_sparse=0.01;
double percent_step=0.01;

void calcDistMatrix1D(mxArray *plhs[], double *M, double *N, int Mrows, int Nrows, int nzmax, 
					  double threshold, double eps)
{
    double dist;
	int  i,j,k;

	// Sparse matrix pointers
	double *sr;
	int *irs, *jcs;

	// Get pointers for the output (sparse) matrix
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);

    // Fill the sparse distance matrix
	k=0;
	for (i=0;i<Nrows;i++) 
	{
        // Column indices
		jcs[i]=k;

        for (j=0;j<Mrows;j++)
		{        
            // Calculated distance dist
            dist=*(N+i)-*(M+j);            

			// If the distance is smaller than the threshold, store it 
			if (abs((int)dist)<=threshold)
			{
				if (dist==0) // Zero distances cannot be saved into a sparse matrix
					dist=eps;

				// Allocate new space if needed
				if(k >= nzmax) 
				{
					// Record old nzmaz
					int oldnzmax = nzmax;

					// Increase allocated space
					percent_sparse += percent_step;
					nzmax=(int)ceil((double)Mrows*(double)Nrows*percent_sparse);

					// Make sure nzmax increases at least by 1
					if(nzmax <= oldnzmax)
						nzmax=oldnzmax+1;

					// Reset storage space for non-zero elements of the sparse matrix
					mxSetNzmax(plhs[0], nzmax);

					// Reallocate memory
					mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
					mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(int)));

					// Get pointers to the new locations
					sr = mxGetPr(plhs[0]);
					irs = mxGetIr(plhs[0]);

				}
					
				// Now it is safe to save the calculated distance
				sr[k]=dist;
				irs[k]=j;
				k++;
			}

        }	
	}
	jcs[Nrows]=k;

	/*
	 *    Reallocate a fresh copy of the output matrix
	 *      Either 0 or k elements have been stored
	 */
	

	if (k>0)
		nzmax=k;
	else
		nzmax=1;

	// Reset storage space for non-zero elements of the sparse matrix
	mxSetNzmax(plhs[0], nzmax);

	// Reallocate memory
	mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(int)));

	// Get pointers to the new locations
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);
}

void calcDistMatrix2D(mxArray *plhs[], double *M, double *N, int Mrows, int Nrows, int nzmax, 
					  double threshold, double eps)
{
    double	mX=0, mY=0;
    double	nX=0, nY=0;
    double  *posM, *posN;
    double dist;
	int  i,j,k;

	// Sparse matrix pointers
	double *sr;
	int *irs, *jcs;
	
	// Get pointers for the output (sparse) matrix
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);

    // Fill the sparse distance matrix
	k=0;
	for (i=0;i<Nrows;i++) 
	{
        // Get source position
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
        
		// Column indices
		jcs[i]=k;

		for (j=0;j<Mrows;j++) 
		{        
            // Get target position
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);
            
			 // Calculate distance dist
            dist=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX));           

			// If the distance is smaller than the threshold, store it 
			if (dist<=threshold)
			{
				if (dist==0) // Zero distances cannot be saved into a sparse matrix
					dist=eps;

				// Allocate new space if needed
				if(k >= nzmax) 
				{	

					// Record old nzmaz
					int oldnzmax = nzmax;

					// Increase allocated space
					percent_sparse += percent_step;
					nzmax=(int)ceil((double)Mrows*(double)Nrows*percent_sparse);

					// Make sure nzmax increases at least by 1
					if(nzmax <= oldnzmax)
						nzmax=oldnzmax+1;

					// Reset storage space for non-zero elements of the sparse matrix
					mxSetNzmax(plhs[0], nzmax);

					// Reallocate memory
					mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
					mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(int)));

					// Get pointers to the new locations
					sr = mxGetPr(plhs[0]);
					irs = mxGetIr(plhs[0]);

				}
				
				// Now it is safe to save the calculated distance
				sr[k]=dist;
				irs[k]=j;
				k++;
			}

        }	
	}
	jcs[Nrows]=k;

	/*
	 *    Reallocate a fresh copy of the output matrix
	 *      Either 0 or k elements have been stored
	 */
	

	if (k>0)
		nzmax=k;
	else
		nzmax=1;

	// Reset storage space for non-zero elements of the sparse matrix
	mxSetNzmax(plhs[0], nzmax);

	// Reallocate memory
	mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(int)));

	// Get pointers to the new locations
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);
}


void calcDistMatrix3D(mxArray *plhs[], double *M, double *N, int Mrows, int Nrows, int nzmax, 
					  double threshold, double eps)
{
    double	mX=0, mY=0, mZ=0;
    double	nX=0, nY=0, nZ=0;
    double  *posM, *posN;
    int  i,j,k;
    double dist;

	// Sparse matrix pointers
	double *sr;
	int *irs, *jcs;

	// Get pointers for the output (sparse) matrix
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);


	// Fill the sparse distance matrix
    k=0;
	for (i=0;i<Nrows;i++) 
	{    
        // Get source position
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
		nZ=*(posN+2*Nrows);

		// Column indices
		jcs[i]=k;
        
        for (j=0;j<Mrows;j++) 
		{       
            // Get target position
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);
			mZ=*(posM+2*Mrows);
            
            // Calculate distance dist
            dist=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX)+(nZ-mZ)*(nZ-mZ));

			// If the distance is smaller than the threshold, store it 
			if (dist<=threshold)
			{
				if (dist==0) // Zero distances cannot be saved into a sparse matrix
					dist=eps;

				// Allocate new space if needed
				if(k >= nzmax)
				{
					// Record old nzmaz
					int oldnzmax = nzmax;

					// Increase allocated space
					percent_sparse += percent_step;
					nzmax=(int)ceil((double)Mrows*(double)Nrows*percent_sparse);

					// Make sure nzmax increases at least by 1
					if(nzmax <= oldnzmax)
						nzmax=oldnzmax+1;

					// Reset storage space for non-zero elements of the sparse matrix
					mxSetNzmax(plhs[0], nzmax);

					// Reallocate memory
					mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
					mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(int)));

					// Get pointers to the new locations
					sr = mxGetPr(plhs[0]);
					irs = mxGetIr(plhs[0]);

				}
				
				// Now it is safe to save the calculated distance
				sr[k]=dist;
				irs[k]=j;
				k++;
			}

        }
	}
	jcs[Nrows]=k;

	/*
	 *    Reallocate a fresh copy of the output matrix
	 *      Either 0 or k elements have been stored
	 */
	

	if (k>0)
		nzmax=k;
	else
		nzmax=1;

	// Reset storage space for non-zero elements of the sparse matrix
	mxSetNzmax(plhs[0], nzmax);

	// Reallocate memory
	mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(int)));

	// Get pointers to the new locations
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);
}


void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    // Variables
	int nzmax; // Number of positions allocated for the sparse matrix
	int cmplx; // We are not considering complex numbers (cmplx is set to 0)
	double threshold;
	double eps;

	// Initialize pointers for 2 inputs and 1 output
    double *M, *N;
	
	// Initialize int variables to store matrix dimensions
	int Mrows,Mcols;
	int Nrows,Ncols;

    // Check that the number of input and output parameters is valid
	if((nrhs<3) || (nrhs>4))
		mexErrMsgTxt("Three input parameters required.");
	if(nlhs > 1)
		mexErrMsgTxt("One output parameter required.");	

	// Input data must be double
	if(!(mxIsDouble(prhs[0]))) 
		mexErrMsgTxt("Input argument must be of type double.");
	
	// Set value for eps if not passed as a parameter
	switch (nrhs)
	{
	case 3: 
		{
			eps=1e-10;	
			break;
		}
	default:
		{
			eps=mxGetScalar(prhs[3]);
			break;
		}
	
	}

	// Read input parameter dimensions
	Mrows=mxGetM(prhs[0]);
	Mcols=mxGetN(prhs[0]);
	Nrows=mxGetM(prhs[1]);
	Ncols=mxGetN(prhs[1]);
	
	// Read parameter threshold
	threshold=mxGetScalar(prhs[2]);

	// Check input parameter dimension
	if ((Mcols>3) || (Ncols>3))
		mexErrMsgTxt("Point coordinates in more than 3 dimensions are not supported.");
	
	if (Mcols!=Ncols)
		mexErrMsgTxt("The points in the coordinate matrices have different number of dimensions.");

    // Assign pointers to each input
	M=mxGetPr(prhs[0]);
	N=mxGetPr(prhs[1]);	
	
	// Create output sparse matrix and set pointers to it
    nzmax=(int)ceil((double)Mrows*(double)Nrows*percent_sparse);
	if (nzmax==0)
		nzmax=1;

	// We don't need to save complex data
	cmplx=0;

	// Create a sparse matrix and assign it to the output parameter
	plhs[0]=mxCreateSparse(Mrows,Nrows,nzmax,cmplx);

   	// Call the distance function
	if (Mcols==1) { calcDistMatrix1D(plhs,M,N,Mrows,Nrows,nzmax,threshold,eps); }
	if (Mcols==2) { calcDistMatrix2D(plhs,M,N,Mrows,Nrows,nzmax,threshold,eps); }
	if (Mcols==3) { calcDistMatrix3D(plhs,M,N,Mrows,Nrows,nzmax,threshold,eps); }

}
