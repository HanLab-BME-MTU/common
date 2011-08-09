/*
 *
 * Compilation:
 * Mac/Linux: mex -I/usr/local/include -I../../mex/include/c++ /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a geodesicDistance.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\..\extern\mex\include\gsl-1.14" -I"..\..\mex\include\c++" "..\..\..\extern\mex\lib\gsl.lib" "..\..\..\extern\mex\lib\cblas.lib" -output geodesicDistance geodesicDistance.cpp
 */


#include <mex.h>

#include <image.hpp>
#include <mx_wrapper.hpp>
#include <wwindow.hpp>

template <unsigned N>
static void dispatch(int nlhs, mxArray *plhs[],
		     int nrhs, const mxArray *prhs[])
{
  //////////////////////////
  // Get input parameters //
  //////////////////////////

  int size[N];
  sizeWrapper<N>::convert(mxGetDimensions(prhs[0]), size);
  unsigned char* ptr = (unsigned char*) mxGetData(prhs[0]);
  image<N, unsigned char> X(size);
  X.fill(ptr);

  ptr = (unsigned char*) mxGetData(prhs[1]);
  image<N, unsigned char> Y(size);
  Y.fill(ptr);

  /////////////////////////
  // Initialize distance //
  /////////////////////////

  image<N, double> D(size);
  D.fill(std::numeric_limits<double>::max());
  D.border_assign(D.margin(), std::numeric_limits<double>::max());

  ///////////////////////////////
  // Compute geodesic distance //
  ///////////////////////////////

  geodesic_distance(X,Y,D);

  //////////
  // Save //
  //////////

  if (nlhs > 0)
    image2mxArray(D, plhs[0]);
}

template <typename W>
static void geodesic_distance(const image<2,unsigned char> & X,
							  const image<2, unsigned char> & Y,
							  image<2,W> & D)
{
	// Initialize D
	for (int i = 0; i < D.width(); ++i)
		for (int j = 0; j < D.height(); ++j)
			if (X(i,j))
				D(i,j) = 0;
	
	const wwindow<2,2,int,W> & win_fwd = chanfrein11_fwd<int,W>();
	const wwindow<2,2,int,W> & win_bkd = chanfrein11_bkd<int,W>();
	
	vector<2,int> p, q;

  // Forward sweep
  for (int i = 0; i < D.width(); ++i)
    for (int j = 0; j < D.height(); ++j)
      {
	p[0] = i; p[1] = j;

	W min = D[p];
		  
	if (min == 0)
	  continue;
	
	for (int r = 0; r < win_fwd.size(); ++r)
	  {
	    q = p + win_fwd.point(r);

	    W cur = D[q] + win_fwd.weight(r);

	    if (cur < min && Y[q])
	      min = cur;
	  }
		  
	D[p] = min;
      }

  // Backward sweep
  for (int i = D.width() - 1; i >= 0; --i)
    for (int j = D.height() - 1; j >= 0; --j)
      {
	p[0] = i; p[1] = j;

	W min = D[p];
	      
	if (min == 0)
	  continue;

	for (int r = 0; r < win_bkd.size(); ++r)
	  {
	    q = p + win_bkd.point(r);

	    W cur = D[q] + win_bkd.weight(r);

	    if (cur < min && Y[q])
	      min = cur;
	  }

	D[p] = min;
      }
}

template <typename W>
static void geodesic_distance(const image<3, unsigned char> & X,
							  const image<3, unsigned char> & Y,
							  image<3, W> & D)
{
	// Initialize D
	for (int i = 0; i < D.width(); ++i)
		for (int j = 0; j < D.height(); ++j)
			for (int k = 0; k < D.depth(); ++k)
				if (X(i,j,k))
					D(i,j,k) = 0;
	
	const wwindow<3,3,int,W> & win_fwd = chanfrein111_fwd<int,W>();
	const wwindow<3,3,int,W> & win_bkd = chanfrein111_bkd<int,W>();	
	
  vector<3,int> p, q;

  // Forward sweep
  for (int i = 0; i < D.width(); ++i)
    for (int j = 0; j < D.height(); ++j)
      for (int k = 0; k < D.depth(); ++k)
      {
	p[0] = i; p[1] = j; p[2] = k;

	W min = D[p];
	      
	if (min == 0)
	  continue;

	for (int r = 0; r < win_fwd.size(); ++r)
	  {
	    q = p + win_fwd.point(r);

	    W cur = D[q] + win_fwd.weight(r);

	    if (cur < min && Y[q])
	      min = cur;
	  }

	D[p] = min;
      }

  // Backward sweep
  for (int i = D.width() - 1; i >= 0; --i)
    for (int j = D.height() - 1; j >= 0; --j)
      for (int k = D.depth() - 1; k >= 0; --k)
      {
	p[0] = i; p[1] = j; p[2] = k;

	W min = D[p];
	      
	if (min == 0)
	  continue;

	for (int r = 0; r < win_bkd.size(); ++r)
	  {
	    q = p + win_bkd.point(r);

	    W cur = D[q] + win_bkd.weight(r);

	    if (cur < min && Y[q])
	      min = cur;
	  }

	D[p] = min;
      }
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  /////////////////////////////////////////////////
  // Check number of input and output parameters //
  /////////////////////////////////////////////////

  if (nrhs < 2)
    mexErrMsgTxt("Two input arguments required.");

  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments");

  ////////////////////////////
  // Check input parameters //
  ////////////////////////////

  int ndim1 = mxGetNumberOfDimensions(prhs[0]);
  int ndim2 = mxGetNumberOfDimensions(prhs[1]);
  
  if (ndim1 != ndim2)
    mexErrMsgTxt("X and Y must have the same number of dimensions.");

  if (ndim1 != 2 && ndim1 != 3)
    mexErrMsgTxt("Invalid number of dimensions (must be 2 or 3).");

  if (!mxIsClass(prhs[0],"logical"))
    mexErrMsgTxt("X is not a logical matrix.");

  if (!mxIsClass(prhs[1],"logical"))
    mexErrMsgTxt("Y is not a logical matrix.");

  //////////////
  // Dispatch //
  //////////////

  switch (ndim1)
    {
    case 2: dispatch<2>(nlhs, plhs, nrhs, prhs); break;
    case 3: dispatch<3>(nlhs, plhs, nrhs, prhs); break;
    }
}
