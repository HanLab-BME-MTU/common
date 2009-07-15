#include <mex.h>

#include <Image.hpp>
#include <nms.hpp>
#include <steerableFiltering.hpp>

////////////////////////////////////////////////////////////////
//							      //
// Input:						      //
//   							      //
// I		input image.				      //
//							      //
// method	string designating which kind of steerable    //
//              filtering should be used. Choice are:         //
//              - '2ndGaussian' 2nd Derivative of a Gaussian  //
//              - 'UnserM1'     1st order Unser edge filter   //
//              - 'UnserM2'     2nd order Unser ridge filter  //
//							      //
// sigmaPSF	standard deviation of the filter.	      //
//   							      //
// Output:						      //
//   							      //
// R		filtering response.			      //
//							      //
// T		optimal angle.				      //
//							      //
// NMS		non-maximal suppression (optional).	      //
//							      //
////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  /////////////////////////////////////////////////
  // Check number of input and output parameters //
  /////////////////////////////////////////////////

  if (nrhs < 3)
    mexErrMsgTxt("Three input arguments required.");

  if (nlhs > 3)
    mexErrMsgTxt("Too many output arguments");

  ////////////////////////////
  // Check input parameters //
  ////////////////////////////

  if (mxGetNumberOfDimensions(prhs[0]) != 2)
    mexErrMsgTxt("Invalid dimension for I argument.");

  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("I is not a double-precision matrix.");

  if (!mxIsChar(prhs[1]))
    mexErrMsgTxt("method argument is not a string.");

  if (mxGetNumberOfElements(prhs[2]) != 1)
    mexErrMsgTxt("Invalid dimension for sigmaPSF argument.");

  //////////////////////////
  // Get input parameters //
  //////////////////////////

  const mwSize* size = mxGetDimensions(prhs[0]);
  double* ptr = mxGetPr(prhs[0]);
  Image<double> I(ptr, size[1], size[0]);

  int strLength = mxGetNumberOfElements(prhs[1]);

  char* methodNameC = new char[strLength];

  mxGetString(prhs[1], methodNameC, strLength + 1);

  std::string methodName(methodNameC, strLength);

  delete[] methodNameC;

  const std::string validMethodNames[] = {
    "2ndGaussian",
    "UnserM1",
    "UnserM2"
  };

  int methodID = 0;

  while (methodID < 3 && methodName != validMethodNames[methodID])
    ++methodID;

  if (methodID == 3)
    mexErrMsgTxt((methodName + " is not a valid method").c_str());

  ptr = mxGetPr(prhs[2]);

  if (*ptr <= 0)
    mexErrMsgTxt("sigmaPSF must be strictly positive.");

  double sigmaPSF = *ptr;
  
  ///////////////////////
  // Compute filtering //
  ///////////////////////
  
  Image<double> R(I.width(), I.height(), 2); // for nms safety
  Image<double> T(I.width(), I.height());

  switch (methodID)
    {
    case 0: sndGaussianFiltering(I, sigmaPSF, R, T); break;
    case 1: unserFilteringM1(I, sigmaPSF, R, T); break;
    case 2: unserFilteringM2(I, sigmaPSF, R, T); break;
    }

  ////////////////////////////
  // Allocate output arrays //
  ////////////////////////////

  if (nlhs > 0) R.convertToMxArray(plhs[0]);
  if (nlhs > 1) T.convertToMxArray(plhs[1]);

  if (nlhs > 2)
    {
      /////////////////////////////////////
      // Compute non-maximal suppression //
      /////////////////////////////////////

      Image<double> NMS(I.width(), I.height());
  
      nms(R, T, NMS);

      NMS.convertToMxArray(plhs[2]);
    }
}
