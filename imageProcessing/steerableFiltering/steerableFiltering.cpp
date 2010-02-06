#include <mex.h>

#include <image.hpp>
#include <mx_wrapper.hpp>
#include <compute_nms.hpp>
#include <steerable_filtering.hpp>

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
    mexErrMsgTxt("Invalid dimension for sigma argument.");

  //////////////////////////
  // Get input parameters //
  //////////////////////////

  int size[2];
  sizeWrapper<2>::convert(mxGetDimensions(prhs[0]), size);
  double* ptr = mxGetPr(prhs[0]);
  image<2, double> ima(size);
  ima.fill(ptr);

  int str_length = mxGetNumberOfElements(prhs[1]);

  char* buf = new char[str_length];

  mxGetString(prhs[1], buf, str_length + 1);

  std::string method_name(buf, str_length);

  delete[] buf;

  const std::string valid_method_names[] = {
    "2ndGaussian",
    "UnserM1",
    "UnserM2",
    "UnserM4"
  };

  int method_id = 0;

  while (method_id < 4 && method_name != valid_method_names[method_id])
    ++method_id;

  if (method_id == 4)
    mexErrMsgTxt((method_name + " is not a valid method").c_str());

  ptr = mxGetPr(prhs[2]);

  if (*ptr <= 0)
    mexErrMsgTxt("sigma must be strictly positive.");

  double sigma = *ptr;
  
  ///////////////////////
  // Compute filtering //
  ///////////////////////
  
  image<2, double> res(ima.size(), 2); // for nms safety
  image<2, double> theta(ima.size());

  switch (method_id)
    {
    case 0: snd_gaussian_filtering(ima, sigma, res, theta); break;
    case 1: unser_m1_filtering(ima, sigma, res, theta); break;
    case 2: unser_m2_filtering(ima, sigma, res, theta); break;
    case 3: unser_m4_filtering(ima, sigma, res, theta); break;
    }
  
  ////////////////////////////
  // Allocate output arrays //
  ////////////////////////////

  if (nlhs > 0) image2mxArray(res, plhs[0]);
  if (nlhs > 1) image2mxArray(theta, plhs[1]);

  if (nlhs > 2)
    {
      /////////////////////////////////////
      // Compute non-maximal suppression //
      /////////////////////////////////////

      image<2, double> nms(ima.size());
  
      compute_nms(res, theta, nms);

      image2mxArray(nms, plhs[2]);
    }
}
