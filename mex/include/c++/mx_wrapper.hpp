# include <mex.h>

# include <image.hpp>

// Matlab mxArray size are row-wise. This trait provides the following
// conversion:
// [nrows, ncols] -> [width, height]
// [nrows, ncols, nslices] -> [width, height, depth]

template <int n>
struct sizeWrapper
{
  static bool eq(const int* s1, const int* s2);
  static void convert(const int* src, int* dst);
};
 
template <>
void sizeWrapper<2>::convert(const int* src, int* dst)
{
  dst[0] = src[1];
  dst[1] = src[0];
}

template <>
bool sizeWrapper<2>::eq(const int* s1, const int* s2)
{
  return s1[0] == s2[1] && s1[1] == s2[0];
}

template <>
void sizeWrapper<3>::convert(const int* src, int* dst)
{
  dst[0] = src[1];
  dst[1] = src[0];
  dst[2] = src[2];
}

template <>
bool sizeWrapper<3>::eq(const int* s1, const int* s2)
{
  return s1[0] == s2[1] && s1[1] == s2[0] && s1[2] == s2[2];
}

// This function copy the content of an image object into a mxArray
// from Matlab.

template <int n, typename T>
void image2mxArray(const image<n, T> & src, mxArray*& dst)
{
  const int* src_size = src.size();
  int src_size2[n];
  sizeWrapper<n>::convert(src_size, src_size2);

  if (dst)
    {
      if (mxGetNumberOfDimensions(dst) != n ||
	  !sizeWrapper<n>::eq(src_size, mxGetDimensions(dst)))
	{
	  mxDestroyArray(dst);
	  dst = mxCreateNumericArray(n, src_size2, mxDOUBLE_CLASS, mxREAL);
	}
    }
  else
    dst = mxCreateNumericArray(n, src_size2, mxDOUBLE_CLASS, mxREAL);
  
  src.raw_data(mxGetPr(dst));
}
