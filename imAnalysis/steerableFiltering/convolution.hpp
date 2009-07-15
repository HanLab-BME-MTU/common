#ifndef   	CONVOLUTION_HPP
# define   	CONVOLUTION_HPP

# include <vector>

# include <matrix.h> // for mxAssert uses

# include <GaussianDerivative1D.hpp>
# include <NonSeparableFilter2D.hpp>
# include <Image.hpp>

//////////////////////////////
// 2D separable convolution //
//////////////////////////////

template <typename T, int N, int M>
inline void convolve(const Image<T> & src,
		     const GaussianDerivative1D<N> & win1,
		     const GaussianDerivative1D<M> & win2,
		     Image<double> & dst,
		     bool mirror = false)
{
  mxAssert(dst.width() == src.width() && dst.height() == src.height(), "");

  int hside1 = win1.size() >> 1;
  int hside2 = win2.size() >> 1;

  int new_margin = std::max(hside1, hside2);

  if (mirror)
    src.borderMirror(new_margin);
  else
    src.borderReplicate(new_margin);

  Image<double> tmp(src.width(), src.height(), src.margin());

  double sum;

  // for all rows
  for (int j = 0; j < src.height(); ++j)
    for (int i = 0; i < src.width(); ++i)
      {
	sum = 0.0;
	for (int k = -hside1; k <= hside1; ++k)
	  sum += win1[k + hside1] * src(i - k, j);
	tmp(i, j) = sum;
      }

  if (mirror)
    tmp.borderMirror(new_margin);
  else
    tmp.borderReplicate(new_margin);

  // for all columns
  for (int i = 0; i < src.width(); ++i)
    for (int j = 0; j < src.height(); ++j)
      {
	sum = 0.0;
	for (int k = -hside2; k <= hside2; ++k)
	  sum += win2[k + hside2] * tmp(i, j - k);
	dst(i, j) = sum;
      }
}

//////////////////////////////////
// 2D non separable convolution //
//////////////////////////////////

template <typename T>
inline void convolve(const Image<T> & src,
		     const NonSeparableFilter2D & win,
		     Image<double> & dst,
		     bool mirror = false)
{
  mxAssert(dst.width() == src.width() && dst.height() == src.height(), "");

  int hside1 = win.width() >> 1;
  int hside2 = win.height() >> 1;

  int new_margin = std::max(hside1, hside2);

  if (mirror)
    src.borderMirror(new_margin);
  else
    src.borderReplicate(new_margin);

  double sum;
  
  for (int i = 0; i < src.width(); ++i)
    for (int j = 0; j < src.height(); ++j)
      {
	sum = 0.0;
	for (int x = -hside1; x <= hside1; ++x)
	  for (int y = -hside2; y <= hside2; ++y)
	    sum += win(x + hside1, y + hside2) * src(i - x, j - y);
	dst(i, j) = sum;
      }  
}

#endif	    /* !CONVOLUTION_HPP */
