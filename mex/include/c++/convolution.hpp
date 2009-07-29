#ifndef   	CONVOLUTION_HPP
# define   	CONVOLUTION_HPP

# include <vector>
# include <cassert>

# include <gaussian_derivative_1d.hpp>
# include <non_separable_filter_2d.hpp>
# include <image.hpp>

//////////////////////////////
// 2D separable convolution //
//////////////////////////////

template <typename T, int N, int M>
inline void convolve(const image<T> & src,
		     const gaussian_derivative_1d<N> & win1,
		     const gaussian_derivative_1d<M> & win2,
		     image<double> & dst,
		     bool mirror = false)
{
  assert(dst.width() == src.width() && dst.height() == src.height());

  int hside1 = win1.size() >> 1;
  int hside2 = win2.size() >> 1;

  int new_margin = std::max(hside1, hside2);

  if (mirror)
    src.border_mirror(new_margin);
  else
    src.border_replicate(new_margin);

  image<double> tmp(src.width(), src.height(), src.margin());

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
    tmp.border_mirror(new_margin);
  else
    tmp.border_replicate(new_margin);

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
inline void convolve(const image<T> & src,
		     const non_separable_filter_2d & win,
		     image<double> & dst,
		     bool mirror = false)
{
  assert(dst.width() == src.width() && dst.height() == src.height());

  int hside1 = win.width() >> 1;
  int hside2 = win.height() >> 1;

  int new_margin = std::max(hside1, hside2);

  if (mirror)
    src.border_mirror(new_margin);
  else
    src.border_replicate(new_margin);

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
