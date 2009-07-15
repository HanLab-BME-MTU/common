#ifndef		NMS_HPP
# define	NMS_HPP

# include <cmath>

# include <matrix.h> // for mxAssert uses

# include <Image.hpp>

template <typename T>
void nms(const Image<T> & src, const Image<double> & angle, Image<double> & dst)
{
  double theta, ct, st, g;

  mxAssert(src.width() == angle.width() && src.width() == dst.width() &&
	   src.height() == angle.height() && src.height() == dst.height(), "");

  src.borderReplicate(src.margin());

  for (int x = 0; x < src.width(); ++x)
    for (int y = 0; y < src.height(); ++y)
      {
	theta = angle(x, y) + M_PI_2;

	ct = cos(theta);
	st = sin(theta);
	
	g = src(x, y);
	
	dst(x, y) = (g <= src(x + ct, y + st) ||
		     g <= src(x - ct, y - st)) ? 0 : g;
      }
}

#endif		/* !NMS_HPP */
