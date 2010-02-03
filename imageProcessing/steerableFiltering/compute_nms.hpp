#ifndef		NMS_HPP
# define	NMS_HPP

# include <cmath>
# include <cassert>

# include <image.hpp>

template <typename T>
void compute_nms(const image<2, T> & src,
		 const image<2, double> & theta,
		 image<2, double> & dst)
{
  double t, ct, st, g;

  assert(src.width() == theta.width() && src.width() == dst.width() &&
	 src.height() == theta.height() && src.height() == dst.height());

  src.border_replicate(src.margin());

  for (int x = 0; x < src.width(); ++x)
    for (int y = 0; y < src.height(); ++y)
      {
	t = theta(x, y) + M_PI_2;

	ct = cos(t);
	st = sin(t);
	
	g = src(x, y);
	
	dst(x, y) = (g <= src(x + ct, y + st) ||
		     g <= src(x - ct, y - st)) ? 0 : g;
      }
}

#endif		/* !NMS_HPP */
