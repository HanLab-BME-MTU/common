#ifndef BRESENHAM_HPP
# define BRESENHAM_HPP

# include <point_2d.hpp>

template <typename T, typename F>
void bresenham(const point_2d<T>& p1,
	       const point_2d<T>& p2,
	       F& func)
{
  int l, q, dqr, dqru;
      
  point_2d<T> p, incr1, incr2;
      
  int dy = p2.y - p1.y;
  int dx = p2.x - p1.x;
      
  incr1.y = dy < 0 ? (dy = -dy, -1) : 1;
  incr1.x = dx < 0 ? (dx = -dx, -1) : 1;
      
  if (dy >= dx)
    {
      dqr = dx << 1;
      dqru = dqr - (dy << 1);
      q = dqr - dy;
      l = dy;
      incr2.y = incr1.y;
    }
  else
    {
      dqr = dy << 1;
      dqru = dqr - (dx << 1);
      q = dqr - dx;
      l = dx;
      incr2.x = incr1.x;
    }
      
  p = p1;
      
  for (int d = l; d >= 0; --d)
    {
      func(p);
	  
      if (q > 0)
	{
	  p += incr1;
	  q += dqru;
	}
      else
	{
	  p += incr2;
	  q += dqr;
	}
    }
}

#endif /* BRESENHAM_HPP */
