#ifndef WINDOW_HPP
# define WINDOW_HPP

# include <algorithm>

# include <vector.hpp>

template <unsigned n, typename T>
class window
{
public:
  window() {}

  window(const vector<2,T> points[])
  {
    std::copy(points, points + n, points_);
  }
  
  const vector<2,T> & point(unsigned i) const { return points_[i]; }

  vector<2,T> & point(unsigned i) { return points_[i]; }

  static unsigned size() { return n; }
  
private:
  vector<2,T> points_[n];
};

template <typename T>
inline
window<4,T> neighb_c4()
{
  static window<4,T> win;
  static bool first = true;
  
  if (first)
    {
      win.point(0)[0] = 1;	win.point(0)[1] = 0;
      win.point(1)[0] = -1;	win.point(1)[1] = 0;
      win.point(2)[0] = 0;	win.point(2)[1] = 1;
      win.point(3)[0] = 0;	win.point(3)[1] = -1;
      first = false;
    }
  return win;
}

template <typename T>
inline
window<8,T> neighb_c8()
{
  static window<8,T> win;
  static bool first = true;

  if (first)
    {
      win.point(0)[0] = -1;	win.point(0)[1] = -1;
      win.point(1)[0] = 0;	win.point(0)[1] = -1;
      win.point(2)[0] = 1;	win.point(0)[1] = -1;
      win.point(3)[0] = -1;	win.point(0)[1] = 0;
      win.point(4)[0] = 1;	win.point(0)[1] = 0;
      win.point(5)[0] = -1;	win.point(0)[1] = 1;
      win.point(6)[0] = 0;	win.point(0)[1] = 1;
      win.point(7)[0] = 1;	win.point(0)[1] = 1;
      first = false;
    }
  
  return win;
}

#endif /* WINDOW_HPP */
