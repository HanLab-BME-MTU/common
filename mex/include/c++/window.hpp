#ifndef WINDOW_HPP
# define WINDOW_HPP

# include <algorithm>

# include <point_2d.hpp>

template <typename T, std::size_t N>
class window
{
public:
  window(const point_2d<T> points[])
  {
    std::copy(points, points + N, points_);
  }
  
  const point_2d<T>& point(std::size_t i) const
  {
    assert(i < N);

    return points_[i];
  }

  static std::size_t card() { return N; }
  
private:
  point_2d<T> points_[N];
};

template <typename T>
const window<T, 4>& neighb_c4()
{
  static const point_2d<T> points[] = {
    point_2d<T>( 1,  0),
    point_2d<T>(-1,  0),
    point_2d<T>( 0,  1),
    point_2d<T>( 0, -1)};
  
  static const window<T, 4> n(points);
  
  return n;
}

template <typename T>
const window<T, 8>& neighb_c8()
{
  static const point_2d<T> points[] = {
    point_2d<T>(-1, -1), point_2d<T>( 0, -1),
    point_2d<T>( 1, -1), point_2d<T>(-1,  0),
    point_2d<T>( 1,  0), point_2d<T>(-1,  1),
    point_2d<T>( 0,  1), point_2d<T>( 1,  1)};
  
  static const window<T, 8> n(points);
  
  return n;
}

#endif /* WINDOW_HPP */
