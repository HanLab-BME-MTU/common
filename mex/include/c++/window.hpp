#ifndef WINDOW_HPP
# define WINDOW_HPP

# include <algorithm>

# include <vector.hpp>

template <typename T, std::size_t N>
class window
{
public:
  window(const vector<2, T> points[])
  {
    std::copy(points, points + N, points_);
  }
  
  const vector<2, T>& point(std::size_t i) const
  {
    assert(i < N);

    return points_[i];
  }

  static std::size_t card() { return N; }
  
private:
  vector<2, T> points_[N];
};

template <typename T>
const window<T, 4> neighb_c4()
{
  vector<2, T> points[4];

  points[0][0] = 1; points[0][1] = 0;
  points[1][0] = -1; points[1][1] = 0;
  points[2][0] = 0; points[2][1] = 1;
  points[3][0] = 0; points[3][1] = -1;

  const window<T, 4> n(points);
  
  return n;
}

template <typename T>
const window<T, 8> neighb_c8()
{
  vector<2, T> points[8];

  points[0][0] = -1; points[0][1] = -1;
  points[1][0] = 0; points[1][1] = -1;
  points[2][0] = 1; points[2][1] = -1;
  points[3][0] = -1; points[3][1] = 0;
  points[4][0] = 1; points[4][1] = 0;
  points[5][0] = -1; points[5][1] = 1;
  points[6][0] = 0; points[6][1] = 1;
  points[7][0] = 1; points[7][1] = 1;

  const window<T, 8> n(points);
  
  return n;
}

#endif /* WINDOW_HPP */
