#ifndef POINT_2D_HPP
# define POINT_2D_HPP

# include <iostream>
# include <cassert>
# include <cmath>
# include <limits>

template <typename T>
struct point_2d
{
  point_2d() : x(0), y(0) {}

  point_2d(T _x, T _y) : x(_x), y(_y) {}

  template <typename U>
  point_2d(const point_2d<U>& rhs)
  {
    x = rhs.x;
    y = rhs.y;		
  }

  template <typename U>
  point_2d & operator=(const point_2d<U>& rhs)
  {
    x = rhs.x;
    y = rhs.y;
    
    return *this;
  }

  // Array-like accessors. Do NOT change the field x, y.
  const T & operator[](int dim) const
  {
    assert(dim >= 0 && dim < 2);

    return reinterpret_cast<const T*>(this)[dim];
  }

  T & operator[](int dim)
  {
    assert(dim >= 0 && dim < 2);

    return reinterpret_cast<T*>(this)[dim];
  }

  bool operator==(const point_2d& rhs) const
  {
    return x == rhs.x && y == rhs.y;
  }

  bool operator!=(const point_2d& rhs) const
  {
    return !this->operator==(rhs);
  }

  point_2d operator+(const point_2d& rhs) const
  {
    return point_2d(x + rhs.x, y + rhs.y);
  }

  point_2d & operator+=(const point_2d& rhs)
  {
    x += rhs.x;
    y += rhs.y;

    return *this;
  }

  point_2d & operator+=(const T& rhs)
  {
    x += rhs;
    y += rhs;

    return *this;
  }

  point_2d operator-(const point_2d& rhs) const
  {
    return point_2d(x - rhs.x, y - rhs.y);
  }

  point_2d & operator-=(const point_2d& rhs)
  {
    x -= rhs.x;
    y -= rhs.y;

    return *this;
  }

  point_2d & operator-=(const T& rhs)
  {
    x -= rhs;
    y -= rhs;

    return *this;
  }

  template <typename U>
  point_2d<U> operator/(const U& v) const
  {
    return point_2d<U>(x / v, y / v);
  }

  double norm2() const { return x * x + y * y; }

  double norm() const { return sqrt(norm2()); }

  T x, y;
};

template <typename T>
std::ostream& operator<<(std::ostream& o, const point_2d<T>& p)
{
  o << '(' << p.x << ',' << p.y << ')';

  return o;
}

template <typename T1, typename T2>
double dist2(const point_2d<T1>& p1, const point_2d<T2>& p2)
{
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;

  return dx * dx + dy * dy;
}
    
template <typename T1, typename T2>
double dist(const point_2d<T1>& p1, const point_2d<T2>& p2)
{
  return sqrt(dist2(p1, p2));
}

template <typename T>
double dot_product(const point_2d<T>& p1, const point_2d<T>& p2)
{
  return p1.x * p2.x + p1.y * p2.y;
}

template <typename T>
bool areLinearlyIndependent(const point_2d<T>& p1, const point_2d<T>& p2)
{
  return fabs(p1.x * p2.y - p1.y * p2.x) > std::numeric_limits<T>::epsilon();
}

#endif /* POINT_2D_HPP */
