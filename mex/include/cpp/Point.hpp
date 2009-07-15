#ifndef POINT_HPP
# define POINT_HPP

# include <iostream>

template <typename T>
struct Point
{
  Point() : x(0), y(0) {}

  Point(T _x, T _y) : x(_x), y(_y) {}

  template <typename U>
  Point& operator=(const Point<U>& rhs)
  {
    x = rhs.x;
    y = rhs.y;
    
    return *this;
  }

  bool operator==(const Point& rhs) const
  {
    return x == rhs.x && y == rhs.y;
  }

  bool operator!=(const Point& rhs) const
  {
    return !this->operator==(rhs);
  }

  Point operator+(const Point& rhs) const
  {
    return Point(x + rhs.x, y + rhs.y);
  }

  Point & operator+=(const Point& rhs)
  {
    x += rhs.x;
    y += rhs.y;

    return *this;
  }

  Point operator-(const Point& rhs) const
  {
    return Point(x - rhs.x, y - rhs.y);
  }

  template <typename U>
  Point<U> operator/(const U& v) const
  {
    return Point<U>(x / v, y / v);
  }

  T x, y;
};

template <typename T>
std::ostream& operator<<(std::ostream& o, const Point<T>& p)
{
  o << '(' << p.x << ',' << p.y << ')';

  return o;
}

template <typename T>
double dist2(const Point<T>& p1, const Point<T>& p2)
{
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;

  return dx * dx + dy * dy;
}
    
template <typename T>
double dist(const Point<T>& p1, const Point<T>& p2)
{
  return sqrt(dist2(p1, p2));
}

template <typename T>
bool areLinearlyIndependent(const Point<T>& p1, const Point<T>& p2)
{
  return fabs(p1.x * p2.y - p1.y * p2.x) > DBL_EPSILON;
}

#endif /* POINT_HPP */
