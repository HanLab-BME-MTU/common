#ifndef WINDOW_HPP
# define WINDOW_HPP

# include <algorithm>

# include <matrix.h> // for mxAssert uses

# include <Point.hpp>

template <typename T, std::size_t Card>
class Window
{
public:
  Window(const Point<T> dps[])
  {
    std::copy(dps, dps + Card, dps_);
  }
  
  const Point<T>& getDPoint(std::size_t i) const
  {
    mxAssert(i < Card, "");

    return dps_[i];
  }

  static std::size_t getCard() { return Card; }
  
private:
  Point<T> dps_[Card];
};

template <typename T>
const Window<T, 4>& neighbC4()
{
  static const Point<T> dps[] = {
    Point<T>( 1,  0),
    Point<T>(-1,  0),
    Point<T>( 0,  1),
    Point<T>( 0, -1)};
  
  static const Window<T, 4> n(dps);
  
  return n;
}

template <typename T>
const Window<T, 8>& neighbC8()
{
  static const Point<T> dps[] = {
    Point<T>(-1, -1), Point<T>( 0, -1),
    Point<T>( 1, -1), Point<T>(-1,  0),
    Point<T>( 1,  0), Point<T>(-1,  1),
    Point<T>( 0,  1), Point<T>( 1,  1)};
  
  static const Window<T, 8> n(dps);
  
  return n;
}

#endif /* WINDOW_HPP */
