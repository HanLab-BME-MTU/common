#ifndef		DISKWINDOW_HPP
# define	DISKWINDOW_HPP

# include <cstring>

# include <matrix.h> // for mxAssert uses

# include <Point.hpp>

template <typename W>
class DiskWindow
{
public:
  DiskWindow(double radius) : radius_(radius)
  {
    mxAssert(radius > 0, "");

    hside_ = (int) rint(radius);

    const double r2 = radius_ * radius_;
    
    for (int x = -hside_; x <= hside_; ++x)
      for (int y = -hside_; y <= hside_; ++y)
	if (x * x + y * y <= r2)
	  dps_.push_back(Point<int>(x, y));

    const int length = 2 * hside_ + 1;

    weights_ = new W*[length];

    for (int x = 0; x < length; ++x)
      weights_[x] = new W[length];
  }

  ~DiskWindow()
  {
    const int length = 2 * hside_ + 1;

    for (int x = 0; x < length; ++x)
      delete[] weights_[x];
    delete[] weights_;
  }

  double radius() const { return radius_; }
  int hside() const { return hside_; }
  std::size_t getCard() const { return dps_.size(); }

  // FIXME: W must be 'castable' into int type.
  void fill(W value)
  {
    const int length = 2 * hside_ + 1;
    
    for (int x = 0; x < length; ++x)
      std::memset(weights_[x], (int) value, sizeof(W) * length);
  }

  bool isInside(const Point<int>& p) const
  {
    const double r2 = radius_ * radius_;

    return p.x * p.x + p.y * p.y <= r2;
  }

  const Point<int>& getDPoint(std::size_t i) const
  {
    mxAssert(i < dps_.size(), "");

    return dps_[i];
  }

  const W & getWeight(std::size_t i) const
  {
    mxAssert(i < dps_.size(), "");
   
    return getWeight(dps_[i]);
  }

  const W & getWeight(const Point<int>& p) const
  {
    return weights_[p.x + hside_][p.y + hside_];
  }

  W & getWeight(std::size_t i)
  {
    mxAssert(i < dps_.size(), "");

    return getWeight(dps_[i]);
  }

  W & getWeight(const Point<int>& p)
  {
    return weights_[p.x + hside_][p.y + hside_];
  }

private:
  double radius_;

  int hside_;

  std::vector<Point<int> > dps_;

  W** weights_;
};

#endif // DISKWINDOW_HPP
