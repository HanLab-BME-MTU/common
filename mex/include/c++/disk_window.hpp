#ifndef		DISK_WINDOW_HPP
# define	DISK_WINDOW_HPP

# include <cstring>
# include <cassert>

# include <point_2d.hpp>

template <typename W>
class disk_window
{
public:
  disk_window(double radius) : radius_(radius)
  {
    assert(radius > 0);

    hside_ = (int) rint(radius);

    const double r2 = radius_ * radius_;
    
    for (int x = -hside_; x <= hside_; ++x)
      for (int y = -hside_; y <= hside_; ++y)
	if (x * x + y * y <= r2)
	  dps_.push_back(point_2d<int>(x, y));

    const int length = 2 * hside_ + 1;

    weights_ = new W*[length];

    for (int x = 0; x < length; ++x)
      weights_[x] = new W[length];
  }

  ~disk_window()
  {
    const int length = 2 * hside_ + 1;

    for (int x = 0; x < length; ++x)
      delete[] weights_[x];
    delete[] weights_;
  }

  double radius() const { return radius_; }
  int hside() const { return hside_; }
  std::size_t card() const { return points_.size(); }

  // FIXME: W must be 'castable' into int type.
  void fill(W value)
  {
    const int length = 2 * hside_ + 1;
    
    for (int x = 0; x < length; ++x)
      std::memset(weights_[x], (int) value, sizeof(W) * length);
  }

  bool contains(const point_2d<int>& p) const
  {
    const double r2 = radius_ * radius_;

    return p.x * p.x + p.y * p.y <= r2;
  }

  const point_2d<int>& point(std::size_t i) const
  {
    assert(i < points_.size());

    return points_[i];
  }

  const W & weight(std::size_t i) const
  {
    assert(i < points_.size());
   
    return weight(points_[i]);
  }

  const W & weight(const point_2d<int>& p) const
  {
    return weights_[p.x + hside_][p.y + hside_];
  }

  W & getWeight(std::size_t i)
  {
    assert(i < points_.size());

    return weight(points_[i]);
  }

  W & getWeight(const point_2d<int>& p)
  {
    return weights_[p.x + hside_][p.y + hside_];
  }

private:
  double radius_;

  int hside_;

  std::vector<point_2d<int> > points_;

  W** weights_;
};

#endif // DISK_WINDOW_HPP
