#ifndef IMAGE_HPP
# define IMAGE_HPP

# include <iostream>
# include <cmath>
# include <cassert>

# include <point_2d.hpp>

# include <mex.h>

namespace detail
{
  template <typename T>
  void allocate_data(T* & data,
		     T** & idx_data,
		     int width,
		     int height,
		     int margin)
  {
    int width_eff = width + (margin << 1);
    int height_eff = height + (margin << 1);

    data = new T[width_eff * height_eff];

    idx_data = new T*[width_eff];

    T* ptr = data + margin;

    for (int x = 0; x < width_eff; ++x)
      {
	idx_data[x] = ptr;
	ptr += height_eff;
      }

    idx_data += margin;
  }

  template <typename T>
  void desallocate_data(T* & data, T** & idx_data, int margin)
  {
    idx_data -= margin;
    delete[] idx_data;
    idx_data = 0;
    delete[] data;
    data = 0;
  }
}

// Important note: image class decribes an array of 2-dimensional
// data which is column-wise stored in memory. Consquently, although
// in Matlab, accessing to the point at the ith row and jth column is
// given by I(i, j), elements in image class are accessible via
// I(j, i) or more traditionnaly I(x, y). We choose this since matlab
// array are also column-wise in memory and the copy between Matlab
// arrays and images class are therefore faster to perform.

template <typename T>
class image
{
public:
  image() :
    width_(0), height_(0), margin_(0), data_(0), idx_data_(0)
  {
  }

  image(int width, int height, int margin = 1) :
    width_(width), height_(height), margin_(margin), data_(0), idx_data_(0)
  {
    detail::allocate_data(data_, idx_data_, width_, height_, margin_);
  }

  // src is a pointer to an array of width * height elements
  // src needs to be column-wise stored image.
  image(T* src, int width, int height, int margin = 1) :
    width_(width), height_(height), margin_(margin), data_(0), idx_data_(0)
  {
    detail::allocate_data(data_, idx_data_, width_, height_, margin_);

    T* ptr = src;

    for (int x = 0; x < width_; ++x)
      {
	std::memcpy(idx_data_[x], ptr, sizeof(T) * height_);
	ptr += height_;
      }
  }

  image(const image & rhs) :
    width_(rhs.width()), height_(rhs.height()), margin_(rhs.margin()),
    data_(0), idx_data_(0)
  {
    detail::allocate_data(data_, idx_data_, width_, height_, margin_);
    
    int width_eff = width_ + (margin_ << 1);
    int height_eff = height_ + (margin_ << 1);

    std::memcpy(data_, rhs.data(), sizeof(T) * width_eff * height_eff);
  }

  ~image()
  {
    detail::desallocate_data(data_, idx_data_, margin_);
  }

  image& operator=(const image& rhs)
  {
    if (&rhs == this)
      return *this;
    
    if (rhs.width() != width_ ||
	rhs.height() != height_ ||
	rhs.margin() != margin_)
      {
	detail::desallocate_data(data_, idx_data_, margin_);
	width_ = rhs.width();
	height_ = rhs.height();
	margin_ = rhs.margin();
	detail::allocate_data(data_, idx_data_, width_, height_, margin_);
      }

    int width_eff = width_ + (margin_ << 1);
    int height_eff = height_ + (margin_ << 1);

    std::memcpy(data_, rhs.data(), sizeof(T) * width_eff * height_eff);

    return *this;
  }

public:
  int width() const { return width_; }
  int height() const { return height_; }
  int margin() const { return margin_; }

  const T* data() const { return data_; }

  const T& operator[](const point_2d<int> & p) const
  {
    assert(contains_large(p.x, p.y));

    return idx_data_[p.x][p.y];
  }

  T& operator[](const point_2d<int> & p)
  {
    assert(contains_large(p.x, p.y));

    return idx_data_[p.x][p.y];
  }

  const T& operator()(int x, int y) const
  {
    assert(contains_large(x, y));

    return idx_data_[x][y];
  }

  T& operator()(int x, int y)
  {
    assert(contains_large(x, y));

    return idx_data_[x][y];
  }

  T operator()(double x, double y) const
  {
    int i = (int) floor(x);
    int j = (int) floor(y);
    
    assert(contains_large(i, j) && contains_large(i + 1, j + 1));

    double dx = x - i;
    double dy = y - j;
    
    double v00 = idx_data_[i][j];
    double v10 = idx_data_[i + 1][j];
    double v01 = idx_data_[i][j + 1];
    double v11 = idx_data_[i + 1][j + 1];
    
    return (dx * (v11 * dy - v10 * (dy - 1.0)) -
	    (dx - 1.0) * (v01 * dy - v00 * (dy - 1.0)));
  }

public:
  bool contains(int x, int y) const
  {
    return (x >= 0 && x < width_ && y >= 0 && y < height_);
  }
  
  bool contains_large(int x, int y) const
  {
    return (x >= -margin_ && x < width_ + margin_ &&
	    y >= -margin_ && y < height_ + margin_);
  }

public:
  bool operator==(const image & rhs) const
  {
    if (this == &rhs)
      return true;

    if (rhs.width() != width_ ||
	rhs.height() != height_ ||
	rhs.margin() != margin_)
      return false;
	
    bool same = true;

    int width_eff = width_ + (margin_ << 1);
    int height_eff = height_ + (margin_ << 1);

    const T* rhs_data = rhs.data();

    for (int i = 0; i < width_eff * height_eff; ++i)
      same &= (data_[i] == rhs_data[i]);

    return same;
  }

  bool operator!=(const image & rhs) const
  {
    return !(this->operator==(rhs));
  }

public:
  void border_replicate(int new_margin) const
  {
    if (new_margin != margin_)
      const_cast<image<T>*>(this)->border_reallocate_and_copy_(new_margin);

    const int imax = width_ - 1;
    const int jmax = height_ - 1;

    for (int i = - margin_; i; ++i)
      for (int j = 0; j <= jmax; ++j)
	{
	  idx_data_[i][j] = idx_data_[0][j];
	  idx_data_[imax - i][j] = idx_data_[imax][j];
	}

    for (int i = - margin_; i <= imax + margin_; ++i)
      for (int j = - margin_; j; ++j)
	{
	  idx_data_[i][j] = idx_data_[i][0];
	  idx_data_[i][jmax - j] = idx_data_[i][jmax];
	}
  }
      
  void border_mirror(int new_margin) const
  {
    if (new_margin != margin_)
      const_cast<image<T>*>(this)->border_reallocate_and_copy_(new_margin);

    const int imax = width_ - 1;
    const int jmax = height_ - 1;

    for (int i = - margin_; i; ++i)
      for (int j = 0; j <= jmax; ++j)
	{
	  idx_data_[i][j] = idx_data_[-i][j];
	  idx_data_[imax - i][j] = idx_data_[imax + i][j];
	}

    for (int i = - margin_; i <= imax + margin_; ++i)
      for (int j = - margin_; j; ++j)
	{
	  idx_data_[i][j] = idx_data_[i][-j];
	  idx_data_[i][jmax - j] = idx_data_[i][jmax + j];
	}
  }

  void border_assign(int new_margin, T val) const
  {
    if (new_margin != margin_)
      const_cast<image<T>*>(this)->border_reallocate_and_copy_(new_margin);

    const int imax = width_ - 1;
    const int jmax = height_ - 1;

    for (int i = - margin_; i; ++i)
      for (int j = 0; j <= jmax; ++j)
	{
	  idx_data_[i][j] = val;
	  idx_data_[imax - i][j] = val;
	}

    for (int i = - margin_; i <= imax + margin_; ++i)
      for (int j = - margin_; j; ++j)
	{
	  idx_data_[i][j] = val;
	  idx_data_[i][jmax - j] = val;
	}
  }

public:
  void image2mxArray(mxArray*& dst) const
  {
    if (dst)
      {
	const mwSize* size = mxGetDimensions(dst);

	if (mxGetNumberOfDimensions(dst) != 2 ||
	    size[0] != height_ ||
	    size[1] != width_)
	  {
	    mxDestroyArray(dst);
	    dst = mxCreateDoubleMatrix(height_, width_, mxREAL);
	  }
      }
    else
      dst = mxCreateDoubleMatrix(height_, width_, mxREAL);
    
    double* ptr = mxGetPr(dst);
    
    for (int x = 0; x < width_; ++x)
      {
	double* offset = ptr + x * height_;
	
	for (int y = 0; y < height_; ++y)
	  offset[y] = idx_data_[x][y];
      }
  }

private:
  void border_reallocate_and_copy_(int new_margin)
  {
    T* data = 0;

    T** idx_data = 0;

    detail::allocate_data(data, idx_data, width_, height_, new_margin);

    for (int x = 0; x < width_; ++x)
      std::memcpy(idx_data[x], & idx_data_[x][0], sizeof(T) * height_);
      
    detail::desallocate_data(data_, idx_data_, margin_);

    margin_ = new_margin;
    data_ = data;
    idx_data_ = idx_data;
  }

private:
  int width_, height_, margin_;

  T* data_;
  T** idx_data_;
};

#endif /* IMAGE_HPP */
