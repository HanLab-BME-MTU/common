#ifndef IMAGE_HPP
# define IMAGE_HPP

# include <iostream>
# include <cmath>

# include <matrix.h> // for mxAssert uses

# include <Point.hpp>

# include <mex.h>

namespace detail
{
  template <typename T>
  void allocateData(T* & data, T** & idxData, int width, int height, int margin)
  {
    int widthEff = width + (margin << 1);
    int heightEff = height + (margin << 1);

    data = new T[widthEff * heightEff];

    idxData = new T*[widthEff];

    T* ptr = data + margin;

    for (int x = 0; x < widthEff; ++x)
      {
	idxData[x] = ptr;
	ptr += heightEff;
      }

    idxData += margin;
  }

  template <typename T>
  void desallocateData(T* & data, T** & idxData, int margin)
  {
    idxData -= margin;
    delete[] idxData;
    idxData = 0;
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
class Image
{
public:
  Image() :
    width_(0), height_(0), margin_(0), data_(0), idxData_(0)
  {
  }

  Image(int width, int height, int margin = 1) :
    width_(width), height_(height), margin_(margin), data_(0), idxData_(0)
  {
    detail::allocateData(data_, idxData_, width_, height_, margin_);
  }

  // src is a pointer to an array of width * height elements
  // src needs to be column-wise stored image.
  Image(T* src, int width, int height, int margin = 1) :
    width_(width), height_(height), margin_(margin), data_(0), idxData_(0)
  {
    detail::allocateData(data_, idxData_, width_, height_, margin_);

    T* ptr = src;

    for (int x = 0; x < width_; ++x)
      {
	std::memcpy(idxData_[x], ptr, sizeof(T) * height_);
	ptr += height_;
      }
  }

  Image(const Image & rhs) :
    width_(rhs.width()), height_(rhs.height()), margin_(rhs.margin()),
    data_(0), idxData_(0)
  {
    detail::allocateData(data_, idxData_, width_, height_, margin_);
    
    int widthEff = width_ + (margin_ << 1);
    int heightEff = height_ + (margin_ << 1);

    std::memcpy(data_, rhs.data(), sizeof(T) * widthEff * heightEff);
  }

  ~Image()
  {
    detail::desallocateData(data_, idxData_, margin_);
  }

  Image& operator=(const Image& rhs)
  {
    if (&rhs == this)
      return *this;
    
    if (rhs.width() != width_ ||
	rhs.height() != height_ ||
	rhs.margin() != margin_)
      {
	detail::desallocateData(data_, idxData_, margin_);
	width_ = rhs.width();
	height_ = rhs.height();
	margin_ = rhs.margin();
	detail::allocateData(data_, idxData_, width_, height_, margin_);
      }

    int widthEff = width_ + (margin_ << 1);
    int heightEff = height_ + (margin_ << 1);

    std::memcpy(data_, rhs.data(), sizeof(T) * widthEff * heightEff);

    return *this;
  }

public:
  int width() const { return width_; }
  int height() const { return height_; }
  int margin() const { return margin_; }

  const T* data() const { return data_; }

  const T& operator[](const Point<int> & p) const
  {
    mxAssert(containsLarge(p.x, p.y), "");

    return idxData_[p.x][p.y];
  }

  T& operator[](const Point<int> & p)
  {
    mxAssert(containsLarge(p.x, p.y), "");

    return idxData_[p.x][p.y];
  }

  const T& operator()(int x, int y) const
  {
    mxAssert(containsLarge(x, y), "");

    return idxData_[x][y];
  }

  T& operator()(int x, int y)
  {
    mxAssert(containsLarge(x, y), "");

    return idxData_[x][y];
  }

  T operator()(double x, double y) const
  {
    int i = (int) floor(x);
    int j = (int) floor(y);
    
    mxAssert(containsLarge(i, j) && containsLarge(i + 1, j + 1), "");

    double dx = x - i;
    double dy = y - j;
    
    double v00 = idxData_[i][j];
    double v10 = idxData_[i + 1][j];
    double v01 = idxData_[i][j + 1];
    double v11 = idxData_[i + 1][j + 1];
    
    return (dx * (v11 * dy - v10 * (dy - 1.0)) -
	    (dx - 1.0) * (v01 * dy - v00 * (dy - 1.0)));
  }

public:
  bool contains(int x, int y) const
  {
    return (x >= 0 && x < width_ && y >= 0 && y < height_);
  }
  
  bool containsLarge(int x, int y) const
  {
    return (x >= -margin_ && x < width_ + margin_ &&
	    y >= -margin_ && y < height_ + margin_);
  }

public:
  bool operator==(const Image & rhs) const
  {
    if (this == &rhs)
      return true;

    if (rhs.width() != width_ ||
	rhs.height() != height_ ||
	rhs.margin() != margin_)
      return false;
	
    bool same = true;

    int widthEff = width_ + (margin_ << 1);
    int heightEff = height_ + (margin_ << 1);

    const T* rhs_data = rhs.data();

    for (int i = 0; i < widthEff * heightEff; ++i)
      same &= (data_[i] == rhs_data[i]);

    return same;
  }

  bool operator!=(const Image & rhs) const
  {
    return !(this->operator==(rhs));
  }

public:
  void borderReplicate(int newMargin) const
  {
    if (newMargin != margin_)
      const_cast<Image<T>*>(this)->borderReallocateAndCopy_(newMargin);

    const int imax = width_ - 1;
    const int jmax = height_ - 1;

    for (int i = - margin_; i; ++i)
      for (int j = 0; j <= jmax; ++j)
	{
	  idxData_[i][j] = idxData_[0][j];
	  idxData_[imax - i][j] = idxData_[imax][j];
	}

    for (int i = - margin_; i <= imax + margin_; ++i)
      for (int j = - margin_; j; ++j)
	{
	  idxData_[i][j] = idxData_[i][0];
	  idxData_[i][jmax - j] = idxData_[i][jmax];
	}
  }
      
  void borderMirror(int newMargin) const
  {
    if (newMargin != margin_)
      const_cast<Image<T>*>(this)->borderReallocateAndCopy_(newMargin);

    const int imax = width_ - 1;
    const int jmax = height_ - 1;

    for (int i = - margin_; i; ++i)
      for (int j = 0; j <= jmax; ++j)
	{
	  idxData_[i][j] = idxData_[-i][j];
	  idxData_[imax - i][j] = idxData_[imax + i][j];
	}

    for (int i = - margin_; i <= imax + margin_; ++i)
      for (int j = - margin_; j; ++j)
	{
	  idxData_[i][j] = idxData_[i][-j];
	  idxData_[i][jmax - j] = idxData_[i][jmax + j];
	}
  }

  void borderAssign(int newMargin, T val) const
  {
    if (newMargin != margin_)
      const_cast<Image<T>*>(this)->borderReallocateAndCopy_(newMargin);

    const int imax = width_ - 1;
    const int jmax = height_ - 1;

    for (int i = - margin_; i; ++i)
      for (int j = 0; j <= jmax; ++j)
	{
	  idxData_[i][j] = val;
	  idxData_[imax - i][j] = val;
	}

    for (int i = - margin_; i <= imax + margin_; ++i)
      for (int j = - margin_; j; ++j)
	{
	  idxData_[i][j] = val;
	  idxData_[i][jmax - j] = val;
	}
  }

public:
  void convertToMxArray(mxArray*& dst) const
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
	  offset[y] = idxData_[x][y];
      }
  }

private:
  void borderReallocateAndCopy_(int newMargin)
  {
    T* data = 0;

    T** idxData = 0;

    detail::allocateData(data, idxData, width_, height_, newMargin);

    for (int x = 0; x < width_; ++x)
      std::memcpy(idxData[x], & idxData_[x][0], sizeof(T) * height_);
      
    detail::desallocateData(data_, idxData_, margin_);

    margin_ = newMargin;
    data_ = data;
    idxData_ = idxData;
  }

private:
  int width_, height_, margin_;

  T* data_;
  T** idxData_;
};

#endif /* IMAGE_HPP */
