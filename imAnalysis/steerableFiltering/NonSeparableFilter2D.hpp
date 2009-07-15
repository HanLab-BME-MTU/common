#ifndef	NONSEPARABLEFILTER2D_HPP
# define NONSEPARABLEFILTER2D_HPP

# include <matrix.h> // for mxAssert uses

class NonSeparableFilter2D
{
public:
  NonSeparableFilter2D(int width, int height) :
    width_(width),
    height_(height)
  {
    data_ = new double*[width_];

    for (int x = 0; x < width_; ++x)
      data_[x] = new double[height_];
  }
  
  ~NonSeparableFilter2D()
  {
    for (int x = 0; x < width_; ++x)
      delete[] data_[x];
    delete[] data_;
  }

public:
  int width() const { return width_; }
  int height() const { return height_; }

  double operator()(int x, int y) const
  {
    mxAssert(x >= 0 && x < width_, "");
    mxAssert(y >= 0 && y < height_, "");

    return data_[x][y];
  }

  double & operator()(int x, int y)
  {
    mxAssert(x >= 0 && x < width_, "");
    mxAssert(y >= 0 && y < height_, "");

    return data_[x][y];
  }

  double normL1() const
  {
    double res = 0;

    for (int x = 0; x < width_; ++x)
      for (int y = 0; y < height_; ++y)
	res += data_[x][y];

    return res;
  }

  double normL2() const
  {
    double res = 0;

    for (int x = 0; x < width_; ++x)
      for (int y = 0; y < height_; ++y)
	res += data_[x][y] * data_[x][y];

    return sqrt(res);
  }

private:
  int width_, height_;

  double** data_;
};

#endif /* !NONSEPARABLEFILTER2D_HPP */
