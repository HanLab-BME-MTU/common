#ifndef	GAUSSIANDERIVATIVE1D_HPP
# define GAUSSIANDERIVATIVE1D_HPP

# include <matrix.h> // for mxAssert uses

template <int N>
struct HermitePolynomial
{
  static double res(double x, double inv_v)
  {
    return inv_v * (x * HermitePolynomial<N - 1>::res(x, inv_v) -
		    HermitePolynomial<N - 1>::derivate(x, inv_v));
  }

  static double derivate(double x, double inv_v)
  {
    return N * HermitePolynomial<N - 1>::res(x, inv_v);
  }
};

template <>
struct HermitePolynomial<0>
{
  static double res(double x, double inv_v) { return 1; }

  static double derivate(double x, double inv_v) { return 0; }
};

template <int N>
class GaussianDerivative1D
{
public:
  GaussianDerivative1D(double sigma) :
    sigma_(sigma),
    size_(0),
    data_(0)
  {
    // Note: sigma should be > 0. However, sigma < 1 yields to a very
    // degrated discretization of the Gaussian filter.

    mxAssert(sigma >= 1, "");

    int hside = ((int) ceil(6 * sigma));

    size_ = 2 * hside + 1;

    data_ = new double[size_];

    double inv_v = 1 / (sigma * sigma);

    double c = (N & 1 ? -1 : 1) / (sqrt(2 * M_PI) * sigma);
    
    for (int x = -hside; x <= hside; ++x)
      data_[x + hside] = c * HermitePolynomial<N>::res(x, inv_v) *
	exp(- 0.5 * x * x * inv_v);
  }
  
  ~GaussianDerivative1D()
  {
    delete[] data_;
  }

public:
  double sigma() const { return sigma_; }

  int size() const { return size_; }

  double operator[](int x) const
  {
    mxAssert(x >= 0 && x < size_, "");

    return data_[x];
  }

  double normL1() const
  {
    double res = 0;

    for (int x = 0; x < size_; ++x)
      res += data_[x];

    return res;
  }

  double normL2() const
  {
    double res = 0;

    for (int x = 0; x < size_; ++x)
      res += data_[x] * data_[x];

    return sqrt(res);
  }

private:
  double sigma_;

  int size_;

  double* data_;
};

typedef GaussianDerivative1D<0> Gaussian;

#endif /* !GAUSSIANDERIVATIVE1D_HPP */
