#ifndef	STEERABLE_FILTERING_HPP
# define STEERABLE_FILTERING_HPP

# include <limits>

# include <image.hpp>
# include <convolution.hpp>
# include <gaussian_derivative_1d.hpp>
# include <solver.hpp>

template <int M>
struct unser_traits
{
  enum { nfilters = (M * (M + 3) / 2) - unser_traits<M - 1>::nfilters };
};

template <>
struct unser_traits<0>
{
  enum { nfilters = 0 };
};

template <int N, int M>
class unser_filtering
{
public:
  static const int nfilters = unser_traits<M>::nfilters;

  unser_filtering(const int size[N]) : res_(size, 2), theta_(size)
  {
    for (int i = 0; i < nfilters; ++i)
      fconvh_[i] = new image<N, double>(size);
  }

  ~unser_filtering()
  {
    for (int i = 0; i < nfilters; ++i)
      delete fconvh_[i];
  }

  const image<N, double> & res() const { return res_; }
  const image<N, double> & theta() const { return theta_; }
  const image<N, double> & fconvh(int i) const { return *fconvh_[i]; }

public:
  void compute(const image<N, double> & ima, double sigma);

private:
  image<N, double> res_, theta_;
  image<N, double> * fconvh_[nfilters];
};

template <>
inline void unser_filtering<2, 1>::compute(const image<2, double> & ima,
					   double sigma)
{
  image<2, double> & fx = *fconvh_[0];
  image<2, double> & fy = *fconvh_[1];

  gaussian g0(sigma);
  gaussian_derivative_1d<1> g1(sigma);

  convolve(ima, g1, g0, fx);
  convolve(ima, g0, g1, fy);
  
  const double a11 = -0.7978845608028653714;

  double r1, r2, t1, t2;

  for (int i = 0; i < ima.width(); ++i)
    for (int j = 0; j < ima.height(); ++j)
      {
 	t1 = atan2(-fx(i, j), fy(i, j));

  	r1 = a11 * (cos(t1) * fy(i, j) - sin(t1) * fx(i, j));

	t2 = t1 + M_PI;

  	r2 = a11 * (cos(t2) * fy(i, j) - sin(t2) * fx(i, j));

	if (r1 > r2)
	  {
	    res_(i, j) = r1;
	    theta_(i, j) = t1;
	  }
	else
	  {
	    res_(i, j) = r2;
	    theta_(i, j) = t2;
	  }
      }
}

template <>
inline void unser_filtering<2, 2>::compute(const image<2, double> & ima,
					   double sigma)
{
  image<2, double> & fxx = *fconvh_[0];
  image<2, double> & fxy = *fconvh_[1];
  image<2, double> & fyy = *fconvh_[2];

  gaussian g0(sigma);
  gaussian_derivative_1d<1> g1(sigma);
  gaussian_derivative_1d<2> g2(sigma);

  convolve(ima, g2, g0, fxx);
  convolve(ima, g1, g1, fxy);
  convolve(ima, g0, g2, fyy);
  
  const double a20 = 0.16286750396763996 * sigma;
  const double a22 = -0.4886025119029199 * sigma;

  double b1, b2, b3, a, b, c, t, ct, st, r;

  solver s(1e-15);

  for (int i = 0; i < ima.width(); ++i)
    for (int j = 0; j < ima.height(); ++j)
      {
	b1 = a22 * fyy(i, j) + a20 * fxx(i, j);
	b2 = 2 * (a20 - a22) * fxy(i, j);
	b3 = a22 * fxx(i, j) + a20 * fyy(i, j);

	a = - b2;
	b = 2 * (b3 - b1);
	c = b2;
	
	if (fabs(a) < s.prec()) a = 0;
	if (fabs(b) < s.prec()) b = 0;
	if (fabs(c) < s.prec()) c = 0;

	// Solve a X^2 + b X + c = 0
	s(a, b, c);

	if (s.nroots() == 0)
	  {
	    // Note: If there is no real solution, we provide 1 real
	    // default solution (theta = 0).
	    
	    res_(i, j) = b1;
	    theta_(i, j) = 0;
	  }
	else
	  if (s.nroots() == 1 && fabs(s.root(0)) < s.prec())
	    {
	      // Note: this is when a and c are null. We test
	      // theta = 0 or pi/2.

	      if (b1 < b3)
		{
		  res_(i, j) = b3;
		  theta_(i, j) = M_PI_2;
		}
	      else
		{
		  res_(i, j) = b1;
		  theta_(i, j) = 0;
		}		  
	    }
	  else
	    {
	      res_(i, j) = -std::numeric_limits<double>::max();
	      
	      for (int k = 0; k < s.nroots(); ++k)
		{
		  t = atan(s.root(k));
		  ct = cos(t);
		  st = sin(t);
		  r = ct * ct * b1 + ct * st * b2 + st * st * b3;

		  if (r > res_(i, j))
		    {
		      res_(i, j) = r;
		      theta_(i, j) = t;
		    }
		}
	    }
      }
}

template <>
inline void unser_filtering<2, 3>::compute(const image<2, double> & ima,
					   double sigma)
{
  // TODO
}

template <>
inline void unser_filtering<2, 4>::compute(const image<2, double> & ima,
					   double sigma)
{
  image<2, double> & fxx = *fconvh_[0];
  image<2, double> & fxy = *fconvh_[1];
  image<2, double> & fyy = *fconvh_[2];
  image<2, double> & fxxxx = *fconvh_[3];
  image<2, double> & fxxxy = *fconvh_[4];
  image<2, double> & fxxyy = *fconvh_[5];
  image<2, double> & fxyyy = *fconvh_[6];
  image<2, double> & fyyyy = *fconvh_[7];

  gaussian g0(sigma);
  gaussian_derivative_1d<1> g1(sigma);
  gaussian_derivative_1d<2> g2(sigma);
  gaussian_derivative_1d<3> g3(sigma);
  gaussian_derivative_1d<4> g4(sigma);

  convolve(ima, g2, g0, fxx);
  convolve(ima, g1, g1, fxy);
  convolve(ima, g0, g2, fyy);
  convolve(ima, g4, g0, fxxxx);
  convolve(ima, g3, g1, fxxxy);
  convolve(ima, g2, g2, fxxyy);
  convolve(ima, g1, g3, fxyyy);
  convolve(ima, g0, g4, fyyyy);
  
  double sigma2 = sigma * sigma;
  double sigma3 = sigma2 * sigma;

  const double a20 = 0.059 * sigma;
  const double a22 = -0.204 * sigma;
  const double a40 = 0.024 * sigma3;
  const double a42 = -0.194 * sigma3;
  const double a44 = 0.063 * sigma3;

  solver s(1e-15);

  double b1, b2, b3, b4, b5, b6, b7, b8, a, b, c, d, e, ct, st, ct2, st2, t, r;

  for (int i = 0; i < ima.width(); ++i)
    for (int j = 0; j < ima.height(); ++j)
      {
	b1 = a22 * fyy(i,j) + a20 * fxx(i,j);
	b2 = 2 * (a20 - a22) * fxy(i,j);
	b3 = a22 * fxx(i,j) + a20 * fyy(i,j);
	b4 = a40 * fxxxx(i,j) + a42 * fxxyy(i,j) + a44 * fyyyy(i,j);
	b5 = 2 * a42 * (fxyyy(i,j) - fxxxy(i,j)) + 4 * (a40 * fxxxy(i,j) -
							a44 * fxyyy(i,j));
	b6 = a42 * (fxxxx(i,j) - 4 * fxxyy(i,j) + fyyyy(i,j)) +
	  6 * a40 * fxxyy(i,j) + 6 * a44 * 
	  fxxyy(i,j);
	b7 = 2 * a42 * (fxxxy(i,j) - fxyyy(i,j)) + 4 * (a40 * fxyyy(i,j) -
							a44 * fxxxy(i,j));
	b8 = a40 * fyyyy(i,j) + a42 * fxxyy(i,j) + a44 * fxxxx(i,j);
		
	a = -b2 - b7;
	b = 2 * (2 * b8 - b6 + b3 - b1);
	c = 3 * (b7 - b5);
	d = 2 * (b6 - 2 * b4 + b3 - b1);
	e = b2 + b5;

	if (fabs(a) < s.prec()) a = 0;
	if (fabs(b) < s.prec()) b = 0;
	if (fabs(c) < s.prec()) c = 0;
	if (fabs(d) < s.prec()) d = 0;
	if (fabs(e) < s.prec()) e = 0;

	// Solve a X^4 + b X^3 + c X^2 + d X + e = 0
	s(a, b, c, d, e);

	if (s.nroots() == 0)
	  {	      
	    // Note: If there is no real solution, we provide 1 real
	    // default solution (theta = 0).
	    
	    res_(i, j) = b1;
	    theta_(i, j) = 0;
	  }
	else
	  if (s.nroots() == 1 && fabs(s.root(0)) < s.prec())
	    {
	      // Note: this is when a, b, c and e are null. We
	      // test theta = 0 or pi/2.
	      
	      if (b3 + b8 < b1 + b4)
		{
		  res_(i, j) = b3 + b8;
		  theta_(i, j) = M_PI_2;
		}
	      else
		{
		  res_(i, j) = b1 + b4;
		  theta_(i, j) = 0;
		}		  
	    }
	  else
	    if (s.nroots() == 3 && c == 0 && e == 0)
	      {
	      }
	    else
	      {
		res_(i, j) = -std::numeric_limits<double>::max();
	      
		for (int k = 0; k < s.nroots(); ++k)
		  {
		    t = atan(s.root(k));

		    ct = cos(t);
		    st = sin(t);

		    ct2 = ct * ct;
		    st2 = st * st;

		    r = ct2 * b1 + ct * st * b2 + st2 * b3 + ct2 * ct2 * b4 +
		      ct2 * ct * st * b5 + ct2 * st2 * b6 + ct * st2 * st * b7
		      + st2 * st2 * b8;
	    
		    if (r > res_(i, j))
		      {
			res_(i, j) = r;
			theta_(i, j) = t;
		      }
		  }
	      }
      }
}

#endif /* !STEERABLE_FILERING_HPP */
