#ifndef	STEERABLE_FILTERING_HPP
# define STEERABLE_FILTERING_HPP

# include <limits>

# include <image.hpp>
# include <convolution.hpp>
# include <gaussian_derivative_1d.hpp>
# include <non_separable_filter_2d.hpp>
# include <solver.hpp>

inline void snd_gaussian_filtering(const image<double> & ima,
				   double sigma,
				   image<double> & res,
				   image<double> & theta)
{
  image<double> fyy0(ima.width(), ima.height());
  image<double> fyy60(ima.width(), ima.height());
  image<double> fyy120(ima.width(), ima.height());

  //////////////////////////////////////////////
  // Compute f * gyy0 (separable convolution) //
  //////////////////////////////////////////////

  gaussian g0(sigma);
  gaussian_derivative_1d<2> g2(sigma);

  convolve(ima, g0, g2, fyy0);

  ///////////////////////////////////////////////////
  // Compute f * gyy60 (non-separable convolution) //
  ///////////////////////////////////////////////////

  int hside = (int) ceil(6 * sigma);

  non_separable_filter_2d gyy60(2 * hside + 1, 2 * hside + 1);

  double t = M_PI / 3;
  double ct = cos(t);
  double st = sin(t);
  double s2 = sigma * sigma;
  double s6 = s2 * s2 * s2;

  for (int x = -hside; x <= hside; ++x)
    for (int y = -hside; y <= hside; ++y)
      {
	double yy = (-x * st + y * ct);

	gyy60(x + hside, y + hside) = 
	  (yy * yy - s2) / (2 * M_PI * s6) * exp(- (x * x + y * y) / (2 * s2));
      }

  convolve(ima, gyy60, fyy60);

  ////////////////////////////////////////////////////
  // Compute f * gyy120 (non-separable convolution) //
  ////////////////////////////////////////////////////

  non_separable_filter_2d gyy120(2 * hside + 1, 2 * hside + 1);

  t = 2 * M_PI / 3;
  ct = cos(t);
  st = sin(t);

  for (int x = -hside; x <= hside; ++x)
    for (int y = -hside; y <= hside; ++y)
      {
	double yy = (-x * st + y * ct);

	gyy120(x + hside, y + hside) = 
	  (yy * yy - s2) / (2 * M_PI * s6) * exp(- (x * x + y * y) / (2 * s2));
      }

  convolve(ima, gyy120, fyy120);

  solver s;

  const double sqrt3_div2 = sqrt(3 / 2);

  double a, b, c, k1, k2, k3, r;

  for (int i = 0; i < ima.width(); ++i)
    for (int j = 0; j < ima.height(); ++j)
      {
	a = sqrt3_div2 * (fyy60(i, j) - fyy120(i, j));
	b = 2 * fyy0(i, j) - fyy60(i, j) - fyy120(i, j);
	c = sqrt3_div2 * (fyy120(i, j) - fyy60(i, j));
	      
	s(a, b, c);

	// Note: If there is no real solution, we provide 1 real
	// default solution (x = 0).
	      
	if (s.nroots() == 0)
	  s(1, 0);
	    
	res(i, j) = -std::numeric_limits<double>::max();
	      
	for (int k = 0; k < s.nroots(); ++k)
	  {
	    t = atan(s.root(k));

	    k1 = 1 + 2 * cos(2 * t);
	    k2 = 1 + 2 * cos(2 * (t - M_PI / 3));
	    k3 = 1 + 2 * cos(2 * (t - 2 * M_PI / 3));

	    // We put a - since we want ima * -Gyy_t.

	    r = -(1.0 / 3) * (k1 * fyy0(i, j) +
			      k2 * fyy60(i, j) +
			      k3 * fyy120(i, j));

	    if (r > res(i, j))
	      {
		res(i, j) = r;
		theta(i, j) = t;
	      }
	  }
      }
}

inline void unser_m1_filtering(const image<double> & ima,
			       double sigma,
			       image<double> & res,
			       image<double> & theta)
{
  image<double> fx(ima.width(), ima.height());
  image<double> fy(ima.width(), ima.height());

  gaussian g0(sigma);
  gaussian_derivative_1d<1> g1(sigma);

  convolve(ima, g1, g0, fx);
  convolve(ima, g0, g1, fy);
  
  const double a11 = -0.7978845608028653714;

  double r1, r2, t1, t2;

  for (int i = 0; i < ima.width(); ++i)
    for (int j = 0; j < ima.height(); ++j)
      {
	//////////////////////////////////////////////////////////////////
        // Bug: this does not yield to a normally distributed response	//
	// with null mean and predictable variance. However, when one	//
	// of r1 or r2 is discarded and - sin(...) is replaced by +	//
	// sin(...), it yields to a normally distributed response with	//
	// null mean and predictable variance.				//
        //////////////////////////////////////////////////////////////////

 	t1 = atan2(-fx(i, j), fy(i, j));

  	r1 = a11 * (cos(t1) * fy(i, j) - sin(t1) * fx(i, j));

	t2 = t1 + M_PI;

  	r2 = a11 * (cos(t2) * fy(i, j) - sin(t2) * fx(i, j));

	if (r1 > r2)
	  {
	    res(i, j) = r1;
	    theta(i, j) = t1;
	  }
	else
	  {
	    res(i, j) = r2;
	    theta(i, j) = t2;
	  }
      }
}

inline void unser_m2_filtering(const image<double> & ima,
			       double sigma,
			       image<double> & res,
			       image<double> & theta)
{
  image<double> fxx(ima.width(), ima.height());
  image<double> fxy(ima.width(), ima.height());
  image<double> fyy(ima.width(), ima.height());

  gaussian g0(sigma);
  gaussian_derivative_1d<1> g1(sigma);
  gaussian_derivative_1d<2> g2(sigma);

  convolve(ima, g2, g0, fxx);
  convolve(ima, g1, g1, fxy);
  convolve(ima, g0, g2, fyy);
  
  const double a20 = 0.16286750396763996 * sigma;	// a20
  const double a22 = -0.4886025119029199 * sigma;	// a22 

  double b1, b2, b3, a, b, c, t, ct, st, r;

  solver s;

  for (int i = 0; i < ima.width(); ++i)
    for (int j = 0; j < ima.height(); ++j)
      {
	b1 = a22 * fyy(i, j) + a20 * fxx(i, j);
	b2 = 2 * (a20 - a22) * fxy(i, j);
	b3 = a22 * fxx(i, j) + a20 * fyy(i, j);

	a = - b2;
	b = 2 * (b3 - b1);
	c = b2;
	      
	s(a, b, c);

	// Note: If there is no real solution, we provide 1 real
	// default solution (x = 0).
	      
	if (s.nroots() == 0)
	  s(1, 0);
	    
	res(i, j) = -std::numeric_limits<double>::max();
	      
	for (int k = 0; k < s.nroots(); ++k)
	  {
	    t = atan(s.root(k));

	    ct = cos(t);
	    st = sin(t);
		  
	    r = ct * ct * b1 + ct * st * b2 + st * st * b3;

	    if (r > res(i, j))
	      {
		res(i, j) = r;
		theta(i, j) = t;
	      }
	  }
      }
}

#endif /* !STEERABLE_FILERING_HPP */
