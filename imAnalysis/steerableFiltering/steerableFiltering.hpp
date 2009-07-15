#ifndef	STEERABLEFILTERING_HPP
# define STEERABLEFILTERING_HPP

# include <Image.hpp>
# include <convolution.hpp>
# include <GaussianDerivative1D.hpp>
# include <NonSeparableFilter2D.hpp>
# include <Solver.hpp>

inline void sndGaussianFiltering(const Image<double> & I,
				 double sigmaPSF,
				 Image<double> & R,
				 Image<double> & T)
{
  Image<double> fyy0(I.width(), I.height());
  Image<double> fyy60(I.width(), I.height());
  Image<double> fyy120(I.width(), I.height());

  //////////////////////////////////////////////
  // Compute f * gyy0 (separable convolution) //
  //////////////////////////////////////////////

  Gaussian g0(sigmaPSF);
  GaussianDerivative1D<2> g2(sigmaPSF);

  convolve(I, g0, g2, fyy0);

  ///////////////////////////////////////////////////
  // Compute f * gyy60 (non-separable convolution) //
  ///////////////////////////////////////////////////

  int hside = (int) ceil(6 * sigmaPSF);

  NonSeparableFilter2D gyy60(2 * hside + 1, 2 * hside + 1);

  double theta = M_PI / 3;
  double ct = cos(theta);
  double st = sin(theta);
  double s2 = sigmaPSF * sigmaPSF;
  double s6 = s2 * s2 * s2;

  for (int x = -hside; x <= hside; ++x)
    for (int y = -hside; y <= hside; ++y)
      {
	double yy = (-x * st + y * ct);

	gyy60(x + hside, y + hside) = 
	  (yy * yy - s2) / (2 * M_PI * s6) * exp(- (x * x + y * y) / (2 * s2));
      }

  convolve(I, gyy60, fyy60);

  ////////////////////////////////////////////////////
  // Compute f * gyy120 (non-separable convolution) //
  ////////////////////////////////////////////////////

  NonSeparableFilter2D gyy120(2 * hside + 1, 2 * hside + 1);

  theta = 2 * M_PI / 3;
  ct = cos(theta);
  st = sin(theta);

  for (int x = -hside; x <= hside; ++x)
    for (int y = -hside; y <= hside; ++y)
      {
	double yy = (-x * st + y * ct);

	gyy120(x + hside, y + hside) = 
	  (yy * yy - s2) / (2 * M_PI * s6) * exp(- (x * x + y * y) / (2 * s2));
      }

  convolve(I, gyy120, fyy120);

  Solver s;

  const double sqrt3Div2 = sqrt(3 / 2);

  double a, b, c, k1, k2, k3, r;

  for (int i = 0; i < I.width(); ++i)
    for (int j = 0; j < I.height(); ++j)
      {
	a = sqrt3Div2 * (fyy60(i, j) - fyy120(i, j));
	b = 2 * fyy0(i, j) - fyy60(i, j) - fyy120(i, j);
	c = sqrt3Div2 * (fyy120(i, j) - fyy60(i, j));
	      
	s(a, b, c);

	// Note: If there is no real solution, we provide 1 real
	// default solution (x = 0).
	      
	if (s.nRoots() == 0)
	  s(1, 0);
	    
	R(i, j) = -DBL_MAX;
	      
	for (int k = 0; k < s.nRoots(); ++k)
	  {
	    theta = atan(s.root(k));

	    k1 = 1 + 2 * cos(2 * theta);
	    k2 = 1 + 2 * cos(2 * (theta - M_PI / 3));
	    k3 = 1 + 2 * cos(2 * (theta - 2 * M_PI / 3));

	    // We put a - since we want I * -Gyy_theta.

	    r = -(1.0 / 3) * (k1 * fyy0(i, j) +
			      k2 * fyy60(i, j) +
			      k3 * fyy120(i, j));

	    if (r > R(i, j))
	      {
		R(i, j) = r;
		T(i, j) = theta;
	      }
	  }
      }
}

inline void unserFilteringM1(const Image<double> & I,
			     double sigmaPSF,
			     Image<double> & R,
			     Image<double> & T)
{
  Image<double> fx(I.width(), I.height());
  Image<double> fy(I.width(), I.height());

  Gaussian g0(sigmaPSF);
  GaussianDerivative1D<1> g1(sigmaPSF);

  convolve(I, g1, g0, fx);
  convolve(I, g0, g1, fy);
  
  const double a11 = -0.7978845608028653714;

  double r1, r2, theta1, theta2;

  for (int i = 0; i < I.width(); ++i)
    for (int j = 0; j < I.height(); ++j)
      {
	//////////////////////////////////////////////////////////////////
        // Bug: this does not yield to a normally distributed response	//
	// with null mean and predictable variance. However, when one	//
	// of r1 or r2 is discarded and - sin(...) is replaced by +	//
	// sin(...), it yields to a normally distributed response with	//
	// null mean and predictable variance.				//
        //////////////////////////////////////////////////////////////////

 	theta1 = atan2(-fx(i, j), fy(i, j));

  	r1 = a11 * (cos(theta1) * fy(i, j) - sin(theta1) * fx(i, j));

	theta2 = theta1 + M_PI;

  	r2 = a11 * (cos(theta2) * fy(i, j) - sin(theta2) * fx(i, j));

	if (r1 > r2)
	  {
	    R(i, j) = r1;
	    T(i, j) = theta1;
	  }
	else
	  {
	    R(i, j) = r2;
	    T(i, j) = theta2;
	  }
      }
}

inline void unserFilteringM2(const Image<double> & I,
			     double sigmaPSF,
			     Image<double> & R,
			     Image<double> & T)
{
  Image<double> fxx(I.width(), I.height());
  Image<double> fxy(I.width(), I.height());
  Image<double> fyy(I.width(), I.height());

  Gaussian g0(sigmaPSF);
  GaussianDerivative1D<1> g1(sigmaPSF);
  GaussianDerivative1D<2> g2(sigmaPSF);

  convolve(I, g2, g0, fxx);
  convolve(I, g1, g1, fxy);
  convolve(I, g0, g2, fyy);
  
  const double a20 = 0.16286750396763996 * sigmaPSF;	// a20
  const double a22 = -0.4886025119029199 * sigmaPSF;	// a22 

  double b1, b2, b3, a, b, c, theta, ct, st, r;

  Solver s;

  for (int i = 0; i < I.width(); ++i)
    for (int j = 0; j < I.height(); ++j)
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
	      
	if (s.nRoots() == 0)
	  s(1, 0);
	    
	R(i, j) = -DBL_MAX;
	      
	for (int k = 0; k < s.nRoots(); ++k)
	  {
	    theta = atan(s.root(k));

	    ct = cos(theta);
	    st = sin(theta);
		  
	    r = ct * ct * b1 + ct * st * b2 + st * st * b3;

	    if (r > R(i, j))
	      {
		R(i, j) = r;
		T(i, j) = theta;
	      }
	  }
      }
}

#endif /* !STEERABLEFILERING_HPP */
