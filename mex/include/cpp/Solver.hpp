#ifndef   	SOLVER_HPP
# define   	SOLVER_HPP

# include <cmath>

# include <matrix.h> // for mxAssert uses

class _Solver
{
public:
  int nRoots() const { return nRoots_; }

  double root(unsigned i) const
  {
    mxAssert(i < 4, "");

    return roots_[i];
  }

protected:
  static double roots_[4];
  static int nRoots_;
};

class Solver : public _Solver
{
public:
  Solver(double prec = DBL_EPSILON) : prec_(prec) {}

  // Resolve a linear equation of form: aX + b = 0 in R
  void operator()(double a, double b)
  {
    nRoots_ = 0;

    if (fabs(a) >= prec_)
      {
	roots_[0] = - b / a;
	nRoots_ = 1;
      }
  }

  // Resolve a quadratic equation of form: aX2 + bX + c = 0 in R
  // ref: "Numerical Recipes in C" p. 183--184
  void operator()(double a, double b, double c)
  {
    nRoots_ = 0;

    if (fabs(a) < prec_)
      {
	this->operator()(b, c);
	return;
      }

    double delta = b * b - 4 * a * c;

    if (fabs(delta) < prec_)
      {
	roots_[0] = - b / (2 * a);
	nRoots_ = 1;
      }
    else
      if (delta > 0)
	{
	  double q = - 0.5 * (b + (b < 0 ? -1 : 1) * sqrt(delta));

	  roots_[0] = q / a;
	  roots_[1] = c / q;
	  nRoots_ = 2;
	}
  }

  // Resolve a cubic equation of form:
  // aX3 + bX2 + cX + d = 0 in R
  void operator()(double a, double b, double c, double d)
  {
    nRoots_ = 0;

    if (fabs(a) < prec_)
      {
	this->operator()(b, c, d);
	return;
      }

    cubic_(b / a, c / a, d / a);
  }

  // Resolve a 4th order equation of form: 
  // aX4 + bX3 + cX2 + dX + e = 0 in R
  void operator()(double a, double b, double c, double d, double e)
  {
    nRoots_ = 0;

    if (fabs(a) < prec_)
      {
	this->operator()(b, c, d, e);
	return;
      }

    quartic_(b / a, c / a, d / a, e / a);
  }

private:
  // ref: "Numerical Recipes in C" p. 184--185
  void cubic_(double a, double b, double c)
  {
    double q = (a * a - 3 * b) / 9;
    double r = (2 * a * a * a - 9 * a * b + 27 * c) / 54;

    double r2 = r * r;
    double q3 = q * q * q;
    double a_3 = a / 3;

    if (r2 < q3)
      {
	double theta = acos(r / sqrt(q3));
	  
	roots_[0] = - 2 * sqrt(q) * cos(theta / 3) - a_3;
	roots_[1] = - 2 * sqrt(q) * cos((theta + 2 * M_PI) / 3) - a_3;
	roots_[2] = - 2 * sqrt(q) * cos((theta - 2 * M_PI) / 3) - a_3;

	nRoots_ = 3;
      }
    else
      {
	double a = - (r < 0 ? -1 : 1) * pow(fabs(r) + sqrt(r2 - q3), 1 / 3.0);

	double b = fabs(a) < prec_ ? 0.0 : q / a;

	roots_[0] = a + b - a_3;

	nRoots_ = 1;
      }
  }

  // ref: http://mathworld.wolfram.com/QuarticEquation.html
  void quartic_(double a3, double a2, double a1, double a0)
  {
    cubic_(-a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);

    double y1 = roots_[0];

    nRoots_ = 0;

    double deltaR = 0.25 * a3 * a3 - a2 + y1;

    if (deltaR < 0)
      return;

    double R = sqrt(deltaR);

    double D, E;
    double deltaD, deltaE;

    if (fabs(R) < prec_)
      {
	deltaD = 0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0);
	deltaE = 0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0);

	if (deltaD >= 0)
	  {
	    D = sqrt(deltaD);
	    roots_[0] = - 0.25 * a3 + 0.5 * D;
	    roots_[1] = - 0.25 * a3 - 0.5 * D;
	    nRoots_ = 2;
	  }

	if (deltaE >= 0)
	  {
	    E = sqrt(deltaE);
	    roots_[nRoots_] = - 0.25 * a3 + 0.5 * E;
	    roots_[nRoots_ + 1] = - 0.25 * a3 - 0.5 * E;
	    nRoots_ += 2;
	  }
      }
    else
      {
	deltaD = 0.75 * a3 * a3 - R * R - 2 * a2 + (a3 * a2 - 2 * a1 - 0.25 * a3 * a3 * a3) / R;

	deltaE = 0.75 * a3 * a3 - R * R - 2 * a2 - (a3 * a2 - 2 * a1 - 0.25 * a3 * a3 * a3) / R;

	if (deltaD >= 0)
	  {
	    D = sqrt(deltaD);
	    roots_[0] = - 0.25 * a3 + 0.5 * R + 0.5 * D;
	    roots_[1] = - 0.25 * a3 + 0.5 * R - 0.5 * D;
	    nRoots_ = 2;
	  }

	if (deltaE >= 0)
	  {
	    E = sqrt(deltaE);
	    roots_[nRoots_] = - 0.25 * a3 - 0.5 * R + 0.5 * E;
	    roots_[nRoots_ + 1] = - 0.25 * a3 - 0.5 * R - 0.5 * E;
	    nRoots_ += 2;
	  }
      }
  }

private:
  double prec_;
};

double _Solver::roots_[4] = {0, 0, 0, 0};
int _Solver::nRoots_ = -1;

#endif	    /* !SOLVER_HPP */
