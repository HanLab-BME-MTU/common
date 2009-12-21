#ifndef QUATERNION_HPP
# define QUATERNION_HPP

# include <vector.hpp>

class quaternion : public vector<4, double>
{
public:
  quaternion() {}
	
  quaternion(double s, double x, double y, double z)
  {
    data_[0] = s;
    data_[1] = x;
    data_[2] = y;
    data_[3] = z;
  }
	
  quaternion(double s, const vector<3, double> & v)
  {
    data_[0] = s;
    data_[1] = v[0];
    data_[2] = v[1];
    data_[3] = v[2];
  }
	
  quaternion(const vector<4, double> & v)
  {
    static_cast<super_type>(*this) = v;
  }
	
  const vector<4, double> & to_vector() const
  {
    return static_cast<super_type>(*this);
  }
	
  double s() const { return data_[0]; }
  double & s() { return data_[0]; }
		
  const vector<3, double> & v() const
  {
    return *(const vector<3, double>*)(const void*)(& this->data_[1]);
  }
	
  vector<3, double> & v()
  {
    return *(vector<3, double>*)(void*)(& this->data_[1]);
  }
	
  bool is_pure() const
  {
    return fabs(data_[0]) < std::numeric_limits<double>::epsilon();
  }
	
  quaternion conj() const
  {
    return quaternion(data_[0], -v());
  }
	
  quaternion inverse() const
  {
    double f = this->norm_l2();
    return quaternion(conj() / (f * f));
  }
	
  vector<3, double> rotate(const vector<3, double> & v)
  {
    return ((*this) * quaternion(0, v) * (*this).inverse()).v();
  }
};

inline
quaternion operator*(const quaternion & lhs, const quaternion & rhs)
{
  quaternion tmp(lhs.s() * rhs.s() - lhs.v() * rhs.v(),
		 vprod(lhs.v(),
		       rhs.v()) + lhs.s() * rhs.v() + rhs.s() * lhs.v());
  return tmp;
}

#endif

