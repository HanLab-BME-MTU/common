#ifndef VECTOR_HPP
# define VECTOR_HPP

template <unsigned n, typename T>
class vector
{
public:
  vector()
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] = 0;
  }
	
  template <typename U>
  vector(const vector<n, U> & rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] = rhs[i];
  }
	
  template <typename U>
  vector & operator=(const vector<n, U> & rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      data_[i] = rhs[i];
    return *this;
  }
	
  const T& operator[](unsigned i) const
  {
    return data_[i];
  }

  T& operator[](unsigned i)
  {
    return data_[i];
  }
	
  unsigned size() const
  {
    return n;
  }
	
  const vector & normalize()
  {
    double n_l2 = 0;
    for (unsigned i = 0; i < n; ++i)
      n_l2 += data_[i] * data_[i];
    n_l2 = sqrt(n_l2);
    for (unsigned i = 0; i < n; ++i)
      data_[i] = T(data_[i] / n_l2);
    return *this;
  }
	
private:
  T data_[n];
};

template <unsigned n, typename T, typename U>
inline
bool operator==(const vector<n, T> & lhs, const vector<n, U> & rhs)
{
  for (unsigned i = 0; i < n; ++i)
    if (lhs[i] != rhs[i])
      return false;
  return true;
}

template <unsigned n, typename T>
inline
vector<n, T> operator+(const vector<n, T> & lhs, const vector<n, T> & rhs)
{
  vector<n, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] + rhs[i];
  return tmp;
}

template <unsigned n, typename T>
inline
vector<n, T> operator-(const vector<n, T> & lhs, const vector<n, T> & rhs)
{
  vector<n, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] - rhs[i];
  return tmp;
}

template <unsigned n, typename T>
inline
vector<n, T> operator*(const vector<n, T> & lhs, const vector<n, T> & rhs)
{
  T tmp = 0;
  for (unsigned i = 0; i < n; ++i)
    tmp += lhs[i] * rhs[i];
  return tmp;
}

template <unsigned n, typename T>
inline
vector<n, T> operator*(const vector<n, T> & lhs, const T & rhs)
{
  vector<n, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] * rhs;
  return tmp;
}

template <unsigned n, typename T>
inline
vector<n, T> operator/(const vector<n, T> & lhs, const T & rhs)
{
  vector<n, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    tmp[i] = lhs[i] / rhs;
  return tmp;
}

template <unsigned n, typename T1, typename T2>
inline
double dist(const vector<n, T1> & lhs, const vector<n, T2> & rhs)
{
  double tmp, d = 0;
  for (unsigned i = 0; i < n; ++i)
    {
      tmp = lhs[i] - rhs[i];
      d += tmp * tmp;
    }
  return sqrt(d);
}

template <unsigned n, typename T>
inline
std::ostream & operator<<(std::ostream & ostr, const vector<n,T> & v)
{
  ostr << '(';
  for (unsigned i = 0; i < n; ++i)
    ostr << v[i] << (i == n - 1 ? ")" : ", ");
  return ostr;
}

#endif
