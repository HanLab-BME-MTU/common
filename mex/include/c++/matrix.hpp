#ifndef MATRIX_HPP
# define MATRIX_HPP

# include <iostream>

# include <vector.hpp>

template <unsigned n, unsigned m, typename T>
class matrix
{
public:
  static const matrix Id;
	
  matrix(T val = 0)
  {
    for (unsigned i = 0; i < n; ++i)
      for (unsigned j = 0; j < m; ++j)
	data_[i][j] = val;
  }
	
  template <typename U>
  matrix(const matrix<n, m, U> & rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      for (unsigned j = 0; j < m; ++j)
	data_[i][j] = rhs(i, j);		
  }
	
  template <typename U>
  matrix & operator=(const matrix<n, m, U> & rhs)
  {
    for (unsigned i = 0; i < n; ++i)
      for (unsigned j = 0; j < m; ++j)
	data_[i][j] = rhs(i, j);
    return *this;		
  }
	
  const T & operator()(unsigned i, unsigned j) const
  {
    return data_[i][j];
  }
	
  T & operator()(unsigned i, unsigned j)
  {
    return data_[i][j];
  }
	
  unsigned size() const { return n * m; }
	
	matrix<m,n,T> transpose() const
	{
		matrix<m, n, T> tmp;
		for (unsigned i = 0; i < n; ++i)
			for (unsigned j = 0; j < m; ++j)
				tmp(j,i) = data_[i][j];
		return tmp;
	}
	
	double trace() const
	{
		double f = 0;
		for (unsigned i = 0; i < n; ++i)
			f += data_[i][i];
		return f;
	}
	
  static matrix identity();
	
private:
  T data_[n][m];
};

template <unsigned n, unsigned m, typename T>
const matrix<n, m, T> matrix<n, m, T>::Id = matrix<n, m, T>::identity();

template <unsigned n, unsigned m, typename T>
inline
matrix<n, m, T> matrix<n, m, T>::identity()
{
  static matrix<n, n, T> id_;
  static bool flower = true;
  if (flower)
    {
      for (unsigned i = 0; i < n; ++i)
	for (unsigned j = 0; j < m; ++j)
	  id_(i, j) = (i == j);
      flower = false;
    }
  return id_;	
}

template <unsigned n, unsigned m, typename T, typename U>
bool
operator==(matrix<n,m,T>& lhs, const matrix<n,m,U>& rhs)
{
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      if (lhs(i, j) != rhs(i, j))
	return false;
  return true;
}

template <unsigned n, unsigned m, typename T>
inline
matrix<n, m, T>
operator+(const matrix<n,m,T>& lhs, const matrix<n,m,T>& rhs)
{
  matrix<n, m, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      tmp(i, j) = lhs(i, j) + rhs(i, j);
  return tmp;
}

template <unsigned n, unsigned m, typename T>
matrix<n, m, T>
operator-(const matrix<n,m,T>& lhs, const matrix<n,m,T>& rhs)
{
  matrix<n,m, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      tmp(i, j) = lhs(i, j) - rhs(i, j);
  return tmp;
}

template <unsigned n, unsigned m, typename T>
matrix<n, m, T>
operator-(const matrix<n,m,T>& rhs)
{
  matrix<n,m, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; i < m; ++i)
      tmp(i, j) = - rhs(i, j);
  return tmp;
}

template <unsigned n, unsigned o, typename T, unsigned m>
matrix<n, m, T>
operator*(const matrix<n,o,T>& lhs, const matrix<o,m,T>& rhs)
{
  matrix<n,m, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      {
	tmp(i, j) = 0;
	for (unsigned k = 0; k < o; ++k)
	  tmp(i, j) += lhs(i, k) * rhs(k, j);
      }
  return tmp;
}

template <unsigned n, unsigned m, typename T>
vector<n, T>
operator*(const matrix<n,m,T>& lhs, const vector<m,T>& rhs)
{
  vector<n, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    {
      T sum = 0;
      for (unsigned j = 0; j < m; ++j)
	sum += lhs(i, j) * rhs[j];
      tmp[i] = sum;
    }
  return tmp;
}

template <unsigned n, unsigned m, typename T>
matrix<n, m, T>
operator*(const matrix<n,m,T>& lhs, T val)
{
  matrix<n, m, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      tmp(i, j) = lhs(i, j) * val;
  return tmp;
}

template <unsigned n, unsigned m, typename T>
inline
matrix<n,m, T>
operator/(const matrix<n,m,T>& lhs, T val)
{
  matrix<n,m, T> tmp;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      tmp(i,j) = lhs(i, j) / val;
  return tmp;
}

template <unsigned n, unsigned m, typename T, typename O>
std::ostream&
operator<<(std::ostream& ostr, const matrix<n,m,T>& v)
{
  for (unsigned i = 0; i < n; ++i)
    {
      ostr << '[';
      for (unsigned j = 0; j < m; ++j)
	ostr << v(i, j) << (j == m - 1 ? "]" : ", ");
      ostr << std::endl;
    }
  return ostr;
}

#endif
