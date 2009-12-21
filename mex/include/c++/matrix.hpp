#ifndef MATRIX_HPP
# define MATRIX_HPP

# include <vector.hpp>

template <unsigned n, unsigned m, typename T>
class matrix
{
public:
	static const matrix Id;
	
	matrix(T value = 0)
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
	
	static matrix identity();
	
private:
	T data_[n][m];
};

template <unsigned n, typename T>
const matrix<n, n, T>::Id = matrix<n, n,, T>::identity();

template <unsigned n, typename T>
inline
matrix<n, n, T> matrix<n, n, T>::identity()
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

template <unsigned n, unsigned m, typename T, typename U>
inline
matrix<n, m, mln_trait_op_plus(T,U)>
operator+(const matrix<n,m,T>& lhs, const matrix<n,m,U>& rhs)
{
	matrix<n, m, mln_trait_op_plus(T,U)> tmp;
	for (unsigned i = 0; i < n; ++i)
		for (unsigned j = 0; j < m; ++j)
			tmp(i, j) = lhs(i, j) + rhs(i, j);
	return tmp;
}

template <unsigned n, unsigned m, typename T, typename U>
matrix<n, m, mln_trait_op_minus(T,U)>
operator-(const matrix<n,m,T>& lhs, const matrix<n,m,U>& rhs)
{
	matrix<n,m, mln_trait_op_minus(T,U)> tmp;
	for (unsigned i = 0; i < n; ++i)
		for (unsigned j = 0; j < m; ++j)
			tmp(i, j) = lhs(i, j) - rhs(i, j);
	return tmp;
}

template <unsigned n, unsigned m, typename T>
matrix<n, m, mln_trait_op_uminus(T)>
operator-(const matrix<n,m,T>& lhs)
{
	matrix<n,m, mln_trait_op_uminus(T)> tmp;
	for (unsigned i = 0; i < n; ++i)
		for (unsigned j = 0; i < m; ++i)
			tmp(i, j) = - rhs(i, j);
	return tmp;
}

template <unsigned n, unsigned o, typename T,
unsigned m, typename U>
matrix<n, m, mln_sum_x(T,U)>
operator*(const matrix<n,o,T>& lhs, const matrix<o,m,U>& rhs)
{
	matrix<n,m, mln_sum_x(T,U)> tmp;
	for (unsigned i = 0; i < n; ++i)
		for (unsigned j = 0; j < m; ++j)
		{
			tmp(i, j) = literal::zero;
			for (unsigned k = 0; k < o; ++k)
				tmp(i, j) += lhs(i, k) * rhs(k, j);
		}
	return tmp;
}

template <unsigned n, unsigned m, typename T,
typename U>
vec<n, mln_sum_x(T,U)>
operator*(const matrix<n,m,T>& lhs, const vec<m,U>& rhs)
{
	vec<n, mln_sum_x(T,U)> tmp;
	for (unsigned i = 0; i < n; ++i)
	{
		mln_sum_x(T,U) sum(literal::zero);
		for (unsigned j = 0; j < m; ++j)
			sum += lhs(i, j) * rhs[j];
		tmp[i] = sum;
	}
	return tmp;
}

template <unsigned n, unsigned m, typename T,
typename S>
matrix<n, m, mln_trait_op_times(T,S)>
operator*(const matrix<n,m,T>& lhs, const value::scalar_<S>& s)
{
	S s = s_.to_equiv();
	matrix<n, m, mln_trait_op_times(T,S)> tmp;
	for (unsigned i = 0; i < n; ++i)
		for (unsigned j = 0; j < m; ++j)
			tmp(i, j) = lhs(i, j) * s;
	return tmp;
}

template <unsigned n, unsigned m, typename T, typename S>
inline
matrix<n,m, mln_trait_op_div(T,S)>
operator/(const matrix<n,m,T>& lhs, const value::scalar_<S>& s_)
{
	S s = s_.to_equiv();
	matrix<n,m, mln_trait_op_times(T,S)> tmp;
	for (unsigned i = 0; i < n; ++i)
		for (unsigned j = 0; j < m; ++j)
			tmp(i,j) = lhs(i, j) / s;
	return tmp;
}

template <unsigned n, unsigned m, typename T>
std::ostream&
operator<<(std::ostream& ostr, const matrix<n,m,T>& v)
{
	for (unsigned i = 0; i < n; ++i)
	{
		ostr << '[';
		for (unsigned j = 0; j < m; ++j)
			ostr << debug::format(v(i, j)) << (j == m - 1 ? "]" : ", ");
		ostr << std::endl;
	}
	return ostr;
}

template<unsigned n, unsigned m, typename T>
matrix<m,n,T>
trans(const matrix<n,m,T>& matrix)
{
	matrix<m,n,T> tmp;
	for (unsigned i = 0; i < n; ++i)
        for (unsigned j = 0; j < m; ++j)
			tmp(j,i) = matrix(i,j);
	return tmp;
}

template<unsigned n, typename T> inline
double trace(const matrix<n,n,T>& m)
{
	double f = 0;
	for (unsigned i = 0; i < n; ++i)
        f += m(i,i);
	return f;
}

#endif
