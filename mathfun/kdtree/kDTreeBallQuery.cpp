# include <mex.h>

#include <list>
#include <map>
#include <vector>

#include <vector.hpp>
#include <kdtree.hpp>

template <unsigned K>
class result_set
{
public:
  typedef typename std::multimap<double, vector<K,double> > result_type;

public:
  result_set()
  {
  }

  void operator()(double key, const vector<K,double> & value)
  {
    m_.insert(std::make_pair(key, value));
  }

  const result_type & set() const { return m_; }

private:
  result_type m_;
};

template <unsigned K>
static void dispatch(int n, int m, double *x_ptr, double *c_ptr, double *d_ptr,
		     int nlhs, mxArray *plhs[])
{
  // Read parameters
  std::vector< vector<K, double> > X;

  vector<K, double> v;

  for (int i = 0; i < n; ++i)
    {
      for (int k = 0; k < K; ++k)
	v[k] = x_ptr[i + (n * k)];
      X.push_back(v);
    }

  std::vector< vector<K, double> > C;
  std::vector<double> R;

  for (int i = 0; i < m; ++i)
    {
      for (int k = 0; k < K; ++k)
	v[k] = c_ptr[i + (m * k)];
      C.push_back(v);
      R.push_back(d_ptr[i]);
    }

  // Build kd-tree
  kdtree<K, double>* tree = new kdtree<K, double>(X.begin(), X.end());

  // Compute queries
  std::list<typename result_set<K>::result_type> res;

  for (int i = 0; i < m; ++i)
    {
      result_set<K> f;

      tree->ball_query(C[i], R[i], f);

      res.push_back(f.set());
    }

  // Write output
  if (nlhs > 0)
    {
    }

  if (nlhs > 1)
    {
    }

  if (nlhs > 2)
    {
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  // Check input/output parameter

  if (nrhs != 3)
    mexErrMsgTxt("Three input arguments required.");

  if (nlhs > 3)
    mexErrMsgTxt("Too many output arguments.");

  int k = mxGetN(prhs[0]);
  int n = mxGetM(prhs[0]);
  int m = mxGetM(prhs[1]);

  if (k != mxGetN(prhs[1]))
    mexErrMsgTxt("X and C must have the same number of columns.");

  if (m != mxGetM(prhs[2]))
    mexErrMsgTxt("C and R must have the same number of rows.");

  double *x_ptr = mxGetPr(prhs[0]);
  double *c_ptr = mxGetPr(prhs[1]);
  double *d_ptr = mxGetPr(prhs[2]);

  switch (k)
    {
    case 2: dispatch<2>(n, m, x_ptr, c_ptr, d_ptr, nlhs, plhs); break;
    case 3: dispatch<3>(n, m, x_ptr, c_ptr, d_ptr, nlhs, plhs); break;
    default: mexErrMsgTxt("Dimension not implemented.");
    }  
}
