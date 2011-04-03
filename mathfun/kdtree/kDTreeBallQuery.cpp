# include <mex.h>

#include <list>
#include <map>
#include <vector>

#include <vector.hpp>
#include <KDTree.hpp>

template <unsigned K>
static void dispatch(int n, int m, double *x_ptr, double *c_ptr, double *d_ptr, int nlhs, mxArray *plhs[])
{
  // Read parameters
	typename KDTree<K, double>::points_type X;

	typename KDTree<K, double>::point_type v;

  for (int i = 0; i < n; ++i)
    {
      for (int k = 0; k < K; ++k)
				v[k] = x_ptr[i + (n * k)];
      X.push_back(v);
    }

	typename KDTree<K, double>::points_type C;

  std::vector<double> R;

  for (int i = 0; i < m; ++i)
    {
      for (int k = 0; k < K; ++k)
				v[k] = c_ptr[i + (m * k)];
      C.push_back(v);
      R.push_back(d_ptr[i]);
    }

  // Build kd-tree
  KDTree<K, double> kdtree(X);

  // Compute queries
  std::list<typename KDTree<K, double>::set_type > res_list;
	typename KDTree<K, double>::set_type res;
	
  for (int i = 0; i < m; ++i)
    {
      kdtree.ball_query(C[i], R[i], res);

      res_list.push_back(res);
			
			res.clear();
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
