#ifndef KDTREE_HPP
# define KDTREE_HPP

# include <vector>
# include <algorithm>

# include <vector.hpp>

template <unsigned n, typename T>
class cmp
{
public:
  cmp(int dim)
  {
    assert(dim >= 0 && dim < n);

    dim_ = dim;
  }

  bool operator()(const vector<n, T> & p1, const vector<n, T> & p2) const
  {
    return p1[dim_] < p2[dim_];
  }

private:
  int dim_;
};

template <unsigned n, typename T>
class kdtree
{
public:
  typedef typename std::vector< vector<n, T> > points_t;
  typedef typename points_t::iterator points_iterator;

  kdtree(const points_iterator first,
	 const points_iterator last,
	 unsigned dim = 0)
    : pt_(0), lt_(0), gt_(0)
  {
    if (last - first > 1)
      {
	// sort the vector of points along the 'dim' coordinate
	sort(first, last, cmp<n, T>(dim));
	// take the median
	points_iterator median = first + (last - first) / 2;

	pt_ = &(*median);

	lt_ = new kdtree(first, median, (dim + 1) & (n - 1));
	gt_ = new kdtree(median, last, (dim + 1) & (n - 1));
      }
    else
      {
	// this kdtree is a leaf
	pt_ = &(*first);
      }
  }

  ~kdtree()
  {
    if (lt_)
      {
	delete lt_;
	delete gt_;
      }
  }

  double closest_point(const vector<n, T> & query,
		       vector<n, T> & closest,
		       unsigned dim = 0) const
  {
    unsigned next_dim = (dim + 1) & (n - 1);

    if (lt_)
      {
	// top-down
	double min_dist = query[dim] < (*pt_)[dim] ?
	  lt_->closest_point(query, closest, next_dim) :
	  gt_->closest_point(query, closest, next_dim);

	// bottom-up

	if (query[dim] < (*pt_)[dim])
	  {
	    // The closest point has been found in the lt_ child. See
	    // whether there isn't any closer point in the gt_ child
	    // within min_dist from the closest point.
	    if ((*pt_)[dim] - query[dim] < min_dist)
	      {
		vector<n, T> tmp_closest;
		double tmp_dist = gt_->closest_point(query,
						     tmp_closest,
						     next_dim);
		if (tmp_dist < min_dist)
		  {
		    closest = tmp_closest;
		    min_dist =  tmp_dist;
		  }
	      }
	  }
	else
	  {
	    // The closest point has been found in the gt_ child. See
	    // whether there isn't any closer point in the lt_ child
	    // within min_dist from the closest point.
	    if (query[dim] - (*pt_)[dim] < min_dist)
	      {
		vector<n, T> tmp_closest;
		double tmp_dist = lt_->closest_point(query,
						     tmp_closest,
						     next_dim);

		if (tmp_dist < min_dist)
		  {
		    closest = tmp_closest;
		    min_dist = tmp_dist;
		  }
	      }
	  }

	return min_dist;
      }
    else
      {
	closest = *pt_;
	return dist(query, *pt_);
      }
  }
  
  // F is a binary operator
  template <typename F>
  void ball_query(const vector<n,T> & center,
		  const double radius,
		  F & f,
		  unsigned dim = 0) const
  {
    unsigned next_dim = (dim + 1) & (n - 1);

    if (lt_)
      {
	double d = center[dim] - (*pt_)[dim];

	if (d > 0)
	  {
	    if (d < radius)
	      lt_->ball_query(center, radius, f, next_dim);
	    gt_->ball_query(center, radius, f, next_dim);
	  }
	else
	  {
	    if (-d < radius)
	      gt_->ball_query(center, radius, f, next_dim);
	    lt_->ball_query(center, radius, f, next_dim);
	  }
      }
    else
      {
	// this kdtree is a leaf
	double d = dist(center, *pt_);

	if (d <= radius)
	  f(d, *pt_);
      }
  }

private:
  const vector<n, T>* pt_;
  kdtree* lt_, * gt_;
};

#endif /* KDTREE_HPP */
