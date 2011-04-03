#ifndef KDTREE_HPP
# define KDTREE_HPP

# include <vector>
# include <algorithm>

# include <vector.hpp>

template <unsigned N, typename T>
class KDTreeVisitor;

template <unsigned N, typename T>
struct KDTraits
{
	typedef typename std::vector< vector<N, T> > points_type;
	typedef typename std::vector< unsigned >::iterator idx_iterator_type;
	typedef std::pair<unsigned, double> pair_type;
};

template <unsigned N, typename T>
class KDTree
{
public:
	virtual ~KDTree() {}
	
	virtual void accept(KDTreeVisitor<N, T> & v) const = 0;
};

# include <KDTreeVisitor.hpp>

template <unsigned N, typename T>
class KDNode : public KDTree<N, T>
{
public:
	KDNode(KDTree<N, T>* lt, KDTree<N, T>* gt, unsigned dim, double value) : lt_(lt), gt_(gt), dim_(dim), value_(value) {}
	
	~KDNode() {
		delete lt_;
		delete gt_;
	}
	
	void accept(KDTreeVisitor<N, T> & v) const { v.visit(*this); }
	
	unsigned dim() const { return dim_; }
	double value() const { return value_; }
	
	const KDTree<N, T> *lt() const { return lt_; }
	const KDTree<N, T> *gt() const { return gt_; }
	
private:
	KDTree<N, T> *lt_, *gt_;
	unsigned dim_;
	double value_;
};

template <unsigned N, typename T>
class KDLeaf : public KDTree<N, T>
{
public:
	KDLeaf(unsigned index) : index_(index) {}
	
	void accept(KDTreeVisitor<N, T> & v) const { v.visit(*this); }

	unsigned index() const { return index_; }
	
private:
	unsigned index_;
};

template <unsigned N, typename T>
class KDTreeBuilder
{
private:
	struct Incr_
	{
		Incr_() { current = 0; }
		
		int operator()() { return current++; }

		int current;		
	};

	struct Cmp_
	{
	public:
		Cmp_() {}
		
		Cmp_(const typename KDTraits<N, T>::points_type & points, unsigned k) : points_(points), k_(k) {}
		
		bool operator()(unsigned i, unsigned j) const
		{
			return points_[i][k_] < points_[j][k_];
		}
		
	private:
		typename KDTraits<N, T>::points_type points_;
		
		unsigned k_;
	};
	
public:
	KDTreeBuilder(const typename KDTraits<N, T>::points_type & points) : points_(points)
	{
		for (int k = 0; k < N; ++k)
			cmp_[k] = Cmp_(points_, k);
	}
	
	KDTree<N, T>* build() const
	{
		if (points_.empty())
			return NULL;
    
		std::vector<unsigned> indices(points_.size());
		std::generate(indices.begin(), indices.end(), Incr_());
    
		return build_(indices.begin(), indices.end());
	}
	
private:
	KDTree<N, T>* build_(const typename KDTraits<N, T>::idx_iterator_type & first,
											 const typename KDTraits<N, T>::idx_iterator_type & last,
											 unsigned k = 0) const
  {
    if (last - first > 1)
    {
			unsigned k2 = (k + 1) & (N - 1);

      // sort the vector of points along the Kth coordinate
      sort(first, last, cmp_[k]);
      
    	// take the median
      typename KDTraits<N, T>::idx_iterator_type median = first + (last - first) / 2;
      
      double value = points_[*median][k];
      
      KDTree<N, T> * lt = build_(first, median, k2);
      KDTree<N, T> * gt = build_(median, last, k2);
      
      return new KDNode<N, T>(lt, gt, k, value);
    }
    else
      return new KDLeaf<N, T>(*first);
  }
  
private:
  typename KDTraits<N, T>::points_type points_;
	
	Cmp_ cmp_[N];
};

#endif /* KDTREE_HPP */
