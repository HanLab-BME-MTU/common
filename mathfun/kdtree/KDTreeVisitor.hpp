#ifndef KDTREEVISITOR_HPP
# define KDTREEVISITOR_HPP

# include <vector>
# include <map>

template <unsigned N, typename T>
class KDNode;

template <unsigned N, typename T>
class KDLeaf;

template <unsigned N, typename T>
class KDTreeVisitor
{
public:
	virtual void visit(const KDNode<N, T> & node) = 0;
	virtual void visit(const KDLeaf<N, T> & leaf) = 0;
};

# include <kdtree.hpp>

template <unsigned N, typename T>
class ClosestPointVisitor : public KDTreeVisitor<N, T>
{
public:
	typedef std::pair<double, unsigned> result_type;
	
public:
	ClosestPointVisitor(const typename KDTraits<N, T>::points_type & points, const vector<N, T> & query)
	: points_(points), query_(query), res_(std::make_pair(0, std::numeric_limits<T>::max())) {}
	
	const result_type & res() const { return res_; }
																				 
	void visit(const KDNode<N, T> & node)
	{
		// top-down
		if (query_[node.dim()] < node.value())
			node.lt()->accept(*this);
		else
			node.gt()->accept(*this);
		
		// bottom-up
		if (query_[node.dim()] < node.value())
		{
			// The closest point has been found in the lt_ child. See
			// whether there isn't any closer point in the gt_ child
			// within min_dist from the closest point.
			if (node.value() - query_[node.dim()] < res_.first)
				node.gt()->accept(*this);
		}
		else
		{
			// The closest point has been found in the gt_ child. See
			// whether there isn't any closer point in the lt_ child
			// within min_dist from the closest point.
			if (query_[node.dim()] - node.value() < res_.first)
				node.lt()->accept(*this);
		}
	}

	void visit(const KDLeaf<N, T> & leaf)
	{
		double d = dist(query_, points_[leaf.index()]);
		
		if (d < res_.first)
			res_ = std::make_pair(d, leaf.index());
	}
	
private:
	const typename KDTraits<N, T>::points_type & points_;
	
	const vector<N, T> & query_;
	
	result_type res_;
};

template <unsigned N, typename T>
class BallQueryVisitor : public KDTreeVisitor<N, T>
{
public:
	typedef std::pair<double, unsigned> pair_type;
	typedef std::multimap<double, unsigned> result_type;
	
public:
	BallQueryVisitor(const typename KDTraits<N, T>::points_type & points, const vector<N, T> & center, double radius)
	: points_(points), center_(center), radius_(radius), res_() {}
	
	const result_type & res() const { return res_; }
	
	void visit(const KDNode<N, T> & node)
	{
		double d = center_[node.dim()] - node.value();

		if (d > 0)
	  {
	    if (d < radius_)
				node.lt()->accept(*this);
			node.gt()->accept(*this);
	  }
		else
	  {
	    if (-d < radius_)
				node.gt()->accept(*this);
			node.lt()->accept(*this);
	  }
	}
	
	void visit(const KDLeaf<N, T> & leaf)
	{
		double d = dist(center_, points_[leaf.index()]);
		
		if (d <= radius_)
			res_.insert(std::make_pair(d, leaf.index()));
	}
	
private:
	const typename KDTraits<N, T>::points_type & points_;
	
	const vector<N, T> & center_;
	
	double radius_;
	
	result_type res_;	
};

#endif /* KDTREEVISITOR_HPP */