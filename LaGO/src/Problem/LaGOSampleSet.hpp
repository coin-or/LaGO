// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOSAMPLESET_HPP_
#define LAGOSAMPLESET_HPP_

#include "LaGObase.hpp"

namespace LaGO {
	
class SamplePoint {
private:
	double two_norm;
	DenseVector point;
	
	static double tol;

public:
	/** Indicates whether this sample points is a startpoint given to the MINLP.
	 * We might want to consider them as some kind of *important* sample points.
	 */ 
	mutable bool is_startpoint;
	/** Value of nonquadratic function in this point.
	 * Infinity if not computed.
	 */
	mutable double funcvalue;

	SamplePoint(const DenseVector& point_, bool is_startpoint_=false)
	: two_norm(point_.twoNorm()), point(point_), is_startpoint(is_startpoint_), funcvalue(getInfinity())
	{ }
	
	double twoNorm() const { return two_norm; }
	const DenseVector& getPoint() const { return point; }
	// cast to const DenseVector&
	operator const DenseVector&() const { return point; }
	
	bool operator<(const SamplePoint& p) const {
		assert(point.getNumElements()==p.getPoint().getNumElements());
		if (two_norm<p.two_norm-tol*(CoinAbs(two_norm)+1)) return true;
		if (two_norm>p.two_norm+tol*(CoinAbs(two_norm)+1)) return false;
		const double* el=point.getElements();
		const double* pel=p.getPoint().getElements();
		for (int i=point.getNumElements(); i>0; --i, ++el, ++pel)
			if (CoinAbs(*el-*pel)>tol*(CoinAbs(*el)+1)) // points are differnt
				return *el<*pel;
		return false; // points are equal
	}
	
}; // class SamplePoint


class SampleSet : public set<SamplePoint> {
public:
	void addVector(const vector<DenseVector>& pointvector, const vector<int>& indices, bool set_startpoint_flag=false);

	void addVector(const vector<DenseVector>& pointvector, bool set_startpoint_flag=false);

//	SampleSet::iterator insert(const SamplePoint& point) {
//		return insert(pair<double, SamplePoint>(point.twoNorm(), point));
//	}

	iterator eraseAndGetNext(const iterator& it) {
		iterator next(it); ++next;
		erase(it);
		return next;
	}
	
}; // class SampleSet
	
} // namespace LaGO

#endif /*LAGOSAMPLESET_HPP_*/
