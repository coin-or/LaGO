// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOCONSTRAINTPROPAGATION_HPP_
#define LAGOCONSTRAINTPROPAGATION_HPP_

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOBoxReduction.hpp"

// from Cgc
#include "DynNet.h"
// from CoinUtils for CoinTriple
#include "CoinSort.hpp"

namespace LaGO {

class ConstraintPropagation {
private:
	MINLPData& data;
public:
	typedef enum {LOWER=1, UPPER=2, ANY=3} BoundType;

	class NodeInfo {
	public:
		int varindex;

		NodeInfo(int varindex_=-1)
		: varindex(varindex_)
		{ }
	}; // NodeInfo
	
	class EdgeInfo {
	public:
		/** A list of (constraint-index, bound type, linear coefficient)
		 */ 
		list<CoinTriple<int, BoundType, double> > coninfo; 

		EdgeInfo() { }

		EdgeInfo(int con_nr, BoundType which_bound, double coeff) {
			coninfo.push_back(CoinTriple<int, BoundType, double>(con_nr, which_bound, coeff));
		}
		
		void add(int con_nr, BoundType which_bound, double coeff) {
			coninfo.push_back(CoinTriple<int, BoundType, double>(con_nr, which_bound, coeff));
		}
	}; // EdgeInfo
	
	typedef Cgc::DynNet<NodeInfo, EdgeInfo> DependencyGraph;

private:
	DependencyGraph depgraph;

	BoxReductionStatistics run(DenseVector& newlow, DenseVector& newup, const DenseVector& oldlow, const DenseVector& oldup, set<pair<const DependencyGraph::Node*, BoundType> >& nodeset);
	bool evalConstraint(interval<double>& val, const MINLPData::Constraint& con, IntervalVector& box, int index);
	
public:
	int print_level;
	int funceval_limit;
	double min_impr;

	ConstraintPropagation(MINLPData& data_)
	: data(data_), depgraph(data_.numVariables(), 3*data_.numVariables()), print_level(1),
	  funceval_limit(10000), min_impr(0.01)
	{ }

	void initDependencyGraph();
	
	BoxReductionStatistics reduceBox();

	BoxReductionStatistics run(DenseVector& newlow, DenseVector& newup, const DenseVector& oldlow, const DenseVector& oldup);

}; // class ConstraintPropagation

} // namespace LaGO

#endif /*LAGOCONSTRAINTPROPAGATION_HPP_*/
