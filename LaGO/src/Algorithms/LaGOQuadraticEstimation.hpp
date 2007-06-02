// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOQUADRATICESTIMATION_HPP_
#define LAGOQUADRATICESTIMATION_HPP_

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOOsiSolver.hpp"

/** Replaces nonquadratic functions by quadratic under- and overestimators.
 */
namespace LaGO {

class QuadraticEstimation {
private:
	MINLPData& data;
	OsiSolver lp;
	double eps;
	
	class SampleSetItem {
	public:
		const SamplePoint& sample_point;
		int colnr;
		int rownr;
		/** Maximal absolute value of the LP coefficients created from this sample point.
		 */  
		double maxcoeff;
		
		SampleSetItem(const SamplePoint& sample_point_, int colnr_, int rownr_, double maxcoeff_)
		: sample_point(sample_point_), colnr(colnr_), rownr(rownr_), maxcoeff(maxcoeff_)
		{ } 		
	};
	
	list<SampleSetItem> sampleset_newsort;
	
	void initLP(BlockFunction& func, const SampleSet::iterator& enforce_tightness);
	SparseVector* constructRow(BlockFunction& func, const SamplePoint& point, int samplepoint_col);
	void updateLP(BlockFunction& func, bool as_underestimator);
	
public:
	QuadraticEstimation(MINLPData& data_)
	: data(data_), eps(1E-4)
	{ }
	
	void computeEstimators();
	
	void computeEstimators(MINLPData::ObjCon& objcon, bool need_lower, bool need_upper);

	void computeEstimator(BlockFunction& func, bool need_lower, bool need_upper);
	
}; // class QuadraticEstimation
	
} // namespace LaGO

#endif /*LAGOQUADRATICESTIMATION_HPP_*/
