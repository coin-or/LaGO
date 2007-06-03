// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOQUADRATICESTIMATION_HPP_
#define LAGOQUADRATICESTIMATION_HPP_

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOOsiSolver.hpp"
#include "IpIpoptApplication.hpp"

/** Replaces nonquadratic functions by quadratic under- and overestimators.
 */
namespace LaGO {
	
class BoxMinimizationProblem;

class QuadraticEstimation {
private:
	MINLPData& data;
	OsiSolver lp;
	int nr_coeff; // number of coefficients of current quadratic function; set by initLP
	int nr_auxvars; // number of sample points that have an extra column in the LP; set by initLP
	Ipopt::IpoptApplication ipopt;
	
	double eps;
	int iter_max;
	
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
	
	/**
	 * @return The column index for the constraint related to the enforce_tightness constraint.
	 */
	int initLP(BlockFunction& func, const SampleSet::iterator& enforce_tightness);
	SparseVector* constructRow(BlockFunction& func, const DenseVector& point, int samplepoint_col);
	void updateLP(BlockFunction& func, bool as_underestimator);
	/** The heart.
	 */
	SmartPtr<QuadraticFunction> getEstimator(BlockFunction& func, bool as_underestimator, const DenseVector& lower, const DenseVector& upper);
	bool doLocMin(SmartPtr<BoxMinimizationProblem>& prob, BlockFunction& func, const SamplePoint& sample_point, double& f_val, double& viol1, double& viol2, double& scale2, SparseVector*& row, bool do_resolve=false);  
		
public:
	QuadraticEstimation(MINLPData& data_);
	
	void computeEstimators();
	
	void computeEstimators(MINLPData::ObjCon& objcon, bool need_lower, bool need_upper);

	void computeEstimator(BlockFunction& func, bool need_lower, bool need_upper);
}; // class QuadraticEstimation
	
} // namespace LaGO

#endif /*LAGOQUADRATICESTIMATION_HPP_*/
