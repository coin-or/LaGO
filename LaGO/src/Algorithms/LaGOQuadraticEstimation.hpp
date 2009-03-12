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
public:
	/** A class for a function and some auxiliary information that is needed to compute a quadratic estimator.  
	 * The reason for the class is that the QuadraticEstimation method should become more independent of LaGO so that it can also be used by other packages which do not want to setup a MINLPData or BlockFunction object. 
	 */ 
	class NonconvexFunction : public ReferencedObject {
	public:
		virtual ~NonconvexFunction() { }
		
		virtual int dim()=0;
		
		virtual SmartPtr<Function> getFunction()=0;
		
		virtual SmartPtr<SparsityGraph> getSparsityGraph()=0;
		
		virtual SampleSet& getSamplePoints()=0;
	};

private:
//	MINLPData& data;
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
	
	class BlockFunctionWrapper : public NonconvexFunction {
	private:
		BlockFunction& blockfunc;
	public:
		BlockFunctionWrapper(BlockFunction& blockfunc_)
		: blockfunc(blockfunc_)
		{ }
		
		int dim() { return blockfunc.indices.size(); } 
		
		SmartPtr<Function> getFunction() { return blockfunc.nonquad; }

		SmartPtr<SparsityGraph> getSparsityGraph() { return blockfunc.sparsitygraph; }
		
		SampleSet& getSamplePoints() { return blockfunc.samplepoints; }		
	};
	
	list<SampleSetItem> sampleset_newsort;

	/**
	 * @return The column index for the constraint related to the enforce_tightness constraint.
	 */
	int initLP(NonconvexFunction& func, const SampleSet::iterator& enforce_tightness);
	/** Constructs a row that can be added to the LP.
	 * @param scale Stores the scaling factor at output.
	 * @param samplepoint_col If not negative, then a coefficient 1 is added for column samplepoint_col.
	 */
	SparseVector* constructRow(NonconvexFunction& func, const DenseVector& point, int samplepoint_col, double& scale);
	void updateLP(bool as_underestimator);
	/** The heart.
	 */
	SmartPtr<QuadraticFunction> getEstimator(NonconvexFunction& func, bool as_underestimator, const DenseVector& lower, const DenseVector& upper);
	bool doLocMin(SmartPtr<BoxMinimizationProblem>& prob, NonconvexFunction& func, const SamplePoint& sample_point, double& f_val, double& viol1, double& viol2, double& scale2, SparseVector*& row, bool do_resolve=false);  

public:
	QuadraticEstimation();
	
	int computeEstimators(MINLPData& data);
	
	/** Computes additional quadratic under- and overestimators for nonquad. terms.
	 * @return The number of computed estimators.
	 */
	int computeImprovingEstimators(MINLPData& data, const DenseVector& refpoint);

	int computeEstimators(MINLPData& data, MINLPData::ObjCon& objcon, bool need_lower, bool need_upper);

	int computeImprovingEstimators(MINLPData& data, MINLPData::ObjCon& objcon, const DenseVector& refpoint, bool need_lower, bool need_upper);

	/** Computes quadratic under- and overestimators of a given nonconvex function.
	 * @param func The nonconvex function.
	 * @param lower Lower bounds on variables.
	 * @param upper Upper bounds on variables.
	 * @param do_lower Indicates whether an underestimator should be computed.
	 * @param do_upper Indicates whether an overestimator should be computed.
	 * @return The quadratic under- and overestimators, if computed. 
	 */
	pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > computeFirstEstimator(NonconvexFunction& func, const DenseVector& lower, const DenseVector& upper, bool do_lower, bool do_upper);

	pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > computeAdditionalEstimator(NonconvexFunction& func, const DenseVector& lower, const DenseVector& upper, SampleSet::iterator enforce_tightness, bool do_lower, bool do_upper);

	/** Computes a quadratic under- or overestimators of a given nonconvex function.
	 * @param func The nonconvex function.
	 * @param lower Lower bounds on variables.
	 * @param upper Upper bounds on variables.
	 * @param refpoint The reference point where the estimator should be exact.
	 * @param as_underestimator Whether an underestimator (=true) or an overestimator (=false) should be computed.
	 * @return The quadratic under- or overestimators, if computed. 
	 */
	SmartPtr<QuadraticFunction> computeEstimator(NonconvexFunction& func, const DenseVector& lower, const DenseVector& upper, const DenseVector& refpoint, bool as_understimator);
	
	void testEstimators(const MINLPData& data) const;
	
	void testEstimator(const Function& orig, const QuadraticFunction& estimator, bool is_underestimator, const DenseVector& lower, const DenseVector& upper) const;   
}; // class QuadraticEstimation
	
} // namespace LaGO

#endif /*LAGOQUADRATICESTIMATION_HPP_*/
