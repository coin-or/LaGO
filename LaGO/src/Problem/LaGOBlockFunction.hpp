// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOBLOCKFUNCTION_HPP_
#define LAGOBLOCKFUNCTION_HPP_

#include "LaGObase.hpp"
#include "LaGOSymSparseMatrix.hpp"
#include "LaGOSparsity.hpp"
#include "LaGOCurvature.hpp"
#include "LaGOSampleSet.hpp"
#include "LaGOQuadraticUnderestimator.hpp"

namespace LaGO {

/** A block from a larger nonlinear function.
 */
class BlockFunction : public Function {
private:
	vector<int> sparsity;
public:
	/** The indices of the variables this function is defined for.
	 */
	vector<int> indices;
	
	/** The nonquadratic part.
	 */
	SmartPtr<Function> nonquad;
	
	/** The quadratic part.
	 */
	SmartPtr<SymSparseMatrix> quad;
	
	SmartPtr<SparsityGraph> sparsitygraph;

	Curvature curvature;
	/** Minimum eigenvalue of matrix defining quadratic part.
	 * +infinity if not known.
	 */
	double quad_mineig;
	/** Maximum eigenvalue of matrix defining quadratic part.
	 * -infinity if not known.
	 */
	double quad_maxeig;
	
	SampleSet samplepoints;

	/** Underestimators of nonquad.
	 */	
	list<QuadraticUnderestimator> underestimators;

	/** Overestimators of nonquad.
	 */
	list<QuadraticUnderestimator> overestimators;

	BlockFunction(int dim)
	: sparsity(dim), curvature(UNKNOWN), quad_mineig(getInfinity()), quad_maxeig(-getInfinity())
	{ for (int i=0; i<dim; ++i) sparsity[i]=i;
	}

	BlockFunction(const vector<int>& indices_)
	: sparsity(indices_.size()), indices(indices_), curvature(UNKNOWN), quad_mineig(getInfinity()), quad_maxeig(-getInfinity())
	{ for (unsigned int i=0; i<sparsity.size(); ++i) sparsity[i]=i;
	}

	double eval(const DenseVector& x) const;
	
	void gradient(DenseVector& grad, const DenseVector& x) const;
	
	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;
	
	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const;
	
	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const;
	
#ifdef COIN_HAS_FILIB
	bool canIntervalEvaluation() const { return (IsNull(quad) || quad->canIntervalEvaluation()) && (IsNull(nonquad) || nonquad->canIntervalEvaluation()); }
	
	interval<double> eval(const IntervalVector& x) const;

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const;
#endif

	/** Indicates whether the function knows about the variables that appear in it.
	 */
	bool haveSparsity() const { return true; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */ 	
	const vector<int>& getSparsity() const { return sparsity; }
	
	void print(ostream& out) const;
}; // class BlockFunction

} // namespace LaGO

#endif /*LAGOBLOCKFUNCTION_HPP_*/
