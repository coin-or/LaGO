// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOQUADRATICFUNCTION_HPP_
#define LAGOQUADRATICFUNCTION_HPP_

#include "LaGObase.hpp"

namespace LaGO {

/** Stores a quadratic function xAx+bx+c.
 */
class QuadraticFunction : public Function {
private:
	vector<int> sparsity;
public:
	SmartPtr<SymSparseMatrix> A;
	SmartPtr<SparseVector> b;
	double constant;

	QuadraticFunction(const SmartPtr<SymSparseMatrix>& A_, const SmartPtr<SparseVector>& b_, double constant_)
	: A(A_), b(b_), constant(constant_)
	{ updateSparsity();
	}
	
	void updateSparsity();
	
	double eval(const DenseVector& x) const {
		return A->xAx(x)+x**b+constant;
	}
	
	void gradient(DenseVector& grad, const DenseVector& x) const {
		A->multVector(grad, x, 2.);
		grad+=*b;
	}
	
	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
		value=eval(x);
		gradient(grad,x);
	}

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
		A->multVector(product, factor, 2.);
	}
	
	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const {
		hessian.clear();
		hessian.add(2., *A);
	}
	
#ifdef COIN_HAS_FILIB
	bool canIntervalEvaluation() const { return A->canIntervalEvaluation(); }
	
	interval<double> eval(const IntervalVector& x) const {
		return A->xAx(x)+x**b+constant;
	}

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
		value=eval(x);
		A->multVector(grad, x, 2.);
		grad+=*b; 
	} 
#endif

	/** Indicates whether the function knows about the variables that appear in it.
	 */
	bool haveSparsity() const { return true; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */ 	
	const vector<int>& getSparsity() const { return sparsity; }
	
	void print(ostream& out) const {
		out << "Quadratic Function: c=" << constant << " b=" << *b;
		if (A->getNumNonzeros()) out << endl << "A: " << *A; 
	}
	
}; // class QuadraticFunction 
	
} // namespace LaGO 

#endif /*LAGOQUADRATICFUNCTION_HPP_*/
