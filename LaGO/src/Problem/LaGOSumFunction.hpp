// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOSUMFUNCTION_HPP_
#define LAGOSUMFUNCTION_HPP_

#include "LaGObase.hpp"

namespace LaGO {

class SumFunction : public Function {
private:
	list<SmartPtr<Function> > functions;
	list<double> factors;
	
	vector<int> sparsity;
	bool have_sparsity;
	
	void setupSparsity();

public:
	SumFunction(double a, const SmartPtr<Function>& f, double b, const SmartPtr<Function>& g);
	
	SumFunction(const list<SmartPtr<Function> >& functions_, const list<double>& factors_=list<double>());

	double eval(const DenseVector& x) const;

	void gradient(DenseVector& grad, const DenseVector& x) const;

	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const;

	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const;

#ifdef COIN_HAS_FILIB
	bool canIntervalEvaluation() const;

	interval<double> eval(const IntervalVector& x) const;

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const;
#endif

	/** Indicates whether the function knows about the variables that appear in it.
	 */
	bool haveSparsity() const { return have_sparsity; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */
	virtual const vector<int>& getSparsity() const { return sparsity; }

	void print(ostream& out) const;

}; // class SumFunction

} // namespace LaGO

#endif /*LAGOSUMFUNCTION_HPP_*/
