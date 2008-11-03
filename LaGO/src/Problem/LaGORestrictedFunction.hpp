// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGORESTRICTEDFUNCTION_HPP_
#define LAGORESTRICTEDFUNCTION_HPP_

#include "LaGObase.hpp"

namespace LaGO {

/** A function restricted to a subset of the variables.
 */
class RestrictedFunction : public Function {
private:
	SmartPtr<const Function> f;
	const vector<int>& indices;
	
	mutable DenseVector fullx;

	vector<int> sparsity;
public:
	RestrictedFunction(const SmartPtr<const Function>& f_, const vector<int>& indices_, const DenseVector& refpoint)
	: f(f_), indices(indices_), fullx(refpoint), sparsity(indices_.size())
	{ assert(IsValid(f));
		for (unsigned int i=0; i<sparsity.size(); ++i) sparsity[i]=i;
	} 

	double eval(const DenseVector& x) const {
		fullx.setElementsOfBlock(x, indices);
		return f->eval(fullx);
	}
	
	void gradient(DenseVector& grad, const DenseVector& x) const;
	
	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const;

	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const;

#ifdef COIN_HAS_FILIB
	bool canIntervalEvaluation() const { return f->canIntervalEvaluation(); }
	
	interval<double> eval(const IntervalVector& x) const;

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const;
#endif
	
	bool haveSparsity() const { return true; } 

	virtual const vector<int>& getSparsity() const { return sparsity; }
	
	void print(ostream& out) const { f->print(out); }	
	
}; // class RestrictedFunction	
	
} // namespace LaGO

#endif /*LAGORESTRICTEDFUNCTION_HPP_*/
