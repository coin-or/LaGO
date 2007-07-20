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
	double a,b;
	SmartPtr<Function> f, g;

public:
	SumFunction(double a_, const SmartPtr<Function>& f_, double b_, const SmartPtr<Function>& g_)
	: a(a_), b(b_), f(f_), g(g_)
	{ }

	double eval(const DenseVector& x) const {
		return a*f->eval(x)+b*g->eval(x);
	}

	void gradient(DenseVector& grad, const DenseVector& x) const {
		f->gradient(grad,x);
		grad.scale(a);
		DenseVector g_grad(grad.getNumElements());
		g->gradient(g_grad,x);
		grad.addVector(b, g_grad);
	}

	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
		f->evalAndGradient(value, grad, x);
		if (a!=1.) { value*=a; grad.scale(a); }
		double g_value;
		DenseVector g_grad(grad.getNumElements());
		g->evalAndGradient(g_value, g_grad, x);
		value+=b*g_value;
		grad.addVector(b, g_grad);
	}

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
		f->hessianVectorProduct(product, x, factor);
		product.scale(a);
		DenseVector g_product(product.getNumElements());
		g->hessianVectorProduct(g_product, x, factor);
		product.addVector(b, g_product);
	}

	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const {
		f->fullHessian(hessian, x);
		hessian.scale(a);
		SymSparseMatrixCreator g_hessian(hessian.getDim());
		g->fullHessian(g_hessian, x);
		hessian.add(b, g_hessian);
	}

#ifdef COIN_HAS_FILIB
	virtual bool canIntervalEvaluation() const { return f->canIntervalEvaluation() && g->canIntervalEvaluation(); }

	interval<double> eval(const IntervalVector& x) const {
		return a*f->eval(x)+b*g->eval(x);
	}

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
		f->evalAndGradient(value, grad, x);
		if (a!=1.) { value*=a; grad.scale(interval<double>(a)); }
		interval<double> g_value;
		IntervalVector g_grad(grad.getNumElements());
		g->evalAndGradient(g_value, g_grad, x);
		value+=b*g_value;
		grad.addVector(b, g_grad);
	}
#endif

	/** Indicates whether the function knows about the variables that appear in it.
	 */
	bool haveSparsity() const { return false; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */
//	virtual const vector<int>& getSparsity() const { throw CoinError("sparsity information not available", "getSparsity()", "Function"); }

	void print(ostream& out) const {
		out << "SumFunction: a=" << a << " b=" << b << endl << "f: " << *f << "g: " << *g;
	}
}; // class SumFunction

} // namespace LaGO

#endif /*LAGOSUMFUNCTION_HPP_*/
