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
//	double a,b;
//	SmartPtr<Function> f, g;
	list<SmartPtr<Function> > functions;
	list<double> factors;
	
	vector<int> sparsity;
	bool have_sparsity;
	
	void setupSparsity();

public:
	SumFunction(double a, const SmartPtr<Function>& f, double b, const SmartPtr<Function>& g);
	
	SumFunction(const list<SmartPtr<Function> >& functions_, const list<double>& factors_=list<double>());

	double eval(const DenseVector& x) const {
		list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
		list<double>::const_iterator it_factor(factors.begin());
		double val=0.;
		
		while (it_func!=functions.end()) {
			val+=*it_factor*(*it_func)->eval(x);
			++it_func; ++it_factor;
		}
		
		return val;
	}

	void gradient(DenseVector& grad, const DenseVector& x) const {
		list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
		list<double>::const_iterator it_factor(factors.begin());

		(*it_func)->gradient(grad,x);
		grad.scale(*it_factor);
		++it_func; ++it_factor;
		
		if (it_func!=functions.end()) {
			DenseVector g_grad(grad.getNumElements());
			do {
				(*it_func)->gradient(g_grad,x);
				grad.addVector(*it_factor, g_grad);
				++it_func; ++it_factor;
			} while(it_func!=functions.end());
		}
	}

	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
		list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
		list<double>::const_iterator it_factor(factors.begin());

		(*it_func)->evalAndGradient(value, grad,x);
		grad.scale(*it_factor);
		value*=*it_factor;
		++it_func; ++it_factor;
		
		if (it_func!=functions.end()) {
			DenseVector g_grad(grad.getNumElements());
			double g_value;
			do {
				(*it_func)->evalAndGradient(g_value, g_grad, x);
				value+=*it_factor*g_value;
				grad.addVector(*it_factor, g_grad);
				++it_func; ++it_factor;
			} while(it_func!=functions.end());
		}
	}

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
		list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
		list<double>::const_iterator it_factor(factors.begin());

		(*it_func)->hessianVectorProduct(product, x, factor);
		product.scale(*it_factor);
		++it_func; ++it_factor;
		
		if (it_func!=functions.end()) {
			DenseVector g_product(product.getNumElements());
			do {
				(*it_func)->hessianVectorProduct(g_product, x, factor);
				product.addVector(*it_factor, g_product);
				++it_func; ++it_factor;
			} while(it_func!=functions.end());
		}
	}

	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const;

#ifdef COIN_HAS_FILIB
	virtual bool canIntervalEvaluation() const {
		bool can_it=true;
		for (list<SmartPtr<Function> >::const_iterator it(functions.begin()); it!=functions.end(); ++it)
			can_it &= (*it)->canIntervalEvaluation();
		return can_it;
	}

	interval<double> eval(const IntervalVector& x) const {
		list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
		list<double>::const_iterator it_factor(factors.begin());
		interval<double> val;
		
		while (it_func!=functions.end()) {
			val+=*it_factor*(*it_func)->eval(x);
			++it_func; ++it_factor;
		}

		return val;
	}

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
		list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
		list<double>::const_iterator it_factor(factors.begin());

		(*it_func)->evalAndGradient(value, grad, x);
		if (*it_factor!=1.) {
			grad.scale(*it_factor);
			value*=*it_factor;
		}
		++it_func; ++it_factor;
		
		if (it_func!=functions.end()) {
			IntervalVector g_grad(grad.getNumElements());
			interval<double> g_value;
			do {
				(*it_func)->evalAndGradient(g_value, g_grad, x);
				value+=*it_factor*g_value;
				grad.addVector(*it_factor, g_grad);
				++it_func; ++it_factor;
			} while(it_func!=functions.end());
		}
	}
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
