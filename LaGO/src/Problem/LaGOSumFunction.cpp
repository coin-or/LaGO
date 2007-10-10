// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOQuadraticOrConvexApproximation.cpp 129 2007-08-23 16:12:56Z stefan $

#include "LaGOSumFunction.hpp"
#include "LaGOSymSparseMatrix.hpp"

namespace LaGO {
	
SumFunction::SumFunction(double a, const SmartPtr<Function>& f, double b, const SmartPtr<Function>& g)
: have_sparsity(true)
{
	if (a && IsValid(f)) {
		functions.push_back(f);
		factors.push_back(a);
		if (!f->haveSparsity()) have_sparsity=false;
	}
	if (b && IsValid(g)) {
		functions.push_back(g);
		factors.push_back(b);
		if (!g->haveSparsity()) have_sparsity=false;
	}
	
	if (have_sparsity) setupSparsity();
}

SumFunction::SumFunction(const list<SmartPtr<Function> >& functions_, const list<double>& factors_)
: functions(functions_), factors(factors_), have_sparsity(true)
{	bool setup_factors=factors.empty(); 
	for (list<SmartPtr<Function> >::iterator it(functions.begin()); it!=functions.end(); ++it) {
		if (setup_factors) factors.push_back(1.);
		if (!(*it)->haveSparsity()) have_sparsity=false;
	}

	if (have_sparsity) setupSparsity();
}

void SumFunction::setupSparsity() {
	if (functions.empty()) return;
	if (++functions.begin()==functions.end()) {
		sparsity=functions.front()->getSparsity();
		return;
	}
	set<int> sparsity_set;
	for (list<SmartPtr<Function> >::iterator it(functions.begin()); it!=functions.end(); ++it) {
		const vector<int>& func_sparsity((*it)->getSparsity());
		for (vector<int>::const_iterator it_(func_sparsity.begin()); it_!=func_sparsity.end(); ++it_)
			sparsity_set.insert(*it_);		
	}
	sparsity.reserve(sparsity_set.size());
	for (set<int>::iterator it(sparsity_set.begin()); it!=sparsity_set.end(); ++it)
		sparsity.push_back(*it);	
}

double SumFunction::eval(const DenseVector& x) const {
	list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
	list<double>::const_iterator it_factor(factors.begin());
	double val=0.;
	
	while (it_func!=functions.end()) {
		val+=*it_factor*(*it_func)->eval(x);
		++it_func; ++it_factor;
	}
	
	return val;
}

void SumFunction::gradient(DenseVector& grad, const DenseVector& x) const {
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

void SumFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
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

void SumFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
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

void SumFunction::fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const {
	list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
	list<double>::const_iterator it_factor(factors.begin());

	(*it_func)->fullHessian(hessian, x);
	hessian.scale(*it_factor);
	++it_func; ++it_factor;
	
	if (it_func!=functions.end()) {
		SymSparseMatrixCreator g_hessian(hessian.getDim());
		do {
			(*it_func)->fullHessian(g_hessian,x);
			hessian.add(*it_factor, g_hessian);
			++it_func; ++it_factor;
		} while(it_func!=functions.end());
	}
}

#ifdef COIN_HAS_FILIB
bool SumFunction::canIntervalEvaluation() const {
	bool can_it=true;
	for (list<SmartPtr<Function> >::const_iterator it(functions.begin()); it!=functions.end(); ++it)
		can_it &= (*it)->canIntervalEvaluation();
	return can_it;
}

interval<double> SumFunction::eval(const IntervalVector& x) const {
	list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
	list<double>::const_iterator it_factor(factors.begin());
	interval<double> val;
	
	while (it_func!=functions.end()) {
		val+=*it_factor*(*it_func)->eval(x);
		++it_func; ++it_factor;
	}

	return val;
}

void SumFunction::evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
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

void SumFunction::print(ostream& out) const {
	out << "SumFunction: " << endl;
	list<SmartPtr<Function> >::const_iterator it_func(functions.begin());
	list<double>::const_iterator it_factor(factors.begin());
	int counter=0;
	while (it_func!=functions.end()) {
		out << counter << ". summand: factor = " << *it_factor << "\t function = " << **it_func << endl;
		++counter;
		++it_func;
		++it_factor;
	}
}

	
} // namespace LaGO
