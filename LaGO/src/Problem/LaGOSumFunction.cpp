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
