// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGORestrictedFunction.hpp"
#include "LaGOSymSparseMatrix.hpp"

namespace LaGO {

void RestrictedFunction::gradient(DenseVector& grad, const DenseVector& x) const {
	fullx.setElementsOfBlock(x, indices);
	DenseVector fullgrad(fullx.getNumElements());
//	clog << "RestrictedFunction::gradient at " << x << endl;

	f->gradient(fullgrad, fullx);

	grad.setToBlock(fullgrad, indices);	
}
	
void RestrictedFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
	fullx.setElementsOfBlock(x, indices);
	DenseVector fullgrad(fullx.getNumElements());

	f->evalAndGradient(value, fullgrad, fullx);
	
	grad.setToBlock(fullgrad, indices);	
}

void RestrictedFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
	fullx.setElementsOfBlock(x, indices);
	DenseVector fullfactor(fullx.getNumElements());
	fullfactor.setElementsOfBlock(factor, indices);
	DenseVector fullproduct(fullx.getNumElements());

	f->hessianVectorProduct(fullproduct, fullx, fullfactor);
	
	product.setToBlock(fullproduct, indices);
}

void RestrictedFunction::fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const {
	fullx.setElementsOfBlock(x, indices);
	SymSparseMatrixCreator fullhess;
	f->fullHessian(fullhess, fullx);
	
	map<int,int> varmap; // reverse of indices mapping
	for (size_t i=0; i<indices.size(); ++i)
		varmap[indices[i]] = i;

	hessian.clear();
	hessian.setDim(x.size());
	
	for (SymSparseMatrixCreator::iterator it(fullhess.begin()); it!=fullhess.end(); ++it) {
		map<int,int>::iterator colit(varmap.find(it->first.first));
		if (colit == varmap.end())
			continue;
		map<int,int>::iterator rowit(varmap.find(it->first.second));
		if (rowit == varmap.end())
			continue;
		hessian.insert(colit->second, rowit->second, it->second);
	}
	
#if 0	
	SymSparseMatrixCreator::iterator it(fullhess.begin());
	for (size_t i=0; i<indices.size(); ++i) {
		assert(i==0 || indices[i-1] < indices[i]); // we require indices to be sorted
		if (it->first.first > indices[i]) continue; // skip indices that have no rows
		while (it!=fullhess.end() && it->first.first < indices[i]) ++it; // skip rows not in restricted function
		if (it==fullhess.end()) break;
		assert(indices[i] == it->first.first);
		for (size_t j=0; j<indices.size(); ++j) {
			if (it->first.second > indices[j]) continue; // skip indices that have no cols
			while (it!=fullhess.end() && it->first.first==indices[i] && it->first.second < indices[j]) ++it; // skip cols not in restricted function
			if (it==fullhess.end()) break;
			if (it->first.first > indices[i]) break;
			assert(indices[i] == it->first.first);
			assert(indices[j] == it->first.second);
			hessian.insert(i, j, it->second);
			cout << indices[i] << ' ' << indices[j] << ": " << it->second << endl;
		}
		if (it==fullhess.end()) break;
	}
#endif
}

#ifdef COIN_HAS_FILIB
interval<double> RestrictedFunction::eval(const IntervalVector& x) const {
	IntervalVector fullintx(fullx, fullx);
	fullintx.setElementsOfBlock(x, indices);
	
	return f->eval(fullintx);
}

void RestrictedFunction::evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
	IntervalVector fullintx(fullx, fullx);
	fullintx.setElementsOfBlock(x, indices);
	IntervalVector fullgrad(fullx.getNumElements());

	f->evalAndGradient(value, fullgrad, fullintx);
	
	grad.setToBlock(fullgrad, indices); 
}
#endif
	
} // namespace LaGO
