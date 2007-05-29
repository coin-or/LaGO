// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOFunction.hpp"
#include "LaGOSymSparseMatrix.hpp"

namespace LaGO {

void Function::fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const {
	vector<int> nonzeros_dummy;
	if (!haveSparsity()) {
		nonzeros_dummy.reserve(x.getNumElements());
		for (int i=0; i<x.getNumElements(); ++i) nonzeros_dummy.push_back(i);
	}
	const vector<int>& nonzeros(haveSparsity() ? getSparsity() : nonzeros_dummy);
	
	hessian.clear();
	hessian.setDim(x.getNumElements());
	
	DenseVector e(x.getNumElements()); // storage for hessian multiplier
	DenseVector hm(x.getNumElements()); // storage for hessian-vector product
	int i_, j_;
	for (unsigned int i=0; i<nonzeros.size(); ++i) {
		i_=nonzeros[i];
		e[i_]=1.;
		hessianVectorProduct(hm, x, e);
		
		for (unsigned int j=0; j<nonzeros.size(); ++j) {
			j_=nonzeros[j];
			if (j_<i_) continue;
			if (hm[j_]) hessian.insert(i_, j_, hm[j_]);			
		}
		
		e[i_]=0.;
	}
	
//	hessian.set(creator);		
}	

} // namespace LaGO
