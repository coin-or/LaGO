// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGODecomposition.hpp"
#include "LaGOSampling.hpp"

namespace LaGO {

void Decomposition::decompose() {
	decompose(data.obj);
	for (int c=0; c<data.numConstraints(); ++c)
		decompose(data.con[c]); 
}

void Decomposition::decompose(MINLPData::ObjCon& objcon) {
	objcon.decompfuncConstant=objcon.origfuncConstant;
	if (IsNull(objcon.origfuncNL)) {
		objcon.decompfuncLin=objcon.origfuncLin;
		return;
	}
	if (IsValid(objcon.origfuncLin))
		objcon.decompfuncLin=new SparseVector(*objcon.origfuncLin);
	else
		objcon.decompfuncLin=new SparseVector();
	
	
	vector<int> nonzeros_dummy; // if the function does not know its sparsity pattern, we assume a dense function
	if (!objcon.origfuncNL->haveSparsity()) {
		nonzeros_dummy.reserve(data.numVariables());
		for (int i=0; i<data.numVariables(); ++i) nonzeros_dummy.push_back(i);
	}
	const vector<int>& nonzeros(objcon.origfuncNL->haveSparsity() ? objcon.origfuncNL->getSparsity() : nonzeros_dummy);
	 
	DenseVector lower(nonzeros.size());
	DenseVector upper(nonzeros.size());
	data.getBox(lower, upper, nonzeros);
	
	// generate sample points (w.r.t. origfuncNL->sparsity)
	list<DenseVector> samplepoints; 
	Sampling sampling;
	sampling.addVector(samplepoints, data.start_points, nonzeros);
	sampling.monteCarlo(samplepoints, lower, upper, 20); 
	
	// compute sparsity graph

	// find connected components
	
	// distinguish between nonquadratic and quadratic components
	
	// decompose function

	
}

	
} // namespace LaGO
