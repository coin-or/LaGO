// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOQuadraticFunction.hpp 128 2007-08-17 15:35:06Z stefan $

#include "LaGOQuadraticFunction.hpp"

namespace LaGO {

void QuadraticFunction::updateSparsity() {
	sparsity.clear();
	set<int> indices;
	for (int i=0; i<b->getNumElements(); ++i)
		indices.insert(b->getIndices()[i]);
	for (int i=0: i<A->getNumNonzeros(); ++i) {
		indices.insert(A->getRowIndices()[i]);
		indices.insert(A->getColIndices()[i]);
	}
	sparsity.reserve(indices.size());
	for (set<int>::iterator it(indices.begin()); it!=indices.end(); ++it)
		sparsity.push_back(*it);
}


	
} // namespace LaGO
