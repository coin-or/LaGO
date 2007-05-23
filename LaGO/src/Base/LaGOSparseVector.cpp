// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSparseVector.hpp" 

namespace LaGO {

SparseVector::SparseVector(const SparseVectorCreator& creator)
: CoinPackedVector(false)
{ if (creator.empty()) return;
	 
	reserve(creator.size());
	SparseVectorCreator::const_iterator it(creator.begin());
	assert(it->first>=0); // we allow only positive indices
	for(; it!=creator.end(); ++it) {
		if (!it->second) continue; // skip zero elements
		insert(it->first, it->second);
	}	
}

}; // namespace LaGO
