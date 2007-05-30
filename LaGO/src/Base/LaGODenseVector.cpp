// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGODenseVector.hpp"
#include "LaGORandomNumber.hpp"

namespace LaGO {
	
void DenseVector::setRandom(const DenseVector& lower, const DenseVector& upper) {
	assert(lower.getNumElements()==upper.getNumElements());
	resize(lower.getNumElements());
	const double* lower_=lower.getElements();
	const double* upper_=upper.getElements();
	double* me=getElements();
	for (int i=getNumElements(); i>0; --i) {
		*me=getRandom(*lower_, *upper_);
		++lower_;
		++upper_;
		++me; 
	}
} 

	
} // namespace LaGO
