// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef SAMPLING_HPP_
#define SAMPLING_HPP_

#include "LaGObase.hpp"

namespace LaGO {

class Sampling {
public:
	void addSet(list<DenseVector>& samplepoints, const set<DenseVector>& pointset, const vector<int>& indices);
	
	void addVector(list<DenseVector>& samplepoints, const vector<DenseVector>& pointvector, const vector<int>& indices);

	void monteCarlo(list<DenseVector>& samplepoints, DenseVector& lower, DenseVector& upper, int nr);
	
	
};
	
} // namespace LaGO

#endif /*SAMPLING_HPP_*/
