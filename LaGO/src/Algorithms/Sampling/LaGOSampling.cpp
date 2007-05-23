// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSampling.hpp"

namespace LaGO {

void Sampling::addSet(list<DenseVector>& samplepoints, const set<DenseVector>& pointset, const vector<int>& indices) {
	for (set<DenseVector>::const_iterator it(pointset.begin()); it!=pointset.end(); ++it) {
		samplepoints.push_back(DenseVector());
		samplepoints.back().setBlock(*it, indices);		
	}
}

void Sampling::addVector(list<DenseVector>& samplepoints, const vector<DenseVector>& pointvector, const vector<int>& indices) {
	for (vector<DenseVector>::const_iterator it(pointvector.begin()); it!=pointvector.end(); ++it) {
		samplepoints.push_back(DenseVector());
		samplepoints.back().setBlock(*it, indices);		
	}
}
	
void Sampling::monteCarlo(list<DenseVector>& samplepoints, DenseVector& lower, DenseVector& upper, int nr) {
	for (int i=0; i<nr; ++i) {
		samplepoints.push_back(DenseVector());
		samplepoints.back().setRandom(lower, upper);
	}
}
	
	
} // namespace LaGO
