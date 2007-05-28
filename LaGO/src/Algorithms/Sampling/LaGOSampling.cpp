// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSampling.hpp"

namespace LaGO {

void Sampling::addSet(list<DenseVector>& samplepoints, const set<DenseVector>& pointset, const vector<int>& indices) {
	for (set<DenseVector>::const_iterator it(pointset.begin()); it!=pointset.end(); ++it) {
		samplepoints.push_back(DenseVector());
		samplepoints.back().setToBlock(*it, indices);		
	}
}

void Sampling::addVector(list<DenseVector>& samplepoints, const vector<DenseVector>& pointvector, const vector<int>& indices) {
	for (vector<DenseVector>::const_iterator it(pointvector.begin()); it!=pointvector.end(); ++it) {
		samplepoints.push_back(DenseVector());
		samplepoints.back().setToBlock(*it, indices);		
	}
}

void Sampling::addVector(list<DenseVector>& samplepoints, const vector<DenseVector>& pointvector) {
	for (vector<DenseVector>::const_iterator it(pointvector.begin()); it!=pointvector.end(); ++it)
		samplepoints.push_back(*it);
}
	
void Sampling::monteCarlo(list<DenseVector>& samplepoints, DenseVector& lower, DenseVector& upper, int nr) {
	for (int i=0; i<nr; ++i) {
		samplepoints.push_back(DenseVector());
		samplepoints.back().setRandom(lower, upper);
	}
}

void Sampling::monteCarlo(list<DenseVector>& samplepoints, const DenseVector& basisvector, const vector<int>& indices, DenseVector& lower, DenseVector& upper, int nr) {
	assert((int)indices.size()==lower.size());
	assert(lower.size()==upper.size());
	for (int i=0; i<nr; ++i) {
		samplepoints.push_back(basisvector);
		for (unsigned int j=0; j<indices.size(); ++j) 
			samplepoints.back()[indices[j]]=getRandom(lower(j), upper(j));
	}
}

	
	
} // namespace LaGO
