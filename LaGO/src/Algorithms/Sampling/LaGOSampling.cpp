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

int Sampling::addVertices(list<DenseVector>& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr) {
	DenseVector x(lower);
	long maxnum=0;

	for (int i=0; i<x.getNumElements(); ++i) {
		if (lower[i]>-getInfinity() && upper[i]<getInfinity()) maxnum++;
		else x[i]=getRandom(lower[i], upper[i]);
	}

	if (!maxnum) return 0; // no bounded variables w.r.t. given bounds

	maxnum=(long)pow(2., (double)maxnum); // the number of vertices
	
	if (maxnum<=0) // we have so many vertices, that we cannot generate them all systematically
		return addRandomVertices(samplepoints, lower, upper, nr);

	long dist=maxnum<nr ? 1 : maxnum/nr; // one plus the number of vertices we skip at each iteration    

	long switch_mask=0;
	int added=0;
	while (added<=nr && added<maxnum) {
		switch_mask+=dist;
		long sm=switch_mask^(switch_mask-dist);
		for (int i=0; i<x.getNumElements(); ++i) {
			if (lower(i)>-getInfinity() && upper(i)<getInfinity()) {
				if (sm%2) x[i]=lower[i]+upper[i]-x[i];
				sm/=2;
			}
		}
		samplepoints.push_back(x);
		added++;
	}

	return added; // should be nearly the same as nr
}
	
int Sampling::addRandomVertices(list<DenseVector>& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr) {
	for(int i=0; i<nr; ++i) {
		samplepoints.push_back(lower);
		DenseVector& x(samplepoints.back());
		bool have_boundedvar=false;
		for (int j=0; j<lower.getNumElements(); ++j) {
			if (lower(j)>-getInfinity() && upper(j)<getInfinity()) {
				if (getRandom(0.,1.)>=.5) x[j]=upper[j];
				have_boundedvar=true;
			} else
				x[j]=getRandom(lower[j], upper[j]);
		}
		if (!have_boundedvar) { // if there is no bounded variable, we would be the same as monte carlo, so we do nothing
			samplepoints.pop_back();
			return 0;
		}		
	}
	return nr; 
}

} // namespace LaGO
