// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSampling.hpp"
#include "LaGOBoxMinimizationProblem.hpp"
#include "IpIpoptApplication.hpp"

namespace LaGO {
	
void Sampling::monteCarlo(SampleSet& samplepoints, DenseVector& lower, DenseVector& upper, int nr) {
	DenseVector point;
	for (int i=0; i<nr; ++i) {
		point.setRandom(lower, upper);
		samplepoints.insert(point);
	}
}

void Sampling::monteCarlo(SampleSet& samplepoints, const DenseVector& basisvector, const vector<int>& indices, DenseVector& lower, DenseVector& upper, int nr) {
	assert((int)indices.size()==lower.size());
	assert(lower.size()==upper.size());
	DenseVector point(basisvector);
	for (int i=0; i<nr; ++i) {
		for (unsigned int j=0; j<indices.size(); ++j) 
			point[indices[j]]=getRandom(lower(j), upper(j));
		samplepoints.insert(point);
	}
}

int Sampling::addVertices(SampleSet& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr) {
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
		samplepoints.insert(x);
		added++;
	}

	return added; // should be nearly the same as nr
}
	
int Sampling::addRandomVertices(SampleSet& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr) {
	DenseVector x(lower.getNumElements());
	for(int i=0; i<nr; ++i) {
		bool have_boundedvar=false;
		for (int j=0; j<lower.getNumElements(); ++j) {
			if (lower(j)>-getInfinity() && upper(j)<getInfinity()) {
				if (getRandom(0.,1.)>=.5) x[j]=upper[j];
				else x[j]=lower[j];
				have_boundedvar=true;
			} else
				x[j]=getRandom(lower[j], upper[j]);
		}
		if (!have_boundedvar) // if there is no bounded variable, we would be the same as monte carlo, so we do nothing
			return 0;
		samplepoints.insert(x);		
	}
	return nr; 
}

SampleSet::iterator Sampling::addMinimizer(SampleSet& samplepoints, const Function& func, const DenseVector& lower, const DenseVector& upper, const SmartPtr<SparsityGraph>& sparsitygraph, const DenseVector* startpoint) {
	SmartPtr<BoxMinimizationProblem> prob(new BoxMinimizationProblem(func, lower, upper, sparsitygraph));
	prob->startpoint=startpoint;
	
	Ipopt::IpoptApplication ipopt;
	ipopt.Options()->SetIntegerValue("print_level", 0);
	ipopt.Initialize("");
//	ipopt.Initialize(); // this reads ipopt.opt
	
	SmartPtr<Ipopt::TNLP> tnlp(GetRawPtr(prob));
	Ipopt::ApplicationReturnStatus ret=ipopt.OptimizeTNLP(tnlp);
	
	switch (ret) {
		case Ipopt::Solve_Succeeded:
		case Ipopt::Solved_To_Acceptable_Level:
		case Ipopt::Maximum_Iterations_Exceeded: {
			SampleSet::iterator it=samplepoints.insert(prob->getSolution()).first;
			it->funcvalue=prob->getOptimalValue();		
			return it;
		}
		default:
			return samplepoints.end();
	}
}

} // namespace LaGO
