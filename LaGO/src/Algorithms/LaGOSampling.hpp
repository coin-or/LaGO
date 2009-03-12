// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef SAMPLING_HPP_
#define SAMPLING_HPP_

#include "LaGObase.hpp"
#include "LaGOSampleSet.hpp"

namespace LaGO {
	
class SparsityGraph;

class Sampling {
public:
	void monteCarlo(SampleSet& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr);

	/** Generates points where some of their components are randomly generated.
	 * @param samplepoints Where to add the generated points.
	 * @param basisvector A vector of full size from which the nonrandom components are taken
	 * @param indices The indices of the random components.
	 * @param lower The lower bounds for the random components (size has to be indices.size()).
	 * @param upper The upper bounds for the random components (size has to be indices.size()).
	 * @param nr The number of points to generate.
	 */
	void monteCarlo(SampleSet& samplepoints, const DenseVector& basisvector, const vector<int>& indices, DenseVector& lower, DenseVector& upper, int nr);
	
	int addVertices(SampleSet& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr);
	int addRandomVertices(SampleSet& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr);

	/** Adds minimizer of a function over a box to sample set.
	 * @return iterator of added sample point, if minimization was successfull. samplepoints.end() otherwise.
	 */	
	SampleSet::iterator addMinimizer(SampleSet& samplepoints, const Function& func, const DenseVector& lower, const DenseVector& upper, const SmartPtr<SparsityGraph>& sparsitygraph=NULL, const DenseVector* startpoint=NULL);
};
	
} // namespace LaGO

#endif /*SAMPLING_HPP_*/
