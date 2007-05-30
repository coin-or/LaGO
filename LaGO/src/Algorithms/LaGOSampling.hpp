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

	void addVector(list<DenseVector>& samplepoints, const vector<DenseVector>& pointvector);

	void monteCarlo(list<DenseVector>& samplepoints, DenseVector& lower, DenseVector& upper, int nr);

	/** Generates points where some of their components are randomly generated.
	 * @param samplepoints Where to add the generated points.
	 * @param basisvector A vector of full size from which the nonrandom components are taken
	 * @param indices The indices of the random components.
	 * @param lower The lower bounds for the random components (size has to be indices.size()).
	 * @param upper The upper bounds for the random components (size has to be indices.size()).
	 * @param nr The number of points to generate.
	 */
	void monteCarlo(list<DenseVector>& samplepoints, const DenseVector& basisvector, const vector<int>& indices, DenseVector& lower, DenseVector& upper, int nr);
	
	int addVertices(list<DenseVector>& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr);
	int addRandomVertices(list<DenseVector>& samplepoints, const DenseVector& lower, const DenseVector& upper, int nr);
};
	
} // namespace LaGO

#endif /*SAMPLING_HPP_*/
