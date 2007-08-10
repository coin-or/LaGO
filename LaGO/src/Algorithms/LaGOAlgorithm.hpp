// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOALGORITHM_HPP_
#define LAGOALGORITHM_HPP_

#include "LaGObase.hpp"

#include "LaGODecomposition.hpp"
#include "LaGOCurvatureCheck.hpp"
#include "LaGOConstraintPropagation.hpp"
#include "LaGOQuadraticEstimation.hpp"
#include "LaGOConvexification.hpp"

namespace LaGO {

/** LaGOs main algorithm: preprocessing and call of branch-and-cut.
 */
class Algorithm {
private:
	MINLPData& data;
	
	Decomposition decomp;
	CurvatureCheck curvcheck;
	ConstraintPropagation conprob;
	QuadraticEstimation quadest;
	Convexification convexify;
	
	void preprocessing();
	
public:
	Algorithm(MINLPData& data_);
	
	void run();	
}; // class Algorithm
	
} // namespace LaGO

#endif /*LAGOALGORITHM_HPP_*/
