// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOAlgorithm.hpp 135 2007-10-09 19:37:41Z stefan $

#ifndef LAGOQUADRATICRELAXTEST_HPP_
#define LAGOQUADRATICRELAXTEST_HPP_

#include "LaGObase.hpp"

#include "LaGODecomposition.hpp"
#include "LaGOCurvatureCheck.hpp"
#include "LaGOConstraintPropagation.hpp"
#include "LaGOQuadraticEstimation.hpp"
#include "LaGOConvexification.hpp"

namespace LaGO {

class QuadraticOrConvexApproximation;
	
class QuadraticRelaxTest {
private:
	MINLPData& data;

	Decomposition decomp;
	CurvatureCheck curvcheck;
	ConstraintPropagation conprob;
	BoxReductionGuessing guess;
	QuadraticEstimation quadest;
	Convexification convexify;
	
	void preprocessing();

	bool solveMINLPRelax(SmartPtr<QuadraticOrConvexApproximation> quad);
	bool solveNLPRelax(SmartPtr<QuadraticOrConvexApproximation> quad);
	
public:
	QuadraticRelaxTest(MINLPData& data_);

	void run();
	
}; // class QuadraticRelaxTest	
	
} // namespace LaGO

#endif /*LAGOQUADRATICRELAXTEST_HPP_*/
