// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOQuadraticFunction.hpp 109 2007-06-03 19:31:14Z stefan $

#ifndef LAGOQUADRATICESTIMATOR_HPP_
#define LAGOQUADRATICESTIMATOR_HPP_

#include "LaGObase.hpp"
#include "LaGOQuadraticFunction.hpp"

namespace LaGO {

class QuadraticEstimator : public ReferencedObject {
public:
	SmartPtr<QuadraticFunction> func;

	Curvature curvature;

	/** Coefficients for convexification / concavification term.
	 */	
	DenseVector alpha;
	
	// store also reference point, ...
	
	QuadraticEstimator(const SmartPtr<QuadraticFunction>& func_)
	: func(func_), curvature(UNKNOWN)
	{ }
	
	
};


} // namespace LaGO

#endif /*LAGOQUADRATICESTIMATOR_HPP_*/
