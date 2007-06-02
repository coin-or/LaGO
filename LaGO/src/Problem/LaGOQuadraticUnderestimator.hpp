// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOQUADRATICUNDERESTIMATOR_HPP_
#define LAGOQUADRATICUNDERESTIMATOR_HPP_

#include "LaGObase.hpp"

namespace LaGO {

/** Stores a quadratic underestimator xAx+bx+c of a function.
 */
class QuadraticUnderestimator {
public:
	SmartPtr<SymSparseMatrix> A;
	SmartPtr<SparseVector> b;
	double constant;

	QuadraticUnderestimator(const SmartPtr<SymSparseMatrix>& A_, const SmartPtr<SparseVector>& b_, double constant_)
	: A(A_), b(b_), constant(constant_)
	{ } 
}; // class QuadraticUnderestimator 
	
} // namespace LaGO 

#endif /*LAGOQUADRATICUNDERESTIMATOR_HPP_*/
