// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOCURVATURECHECK_HPP_
#define LAGOCURVATURECHECK_HPP_

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"

namespace LaGO {

/** Determines the curvature of a function or functions in a MINLP (convex, concave, indefinite).
 */
class CurvatureCheck {
private:
	MINLPData& data;
	
	static double eigenvalue_tolerance;
	
public:
	CurvatureCheck(MINLPData& data_)
	: data(data_)
	{ }
	
	void computeCurvature();
	
	void computeCurvature(MINLPData::ObjCon& objcon);

	void computeCurvature(BlockFunction& func);
	
}; // class CurvatureCheck
	
} // namespace LaGO

#endif /*LAGOCURVATURECHECK_HPP_*/
