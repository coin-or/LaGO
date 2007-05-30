// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOCONSTRAINTPROPAGATION_HPP_
#define LAGOCONSTRAINTPROPAGATION_HPP_

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOBoxReduction.hpp"

namespace LaGO {

class ConstraintPropagation {
private:
	MINLPData& data;
public:
	ConstraintPropagation(MINLPData& data_)
	: data(data_)
	{ }
	
	BoxReductionStatistics reduceBox();

}; // class ConstraintPropagation
	
} // namespace LaGO

#endif /*LAGOCONSTRAINTPROPAGATION_HPP_*/
