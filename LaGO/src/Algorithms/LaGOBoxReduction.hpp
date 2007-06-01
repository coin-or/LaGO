// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOBOXREDUCTION_HPP_
#define LAGOBOXREDUCTION_HPP_

#include "LaGObase.hpp"

namespace LaGO {
	
class MINLPData;

class BoxReductionStatistics {
public:
	bool empty_box; // this indicates infeasiblity
	
	double avg_reduction;
	double max_reduction;
	int nr_fixed_var;
	int shrinked_integer_var;
	
	BoxReductionStatistics()
	: empty_box(false), avg_reduction(0.), max_reduction(0.), nr_fixed_var(0), shrinked_integer_var(0)
	{ }
	
	static void printBox(ostream& out, const MINLPData& data);
	
}; // class BoxReductionStatistics

ostream& operator<<(ostream& out, const BoxReductionStatistics& statistics);

} // namespace LaGO

#endif /*LAGOBOXREDUCTION_HPP_*/
