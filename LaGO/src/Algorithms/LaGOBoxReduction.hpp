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

/** Sets unbounded nonlinear variables to a guessed bound.
 */
class BoxReductionGuessing {
private:
	MINLPData& data;
public:
	/** Tells for each variable which bounds has been guessed.
	 * 0 if no bound guessed, 1 if lower has been guessed, 2 if upper has been guessed, and 3 if both have been guessed
	 */
	vector<int> bound_is_guessed;
	bool nr_guessed_bounds;

	BoxReductionGuessing(MINLPData& data_)
	: data(data_), nr_guessed_bounds(0)
	{ }

	/** Guesses bounds of unbounded nonlinear variables. 
	 */	
	void guessBounds();

	/** Info about the number of variables which bounds have been guessed.
	 * @return the number of variables where at least one of its bounds was guessed in the last call to guessBounds().
	 */	
	int numGuessedBounds() { return nr_guessed_bounds; }

};

} // namespace LaGO

#endif /*LAGOBOXREDUCTION_HPP_*/
