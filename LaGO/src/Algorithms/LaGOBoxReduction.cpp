// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOBoxReduction.hpp 106 2007-05-31 20:31:22Z stefan $

#include "LaGOBoxReduction.hpp"
#include "LaGOMINLPData.hpp"

namespace LaGO {

ostream& operator<<(ostream& out, const BoxReductionStatistics& statistics) {
	if (statistics.empty_box) out << "empty box!\t";
	out << "avg. reduction: " << statistics.avg_reduction << " max. reduction: " << statistics.max_reduction;
	if (statistics.nr_fixed_var) out << " fixed var.: " << statistics.nr_fixed_var;
	if (statistics.shrinked_integer_var) out << " reduced discrete var.: " << statistics.shrinked_integer_var;
	out << endl;
	return out;	
}

void BoxReductionStatistics::printBox(ostream& out, const MINLPData& data) {
	for (int i=0; i<data.numVariables(); ++i) {
		const MINLPData::Variable& var(data.getVariable(i));
		out << var.getName() << "\t [" << var.getLower() << ", " << var.getUpper() << ']' << endl;
	} 
}


int BoxReductionGuessing::guessBounds() {
	if (bound_is_guessed.size()) bound_is_guessed.clear();
	bound_is_guessed.resize(data.numVariables(), 0);
	
	int guessed=0;
	bool missing_bounds=false;
	double min_lower=-1000.; // so that unknown bounds are set to at least 10000
	double max_upper=1000.;
	for (int i=0; i<data.numVariables(); ++i) {
		const MINLPData::Variable& var(data.getVariable(i));
		
		if (var.getLower()>-getInfinity()) {
			if (var.getLower()<min_lower) min_lower=var.getLower();
		} else if (var.isNonlinear()) missing_bounds=true;
		
		if (var.getUpper()<getInfinity()) {
			if (var.getUpper()>max_upper) max_upper=var.getUpper();
		} else if (var.isNonlinear()) missing_bounds=true;
	}
	
	if (missing_bounds)
		for (int i=0; i<data.numVariables(); ++i) {
			const MINLPData::Variable& var(data.getVariable(i));
			if (!var.isNonlinear()) continue; // skip linear variables
			if (var.getLower()<=-getInfinity() || var.getUpper()>=getInfinity()) ++guessed;
			if (var.getLower()<=-getInfinity()) {
				bound_is_guessed[i]+=1;
				data.var[i].lower=10*min_lower;
			}
			if (var.getUpper()>= getInfinity()) {
				data.var[i].upper=10*max_upper;
				bound_is_guessed[i]+=2;
			}
		}

  return guessed;
}


} // namespace LaGO
