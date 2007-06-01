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


} // namespace LaGO
