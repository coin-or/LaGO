// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOAlgorithm.hpp"

namespace LaGO {

Algorithm::Algorithm(MINLPData& data_)
: data(data_), decomp(data_), curvcheck(data_), conprob(data_)
{
#ifdef COIN_HAS_FILIB
	filib::fp_traits<double>::setup();
#endif
}

	
void Algorithm::preprocessing() {
	decomp.decompose();
	
//	conprob.print_level=5;
	conprob.initDependencyGraph();
	BoxReductionStatistics statistics=conprob.reduceBox();
	cout << statistics;
	if (statistics.empty_box) {
		cout << "Problem infeasible." << endl;
		return;
	}
//	BoxReductionStatistics::printBox(cout, data);
	
	BoxReductionGuessing guess(data);
	int guess_nr=guess.guessBounds();
	cout << "Variables with guessed bounds: " << guess_nr << endl; 
	
	curvcheck.computeCurvature();
	
	quadest.computeEstimators(data);
	

//	cout << data;
}
	
void Algorithm::run() {
	preprocessing();
}
	
} // namespace LaGO
