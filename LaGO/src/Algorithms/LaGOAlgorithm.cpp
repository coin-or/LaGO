// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOAlgorithm.hpp"

namespace LaGO {

Algorithm::Algorithm(MINLPData& data_)
: data(data_), decomp(data_), curvcheck(data_), conprob(data_)
{ }

	
void Algorithm::preprocessing() {
	decomp.decompose();
	
	conprob.reduceBox();
	
	curvcheck.computeCurvature();
	

	cout << data;
};
	
void Algorithm::run() {
	preprocessing();
};
	
} // namespace LaGO
