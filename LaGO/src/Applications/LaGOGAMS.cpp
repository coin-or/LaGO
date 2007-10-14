// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOGamsReader.hpp"
#include "LaGOGamsSolver.hpp"
#include "LaGOAlgorithm.hpp"
#include "LaGOQuadraticRelaxTest.hpp"

using namespace LaGO;

int main(int argc, char** argv) {
	cout << "LaGO " << LAGOVERSION() << endl;
	if (argc<2) {
		cout << "usage: " << argv[0] << " <gams-control-file>" << endl;
		exit(EXIT_FAILURE);
	}

	GamsReader interface;
	SmartPtr<MINLPData> prob(interface.getProblem(argv[1]));
	
//	Algorithm alg(*prob);
	QuadraticRelaxTest alg(*prob);
	alg.run();
	
	// run NLP solver if we have a startpoint
	if (alg.solution_candidate.size()) {
		DenseVector lower, upper;
		prob->getBox(lower, upper);

		GamsSolver solver(interface, prob->getDiscrVariables());
		solver.setBounds(lower, upper);
		solver.setStartPoint(alg.solution_candidate);
		bool success=solver.callSolver();
		if (success) {
			clog << "GamsSolver returned with model status " << solver.getModelStatus() << " and solver status " << solver.getSolverStatus() << endl;
			clog << "GamsSolver reports optimal value " << solver.getOptimalValue() << endl;
		} else {
			cerr << "Failure in calling GamsSolver." << endl;
		}
	}
	
	// write solution file

	cout << "LaGO finished." << endl;
	return EXIT_SUCCESS;	
} // main
