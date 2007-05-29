// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOGamsReader.hpp"

#include "LaGODecomposition.hpp"
#include "LaGOCurvatureCheck.hpp"

using namespace LaGO;

int main(int argc, char** argv) {
	cout << "LaGO " << LAGOVERSION() << endl;
	if (argc<2) {
		cout << "usage: " << argv[0] << " <gams-control-file>" << endl;
		exit(EXIT_FAILURE);
	}

	GamsReader interface;
	SmartPtr<MINLPData> prob(interface.getProblem(argv[1]));
	
	Decomposition decomp(*prob);
	decomp.decompose();
	
	CurvatureCheck curv(*prob);
	curv.computeCurvature();

	cout << *prob;

	cout << "LaGO finished." << endl;
	return EXIT_SUCCESS;	
} // main
