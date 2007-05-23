// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOConfig.h 94 2007-05-21 13:54:40Z stefan $

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOGamsReader.hpp"

using namespace LaGO;

int main(int argc, char** argv) {
	cout << "LaGO " << LAGOVERSION() << endl;
	if (argc<2) {
		cout << "usage: " << argv[0] << " <gams-control-file>" << endl;
		exit(EXIT_FAILURE);
	}

	GamsReader interface;
	SmartPtr<MINLPData> prob(interface.getProblem(argv[1]));

	

	cout << "LaGO finished." << endl;
	return EXIT_SUCCESS;	
} // main
