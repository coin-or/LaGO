// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOCurvatureCheck.hpp"
#include "LaGOSampling.hpp"

namespace LaGO {

double CurvatureCheck::eigenvalue_tolerance=1E-6;

void CurvatureCheck::computeCurvature() {
	computeCurvature(data.obj);
	for (int c=0; c<data.numConstraints(); ++c)
		computeCurvature(data.con[c]);
}

void CurvatureCheck::computeCurvature(MINLPData::ObjCon& objcon) {
	// assert that function has been decomposed if nonlinear
	assert(IsNull(objcon.origfuncNL) || !objcon.decompfuncNL.empty());
	
	for (vector<SmartPtr<BlockFunction> >::iterator it(objcon.decompfuncNL.begin()); it!=objcon.decompfuncNL.end(); ++it)
		computeCurvature(**it);
}

void CurvatureCheck::computeCurvature(BlockFunction& blockfunc) {
	// need sample points if nonquadratic
	if (IsValid(blockfunc.nonquad) && blockfunc.samplepoints.empty()) {
		DenseVector lower, upper;
		data.getBox(lower, upper, blockfunc.indices);
		Sampling sampling;
		sampling.addVector(blockfunc.samplepoints, data.start_points, blockfunc.indices);
		sampling.monteCarlo(blockfunc.samplepoints, lower, upper, 20);
		sampling.addVertices(blockfunc.samplepoints, lower, upper, 20);
		clog << "Sampleset of size " << blockfunc.samplepoints.size() << " created." << endl; 
	}

	double mineig, maxeig;
	if (IsNull(blockfunc.nonquad)) { // quadratic function
		if (!blockfunc.quad->computeMinMaxEigenValue(mineig, maxeig)) {
			cerr << "CurvatureCheck::computeCurvature: Error computing eigenvalue for quadratic function." << endl;
			exit(EXIT_FAILURE);
		}
	} else {
		// run through sample points
		// compute hessian, add blockfunc.quad, create SparseSymMatrix
		// compute eigenvalues
		
	}

	clog << "mineig: " << mineig << "\t maxeig: " << maxeig << endl;
	// set curvature
}
	
} // namespace LaGO
