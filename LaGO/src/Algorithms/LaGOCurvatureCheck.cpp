// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOCurvatureCheck.hpp"
#include "LaGOSampling.hpp"

namespace LaGO {

double CurvatureCheck::eigenvalue_tolerance=1E-8;

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
		sampling.addVertices(blockfunc.samplepoints, lower, upper, 256);
//		clog << "Sampleset of size " << blockfunc.samplepoints.size() << " created." << endl; 
	}

	double mineig=getInfinity(), maxeig=-getInfinity();
	if (IsNull(blockfunc.nonquad)) { // quadratic function
		if (!blockfunc.quad->computeMinMaxEigenValue(mineig, maxeig)) {
			cerr << "CurvatureCheck::computeCurvature: Error computing eigenvalue for quadratic function." << endl;
			exit(EXIT_FAILURE);
		}
		blockfunc.quad_mineig=mineig;
		blockfunc.quad_maxeig=maxeig;
		mineig*=2; // we actually wanted the eigenvalue of the hessian
		maxeig*=2;
	} else {
		list<DenseVector>::iterator it_sp(blockfunc.samplepoints.begin());
		SymSparseMatrixCreator hessian_creator;			
		int eigval_successes=0; // number of successful eigenvalue computations
		while (it_sp!=blockfunc.samplepoints.end()) {
			try {
				blockfunc.nonquad->fullHessian(hessian_creator, *it_sp);
			} catch (FunctionEvaluationError error) {
				cerr << "CurvatureCheck::computeCurvature: skip sample point due to " << error << endl;
				it_sp=blockfunc.samplepoints.erase(it_sp);
				continue;								
			}
			if (IsValid(blockfunc.quad)) // add quadratic term
				hessian_creator.add(2, *blockfunc.quad);
				
			SymSparseMatrix hessian(hessian_creator);
			double mineig_, maxeig_;
			if (!hessian.computeMinMaxEigenValue(mineig_, maxeig_)) {
				cerr << "CurvatureCheck::computeCurvature: Error computing eigenvalue of hessian of nonquadratic function." << endl;
				continue;
			}
			++eigval_successes;
			if (mineig_<mineig) mineig=mineig_;
			if (maxeig_>maxeig) maxeig=maxeig_;
						
			++it_sp;
		};
		if (blockfunc.samplepoints.empty() || ++blockfunc.samplepoints.begin()==blockfunc.samplepoints.end()) {
			cerr << "Not enough sample points for reliable computation of curvature type." << endl;
			exit(EXIT_FAILURE);   		
		}		
	}
	
	if (mineig<-eigenvalue_tolerance && maxeig>eigenvalue_tolerance) blockfunc.curvature=INDEFINITE;
	else if (mineig<-eigenvalue_tolerance) blockfunc.curvature=CONCAVE; // so max<=eigtol
	else if (maxeig>eigenvalue_tolerance) blockfunc.curvature=CONVEX; // so min>=-eigtol
	else blockfunc.curvature=CONVEXCONCAVE;

//	clog << "mineig: " << mineig << "\t maxeig: " << maxeig << "\t curvature: " << blockfunc.curvature << endl;
}
	
} // namespace LaGO
