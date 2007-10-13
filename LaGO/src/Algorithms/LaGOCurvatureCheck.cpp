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
	/*TODO: objcon.origfuncNL could have contained an actually linear function (or at least LaGO::Decomposition thought so). 
	 * In this case objcon.decompfuncNL will be empty. (e.g. for oil2) */ 
	
	for (vector<SmartPtr<BlockFunction> >::iterator it(objcon.decompfuncNL.begin()); it!=objcon.decompfuncNL.end(); ++it)
		computeCurvature(**it);
}

void CurvatureCheck::computeCurvature(BlockFunction& blockfunc) {
	// need sample points if nonquadratic
	if (IsValid(blockfunc.nonquad) && blockfunc.samplepoints.empty()) {
		blockfunc.samplepoints.addVector(data.start_points, blockfunc.indices, true);
		DenseVector lower, upper;
		data.getBox(lower, upper, blockfunc.indices);

		bool allfixed=true;
		for (int i=0; i<(int)blockfunc.indices.size() && allfixed; ++i)
			if (lower.getElements()[i]!=upper.getElements()[i]) allfixed=false;
		if (allfixed) {
			clog << "CurvatureCheck: All variables in block of function " << *blockfunc.nonquad << " fixed. Consider it as constant." << endl;
			blockfunc.curvature=CONVEXCONCAVE;
			if (IsValid(blockfunc.quad)) {
				blockfunc.quad_mineig=0.;
				blockfunc.quad_maxeig=0.;
			}
			return;
		}
			
		Sampling sampling;
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
//		cout << "Checking curvature of " << *blockfunc.nonquad << " using " << blockfunc.samplepoints.size() << " sample points." << endl;
		SampleSet::iterator it_sp(blockfunc.samplepoints.begin());
		SymSparseMatrixCreator hessian_creator;			
		int eigval_successes=0; // number of successful eigenvalue computations
		while (it_sp!=blockfunc.samplepoints.end()) {
			try {
				blockfunc.nonquad->fullHessian(hessian_creator, *it_sp);
			} catch (FunctionEvaluationError error) {
				cerr << "CurvatureCheck::computeCurvature: skip sample point due to " << error << endl;
				it_sp=blockfunc.samplepoints.eraseAndGetNext(it_sp);
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
		if (blockfunc.samplepoints.size()<2) {
			cerr << "Not enough sample points for reliable computation of curvature type of block of function " << *blockfunc.nonquad << endl;
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
