// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOQuadraticFunction.hpp 109 2007-06-03 19:31:14Z stefan $

#include "LaGOConvexification.hpp"

#include <IpLapack.hpp>

namespace LaGO {

Convexification::Convexification(MINLPData& data_)
: data(data_)
{ }

void Convexification::convexify() {
	convexify(data.obj, true, false);
	for (int c=0; c<data.numConstraints(); ++c)
		convexify(data.con[c], data.con[c].upper<getInfinity(), data.con[c].lower>-getInfinity());
}
	
void Convexification::convexify(MINLPData::ObjCon& objcon, bool need_lower, bool need_upper) {
	for (unsigned int k=0; k<objcon.decompfuncNL.size(); ++k) {
		BlockFunction& func(*objcon.decompfuncNL[k]);

		bool do_lower=need_lower && !(func.curvature&CONVEX);
		bool do_upper=need_upper && !(func.curvature&CONCAVE);
		clog << objcon.name << " block " << k << " is " << func.curvature << '\t' << do_lower << do_upper << endl;
		
		if (do_lower || do_upper) {
			clog << "Convexify " << objcon.name << " block " << k << ": ";
			convexify(func, do_lower, do_upper);
			clog << endl;
		}
	}	
}

void Convexification::convexify(BlockFunction& func, bool do_lower, bool do_upper) {
	DenseVector diam;
	data.getBoxDiameter(diam, func.indices);

	//TODO: we could also convexify quad. underestimators together with quadratic terms
	//however, with the GAMS interface we do not have such constellation anyway
	if (IsValid(func.quad)) {
		convexify(*func.quad, diam,
			(do_lower && func.quad_mineig<0) ? &func.alpha_convexify : NULL,
			(do_upper && func.quad_maxeig>0) ? &func.alpha_concavify : NULL
		);
	}
	
	if (do_lower && IsValid(func.nonquad)) {
		assert(!func.underestimators.empty());
		for (list<SmartPtr<QuadraticEstimator> >::iterator it(func.underestimators.begin()); it!=func.underestimators.end(); ++it)
			convexify(*(*it)->func->A, diam, &(*it)->alpha, NULL); 
	}
	
	if (do_upper && IsValid(func.nonquad)) {
		assert(!func.overestimators.empty());
		for (list<SmartPtr<QuadraticEstimator> >::iterator it(func.overestimators.begin()); it!=func.overestimators.end(); ++it)
			convexify(*(*it)->func->A, diam, NULL, &(*it)->alpha); 
	}
	
}

void Convexification::convexify(SymSparseMatrix& matrix, const DenseVector& diam, DenseVector* alpha_convexify, DenseVector* alpha_concavify) {
  // lower triangular of quadratic term matrix, scaled by box diameter
  double* matrix_ = new double [matrix.getNumCols() * matrix.getNumCols()];
  CoinZeroN(matrix_, matrix.getNumCols() * matrix.getNumCols());

  for (int i=0; i<matrix.getNumNonzeros(); ++i) {
    int row = matrix.getRowIndices()[i];
    int col = matrix.getColIndices()[i];
    // compute value of matrix entry = q_ij * (u_i-l_i) * (u_j-l_j)
    if (row < col)
    	matrix_[row * matrix.getNumCols() + col] = matrix.getValues()[i] * diam[row] * diam[col];
    else
    	matrix_[col * matrix.getNumCols() + row] = matrix.getValues()[i] * diam[row] * diam[col];
  }

	double* eigval=new double[matrix.getNumCols()];
	int info;
  Ipopt::IpLapackDsyev(false, matrix.getNumCols(), matrix_, matrix.getNumCols(), eigval, info);
  
	if (info!=0) {
		cerr << "Convexification::convexify: Error computing eigenvalue of " << matrix << endl;
		exit(EXIT_FAILURE); 
	}

	if (alpha_convexify && eigval[0]<0) { // want convexification and matrix is not positive semidefinite
		clog << eigval[0] << ' ';
		alpha_convexify->resize(matrix.getNumCols());
		for (int i=0; i<matrix.getNumCols(); ++i)
			if (diam[i]==0.) (*alpha_convexify)[i]=0.;
			else (*alpha_convexify)[i]=eigval[0]/(diam[i]*diam[i]);
  }
	
	if (alpha_concavify && eigval[matrix.getNumCols()-1]>0) { // want concavification and matrix is not negative semidefinite
		clog << eigval[matrix.getNumCols()-1] << ' ';
		alpha_concavify->resize(matrix.getNumCols());
		for (int i=0; i<matrix.getNumCols(); ++i)
			if (diam[i]==0.) (*alpha_concavify)[i]=0.;
			else (*alpha_concavify)[i]=eigval[matrix.getNumCols()-1]/(diam[i]*diam[i]);
	}

  delete[] matrix_;
  delete[] eigval;
}	
	
} // namespace LaGO
