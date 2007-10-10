// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOBlockFunction.hpp"

namespace LaGO {
	
double BlockFunction::eval(const DenseVector& x) const {
	assert(x.getNumElements()==(int)indices.size());
	double val=0;
	if (IsValid(quad)) val+=quad->xAx(x);
	if (IsValid(nonquad)) val+=nonquad->eval(x);
	return val;
}

void BlockFunction::gradient(DenseVector& grad, const DenseVector& x) const {
	assert(x.getNumElements()==(int)indices.size());
	assert(grad.getNumElements()==(int)indices.size());
	if (IsValid(nonquad)) {
		nonquad->gradient(grad, x);
		if (IsValid(quad))
			quad->addMultVector(grad, x, 2.); 
	} else if (IsValid(quad))
		quad->multVector(grad, x, 2.);
	else grad.clear();
}

void BlockFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
	assert(x.getNumElements()==(int)indices.size());
	assert(grad.getNumElements()==(int)indices.size());
	if (IsValid(nonquad)) {
		nonquad->evalAndGradient(value, grad, x);
		if (IsValid(quad)) {
			quad->addMultVector(grad, x, 2.);
			value+=quad->xAx(x);
		} 
	} else if (IsValid(quad)) {
		quad->multVector(grad, x, 2.);
		value=.5*(grad*x);
	} else {
		grad.clear();
		value=0.;
	}
}

void BlockFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
	assert(x.getNumElements()==(int)indices.size());
	assert(product.getNumElements()==(int)indices.size());
	assert(factor.getNumElements()==(int)indices.size());
	if (IsValid(nonquad)) {
		nonquad->hessianVectorProduct(product, x, factor);
		if (IsValid(quad)) quad->addMultVector(product, factor, 2.);		
	} else if (IsValid(quad)) quad->multVector(product, factor, 2.);
	else product.clear();	
}

void BlockFunction::fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const {
	if (IsValid(nonquad)) nonquad->fullHessian(hessian, x);
	else hessian.clear();
	if (IsValid(quad))
		hessian.add(2., *quad);
}

#ifdef COIN_HAS_FILIB
interval<double> BlockFunction::eval(const IntervalVector& x) const {
	interval<double> val;
	if (IsValid(nonquad))	val=nonquad->eval(x);
//	if (IsValid(quad)) val+=quad->xAx(x);
	return val;
}

void BlockFunction::evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
	assert(x.getNumElements()==(int)indices.size());
	assert(grad.getNumElements()==(int)indices.size());
	if (IsValid(nonquad)) {
		nonquad->evalAndGradient(value, grad, x);
		if (IsValid(quad)) {
//			quad->addMultVector(grad, x, 2.);
//			value+=quad->xAx(x);
		} 
	} else if (IsValid(quad)) {
//		quad->multVector(grad, x, 2.);
//		value=quad->xAx(x);
	} else {
		grad.clear();
		value=interval<double>(0.);
	}
}
#endif // COIN_HAS_FILIB

double BlockFunction::evalUnderEstimator(DenseVector& x) const {
	double val=-getInfinity();
	for (list<SmartPtr<QuadraticEstimator> >::const_iterator it_est(underestimators.begin()); it_est!=underestimators.end(); ++it_est) {
		double estval=(*it_est)->func->eval(x);
		if (estval>val) val=estval;
	}
	
	return val;
}

double BlockFunction::evalOverEstimator(DenseVector& x) const {
	double val=getInfinity();
	for (list<SmartPtr<QuadraticEstimator> >::const_iterator it_est(overestimators.begin()); it_est!=overestimators.end(); ++it_est) {
		double estval=(*it_est)->func->eval(x);
		if (estval<val) val=estval;
	}

	return val;
}

void BlockFunction::print(ostream& out) const {
	out << "Curvature: " << curvature << " Function: ";
	if (IsValid(quad)) out << *quad << endl;
	if (IsValid(nonquad)) out << *nonquad << endl;
	if (IsValid(sparsitygraph)) out << "SparsityGraph: " << *sparsitygraph;
}

} // namespace LaGO
