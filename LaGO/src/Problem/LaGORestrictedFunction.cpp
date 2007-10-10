// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGORestrictedFunction.hpp"

namespace LaGO {

void RestrictedFunction::gradient(DenseVector& grad, const DenseVector& x) const {
	fullx.setElementsOfBlock(x, indices);
	DenseVector fullgrad(fullx.getNumElements());
//	clog << "RestrictedFunction::gradient at " << x << endl;

	f->gradient(fullgrad, fullx);

	grad.setToBlock(fullgrad, indices);	
}
	
void RestrictedFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
	fullx.setElementsOfBlock(x, indices);
	DenseVector fullgrad(fullx.getNumElements());

	f->evalAndGradient(value, fullgrad, fullx);
	
	grad.setToBlock(fullgrad, indices);	
}

void RestrictedFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
	fullx.setElementsOfBlock(x, indices);
	DenseVector fullfactor(fullx.getNumElements());
	fullfactor.setElementsOfBlock(factor, indices);
	DenseVector fullproduct(fullx.getNumElements());

	f->hessianVectorProduct(fullproduct, fullx, fullfactor);
	
	product.setToBlock(fullproduct, indices);
}

#ifdef COIN_HAS_FILIB
interval<double> RestrictedFunction::eval(const IntervalVector& x) const {
	IntervalVector fullintx(fullx, fullx);
	fullintx.setElementsOfBlock(x, indices);
	
	return f->eval(fullintx);
}

void RestrictedFunction::evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
	IntervalVector fullintx(fullx, fullx);
	fullintx.setElementsOfBlock(x, indices);
	IntervalVector fullgrad(fullx.getNumElements());

	f->evalAndGradient(value, fullgrad, fullintx);
	
	grad.setToBlock(fullgrad, indices); 
}
#endif
	
} // namespace LaGO
