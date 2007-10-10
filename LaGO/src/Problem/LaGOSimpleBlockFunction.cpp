// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSimpleBlockFunction.hpp"

namespace LaGO {

double SimpleBlockFunction::eval(const DenseVector& x) const {
	my_x.setToBlock(x, indices);
//	clog << "SimpleBlockFunction of " << *f << " value at " << my_x << " is " << f->eval(my_x) << endl;
	return f->eval(my_x);
}

void SimpleBlockFunction::gradient(DenseVector& grad, const DenseVector& x) const {
	my_x.setToBlock(x, indices);
	DenseVector my_grad(my_x.getNumElements());

	f->gradient(my_grad, my_x);

//	clog << "SimpleBlockFunction of " << *f << " gradient at " << my_x;
//	clog << "\t is " << my_grad << endl;
	
	grad.clear();
	grad.setElementsOfBlock(my_grad, indices);	
}
	
void SimpleBlockFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
	my_x.setToBlock(x, indices);
	DenseVector my_grad(my_x.getNumElements());

	f->evalAndGradient(value, my_grad, my_x);
	
	grad.clear();
	grad.setElementsOfBlock(my_grad, indices);	
}

void SimpleBlockFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
	my_x.setToBlock(x, indices);
	DenseVector my_factor(my_x.getNumElements());
	my_factor.setToBlock(factor, indices);
	DenseVector my_product(my_x.getNumElements());

	f->hessianVectorProduct(my_product, my_x, my_factor);
	
	product.clear();
	product.setElementsOfBlock(my_product, indices);
}

#ifdef COIN_HAS_FILIB
interval<double> SimpleBlockFunction::eval(const IntervalVector& x) const {
	IntervalVector my_intx(my_x.getNumElements());
	my_intx.setToBlock(x, indices);
	
	return f->eval(my_intx);
}

void SimpleBlockFunction::evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
	IntervalVector my_intx(my_x.getNumElements());
	my_intx.setToBlock(x, indices);
	IntervalVector my_grad(my_x.getNumElements());

	f->evalAndGradient(value, my_grad, my_intx);
	
	grad.clear();
	grad.setElementsOfBlock(my_grad, indices); 
}
#endif

	
} // namespace LaGO
