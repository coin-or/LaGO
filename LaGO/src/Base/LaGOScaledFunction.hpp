// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOSCALEDFUNCTION_HPP_
#define LAGOSCALEDFUNCTION_HPP_

namespace LaGO {

/** Abstract base class for a function multiplied by a constant
 */	
class ScaledFunction : public Function {
	SmartPtr<Function> func;
	double factor;
public:
	ScaledFunction(const SmartPtr<Function>& func_, double factor_)
	: func(func_), factor(factor_)
	{ }
	
	double eval(const DenseVector& x) const { return factor*func->eval(x); }
	
	void gradient(DenseVector& grad, const DenseVector& x) const {
		func->gradient(grad, x);
		grad*=factor;
	}
	
	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) {
		func->evalAndGradient(value, grad, x);
		value*=factor;
		grad*=factor;
	}

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
		func->hessianVectorProduct(product, x, factor);
		product*=this->factor;
	}
	
	/** Indicates whether the function knows about the variables that appear in it.
	 */
	bool haveSparsity() const { return func->haveSparsity(); }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */ 	
	const vector<int>& getSparsity() const { return func->getSparsity(); }
	
	void print(ostream& out) const {
		out << factor << '*';
		func->print(out);
	}
}; // class ScaledFunction

} // namespace LaGO

#endif /*LAGOSCALEDFUNCTION_HPP_*/
