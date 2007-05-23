// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOFUNCTION_HPP_
#define LAGOFUNCTION_HPP_

#include "LaGObase.hpp"
#include "CoinError.hpp"

namespace LaGO {
	
class FunctionEvaluationError : public CoinError {
public:
	FunctionEvaluationError(const string& message, const string& methodName, const string& className)
	: CoinError(message, methodName, className)
	{ }

	friend ostream& operator<<(ostream& out, FunctionEvaluationError& error);
}; // class FunctionEvaluationError 

/** Abstract base class for a function.
 */	
class Function : public ReferencedObject {
public:	
	virtual double eval(const DenseVector& x) const=0;
	
	virtual void gradient(DenseVector& grad, const DenseVector& x) const=0;
	
	inline virtual void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;

	virtual void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const=0;
	
	/** Indicates whether the function knows about the variables that appear in it.
	 */
	virtual bool haveSparsity() const { return false; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */ 	
	virtual const vector<int>& getSparsity() const { throw CoinError("sparsity information not available", "getSparsity()", "Function"); }	

//	virtual void fullHessian(const DenseVector& x) throw FunctionEvaluationError const=0;	
}; // class Function

void Function::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
	value=eval(x);
	gradient(grad,x);
};


/** Abstract base class for a function.
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

//	virtual void fullHessian(const DenseVector& x) throw FunctionEvaluationError const=0;	
}; // class Function

	
} // namespace LaGO

#endif /*LAGOFUNCTION_HPP_*/
