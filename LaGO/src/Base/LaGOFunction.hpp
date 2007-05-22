// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOConfig.h 94 2007-05-21 13:54:40Z stefan $

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

//	virtual void fullHessian(const DenseVector& x) throw FunctionEvaluationError const=0;	
}; // class Function

void Function::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
	value=eval(x);
	gradient(grad,x);
};

	
} // namespace LaGO

#endif /*LAGOFUNCTION_HPP_*/
