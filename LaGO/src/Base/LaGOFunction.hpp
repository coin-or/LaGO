// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOFUNCTION_HPP_
#define LAGOFUNCTION_HPP_

#include "LaGObase.hpp"
#include "CoinError.hpp"

namespace LaGO {

class DenseVector;
class IntervalVector;
class SymSparseMatrixCreator;
	
class FunctionEvaluationError : public CoinError {
public:
	FunctionEvaluationError(const string& message, const string& className, const string& methodName)
	: CoinError(message, methodName, className)
	{ }

	friend ostream& operator<<(ostream& out, FunctionEvaluationError& error) {
		out << error.className() << "::" << error.methodName() << ": " << error.message();
		return out;
	}
}; // class FunctionEvaluationError 

/** Abstract base class for a function.
 */	
class Function : public ReferencedObject {
public:
	virtual double eval(const DenseVector& x) const=0;
	
	virtual void gradient(DenseVector& grad, const DenseVector& x) const=0;
	
	virtual void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
		value=eval(x);
		gradient(grad,x);
	}

	virtual void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const=0;
	
	virtual void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const;
	
#ifdef COIN_HAS_FILIB
	virtual bool canIntervalEvaluation() const { return false; }
	
	virtual interval<double> eval(const IntervalVector& x) const { throw CoinError("interval evaluation not possible", "eval(const IntervalVector&)", "Function"); }

	virtual void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const { throw CoinError("interval evaluation not possible", "evalAndGradient(interval<double>& value, IntervalVector&, const IntervalVector&)", "Function"); }
#endif

	/** Indicates whether the function knows about the variables that appear in it.
	 */
	virtual bool haveSparsity() const { return false; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */ 	
	virtual const vector<int>& getSparsity() const { throw CoinError("sparsity information not available", "getSparsity()", "Function"); }
	
	virtual void print(ostream& out) const=0;	

	friend ostream& operator<<(ostream& out, const Function& f) { f.print(out); return out; }
}; // class Function
	
} // namespace LaGO

#endif /*LAGOFUNCTION_HPP_*/
