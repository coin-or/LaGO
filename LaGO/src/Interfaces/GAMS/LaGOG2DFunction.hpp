// Copyright (C) Stefan Vigerske 2009
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOG2DFUNCTION_HPP_
#define LAGOG2DFUNCTION_HPP_

// filib uses MIN and MAX
#ifdef MIN
#undef MIN
#endif
#ifdef MAX
#undef MAX
#endif
#include "LaGObase.hpp"
#include "LaGOFunction.hpp"

#ifdef COIN_HAS_COUENNE
class expression;
class exprVar;
class Domain;
#endif

namespace LaGO {

class G2DFunction : public Function {
private:
	int n_instr;
	unsigned int* instr;
	double* constants;
	
	double* v;
	double* vbar;
	
	vector<int> sparsity;
	
public:
	double workfactor;
	
	G2DFunction(int n_instr_, unsigned int* instr_, double* constants_);
	
	~G2DFunction();
	
	double eval(const DenseVector& x) const;

	void gradient(DenseVector& grad, const DenseVector& x) const;
	
	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const;
	
	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const;
	
#ifdef COIN_HAS_FILIB
	bool canIntervalEvaluation() const { return true; }
	
	interval<double> eval(const IntervalVector& x) const;

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const;
#endif

	/** Indicates whether the function knows about the variables that appear in it.
	 */
	bool haveSparsity() const { return true; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */ 	
	const vector<int>& getSparsity() const { return sparsity; }
	
#ifdef COIN_HAS_COUENNE
	expression* getAsCouenneExpression(std::vector<exprVar*>& vars, Domain* domain) const;
#endif
	
	void print(ostream& out) const;
};

} /* namespace LaGO */

#endif /*LAGOG2DFUNCTION_HPP_*/
