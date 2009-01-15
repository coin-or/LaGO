// Copyright (C) Stefan Vigerske 2009
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOOSFUNCTION_HPP_
#define LAGOOSFUNCTION_HPP_

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

class OSInstance;
class OSExpressionTree;
class OSnLNode;

namespace LaGO {

class OSFunction : public Function {
private:
	OSInstance* osinstance;
	int osconidx;
	
	OSExpressionTree* exptree;

	std::vector<int> sparsity;
	
#ifdef COIN_HAS_COUENNE
	expression* createCouenneExpression(OSnLNode* node, std::vector<exprVar*>& vars, Domain* domain) const;
#endif
	
public:
	OSFunction(OSInstance* osinstance_, int osconidx_);
	
	~OSFunction();
	
	double eval(const DenseVector& x) const;

	void gradient(DenseVector& grad, const DenseVector& x) const;
	
	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const;
	
	void fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const;
	
	/** Indicates whether the function knows about the variables that appear in it.
	 */
	bool haveSparsity() const { return true; }

	/** Returns a list of variable indices that appear in this function.
	 * You can only rely on the result of this function if haveSparsity() returns true.
	 */ 	
	const std::vector<int>& getSparsity() const { return sparsity; }
	
#ifdef COIN_HAS_COUENNE
	expression* getAsCouenneExpression(std::vector<exprVar*>& vars, Domain* domain) const;
#endif
	
	void print(ostream& out) const;
};

} /* namespace LaGO */

#endif /*LAGOOSFUNCTION_HPP_*/
