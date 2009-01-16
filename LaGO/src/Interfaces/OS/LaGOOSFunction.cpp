// Copyright (C) Stefan Vigerske 2009
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOConfig.h"
// have to be above LaGObase because LaGO gets the std namespace, giving a conflict with Couenne's unary_function
#ifdef COIN_HAS_COUENNE
#include "exprAbs.hpp"
#include "exprConst.hpp"
#include "exprCos.hpp"
#include "exprDiv.hpp"
#include "exprExp.hpp"
#include "exprInv.hpp"
#include "exprLog.hpp"
#include "exprMax.hpp"
#include "exprMin.hpp"
#include "exprMul.hpp"
#include "exprOpp.hpp"
#include "exprPow.hpp"
#include "exprSin.hpp"
#include "exprSub.hpp"
#include "exprSum.hpp"
#include "exprVar.hpp"
#endif

#include "LaGOOSFunction.hpp"
#include "LaGOSymSparseMatrix.hpp"

#include <sstream>

#include "OSInstance.h"
#include "OSErrorClass.h"

namespace LaGO {

OSFunction::OSFunction(OSInstance* osinstance_, int osconidx_)
: osinstance(osinstance_), osconidx(osconidx_), exptree(osinstance->getNonlinearExpressionTreeMod(osconidx))
{
	assert(exptree);
	map<int, int>* varmap = exptree->getVariableIndiciesMap();
	
	sparsity.reserve(varmap->size());
	for (map<int,int>::iterator it(varmap->begin()); it != varmap->end(); ++it)
		sparsity.push_back(it->first);
}

OSFunction::~OSFunction()
{
}

double OSFunction::eval(const DenseVector& x) const
{
	try
	{
		return exptree->calculateFunction(const_cast<double*>(x.getElements()), true);
	}
	catch (const ErrorClass& error)
	{
		throw FunctionEvaluationError(error.errormsg, "OSFunction", "eval");
	}
	
//	return 0.;
}

void OSFunction::gradient(DenseVector& grad, const DenseVector& x) const
{
//	cout << "x: " << x << endl;
	grad.clear();
	if (osconidx >= 0)
	{ // be constraint
		::SparseVector* vec;
		try
		{
			vec = osinstance->calculateConstraintFunctionGradient(const_cast<double*>(x.getElements()), osconidx, true);
		}
		catch (const ErrorClass& error)
		{
			throw FunctionEvaluationError(error.errormsg, "OSFunction", "gradient");
		}

		assert(vec);
		for (int i = 0; i < vec->number; ++i) {
			if (!(vec->values[i] == vec->values[i]))
				throw FunctionEvaluationError("nan in gradient", "OSFunction", "gradient");
			grad[vec->indexes[i]] = vec->values[i];
		}

		// OS puts the coefficients from the linear variables into the gradient
		// we don't want this, so we delete them again
		SparseMatrix* lincoeff = osinstance->getLinearConstraintCoefficientsInRowMajor();
		for (int j = lincoeff->starts[osconidx]; j < lincoeff->starts[osconidx+1]; ++j)
			grad[lincoeff->indexes[j]] -= lincoeff->values[j];
	}
	else
	{ // be objective
		double* osgrad;
		try
		{
			osgrad = osinstance->calculateObjectiveFunctionGradient(const_cast<double*>(x.getElements()), osconidx, true);
		}
		catch (const ErrorClass& error)
		{
			throw FunctionEvaluationError(error.errormsg, "OSFunction", "gradient");
		}
		CoinMemcpyN(osgrad, osinstance->getVariableNumber(), grad.getElements());

//		cout << "grad before: " << grad << endl;
		
		// OS puts the coefficients from the linear variables into the gradient
		// we don't want this, so we delete them again
		::SparseVector* lincoeff = osinstance->getObjectiveCoefficients()[0];
		for (int j = 0; j < lincoeff->number; ++j)
			grad[lincoeff->indexes[j]] -= lincoeff->values[j];		
	}
//	cout << "grad: " << grad << endl;

#ifndef NDEBUG
	map<int, int>* varmap = exptree->getVariableIndiciesMap();
	for (int i = 0; i < grad.getNumElements(); ++i)
		assert(varmap->count(i) || fabs(grad[i]) < 1e-12); // variables not in expression tree should have no gradient
#endif	
}

void OSFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const
{
	gradient(grad, x);
	value = exptree->calculateFunction(const_cast<double*>(x.getElements()), false);
}

void OSFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const
{
	SparseHessianMatrix* hessian;
	try
	{
		hessian = osinstance->calculateHessian(const_cast<double*>(x.getElements()), osconidx, true);
		assert(hessian);
	}
	catch (const ErrorClass& error)
	{
		throw FunctionEvaluationError(error.errormsg, "OSFunction", "hessianVectorProduct");
	}
	
	product.clear();
	for (int i = 0; i < hessian->hessDimension; ++i)
	{
		if (!(hessian->hessValues[i] == hessian->hessValues[i]))
			throw FunctionEvaluationError("nan in hessian", "OSFunction", "gradient");
		product[hessian->hessRowIdx[i]] += hessian->hessValues[i] * factor[hessian->hessColIdx[i]];
	}
}

void OSFunction::fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const
{
	SparseHessianMatrix* oshessian;
	try
	{
		oshessian = osinstance->calculateHessian(const_cast<double*>(x.getElements()), osconidx, true);
		assert(oshessian);
	}
	catch (const ErrorClass& error)
	{
		throw FunctionEvaluationError(error.errormsg, "OSFunction", "hessianVectorProduct");
	}

	for (int i = 0; i < oshessian->hessDimension; ++i)
	{
		if (!(oshessian->hessValues[i] == oshessian->hessValues[i]))
			throw FunctionEvaluationError("nan in hessian", "OSFunction", "gradient");
		hessian.insert(oshessian->hessColIdx[i], oshessian->hessRowIdx[i], oshessian->hessValues[i]);
	}
}


#ifdef COIN_HAS_COUENNE
expression* OSFunction::getAsCouenneExpression(std::vector<exprVar*>& vars, Domain* domain) const
{
	return createCouenneExpression(exptree->m_treeRoot, vars, domain);
}

expression* OSFunction::createCouenneExpression(OSnLNode* node, std::vector<exprVar*>& vars, Domain* domain) const {
  switch (node->inodeInt) {
     case OS_PLUS :
    	 return new exprSum(createCouenneExpression(node->m_mChildren[0], vars, domain), createCouenneExpression(node->m_mChildren[1], vars, domain));
     case OS_SUM :
    	 switch (node->inumberOfChildren==0) {
    		 case 0:
    			 return new exprConst(0.);
    		 case 1:
    			 return createCouenneExpression(node->m_mChildren[0], vars, domain);
    		 default:
    			 expression** args = new expression*[node->inumberOfChildren];
           for(int i = 1; i < node->inumberOfChildren; ++i)
          	 args[i] = createCouenneExpression(node->m_mChildren[i], vars, domain);
           expression* base = new exprSum(args, node->inumberOfChildren);
           delete[] args;
           return base;
        }
     case OS_MINUS :
    	 return new exprSub(createCouenneExpression(node->m_mChildren[0], vars, domain), createCouenneExpression(node->m_mChildren[1], vars, domain));
     case OS_NEGATE :
    	 return new exprOpp(createCouenneExpression(node->m_mChildren[0], vars, domain));
     case OS_TIMES :
    	 return new exprMul(createCouenneExpression(node->m_mChildren[0], vars, domain), createCouenneExpression(node->m_mChildren[1], vars, domain));
     case OS_DIVIDE :
    	 // couenne does not like expressions of the form exp1/exp2 with exp1 a constant, so we write them as exp1 * 1/exp2
    	 if (node->m_mChildren[1]->inodeInt == OS_NUMBER)
      	 return new exprMul(createCouenneExpression(node->m_mChildren[0], vars, domain), new exprInv(createCouenneExpression(node->m_mChildren[1], vars, domain)));
    	 else
    		 return new exprDiv(createCouenneExpression(node->m_mChildren[0], vars, domain), createCouenneExpression(node->m_mChildren[1], vars, domain));
     case OS_POWER :
    	 // couenne does not like expressions of the form exp1 ^ exp2 with exp2 not a constant, so we write them as exp(log(exp1)*exp2)
    	 if (node->m_mChildren[1]->inodeInt != OS_NUMBER)
    		 return new exprExp(new exprMul(new exprLog(createCouenneExpression(node->m_mChildren[0], vars, domain)), createCouenneExpression(node->m_mChildren[1], vars, domain)));
    	 else
    		 return new exprPow(createCouenneExpression(node->m_mChildren[0], vars, domain), createCouenneExpression(node->m_mChildren[1], vars, domain));
     case OS_PRODUCT:
    	 switch (node->inumberOfChildren==0) {
    		 case 0:
    			 return new exprConst(1.);
    		 case 1:
    			 return createCouenneExpression(node->m_mChildren[0], vars, domain);
    		 default:
    			 expression** args = new expression*[node->inumberOfChildren];
           for(int i = 1; i < node->inumberOfChildren; ++i)
          	 args[i] = createCouenneExpression(node->m_mChildren[i], vars, domain);
           expression* base = new exprMul(args, node->inumberOfChildren);
           delete[] args;
           return base;
        }
     case OS_ABS :
    	 return new exprAbs(createCouenneExpression(node->m_mChildren[0], vars, domain));
     case OS_SQUARE :
    	 return new exprPow(createCouenneExpression(node->m_mChildren[0], vars, domain), new exprConst(2.));
     case OS_SQRT :
    	 return new exprPow(createCouenneExpression(node->m_mChildren[0], vars, domain), new exprConst(0.5));
     case OS_LN :
    	 return new exprLog(createCouenneExpression(node->m_mChildren[0], vars, domain));
     case OS_EXP :
    	 return new exprExp(createCouenneExpression(node->m_mChildren[0], vars, domain));
     case OS_SIN :
    	 return new exprSin(createCouenneExpression(node->m_mChildren[0], vars, domain));
     case OS_COS :
    	 return new exprCos(createCouenneExpression(node->m_mChildren[0], vars, domain));
     case OS_MIN :
    	 switch (node->inumberOfChildren==0) {
    		 case 0:
    			 return new exprConst(0.);
    		 case 1:
    			 return createCouenneExpression(node->m_mChildren[0], vars, domain);
    		 default:
    			 expression** args = new expression*[node->inumberOfChildren];
           for(int i=1;i<node->inumberOfChildren;i++)
          	 args[i] = createCouenneExpression(node->m_mChildren[i], vars, domain);
           expression* base = new exprMin(args, node->inumberOfChildren);
           delete[] args;
           return base;
        }
     case OS_MAX :
    	 switch (node->inumberOfChildren==0) {
    		 case 0:
    			 return new exprConst(0.);
    		 case 1:
    			 return createCouenneExpression(node->m_mChildren[0], vars, domain);
    		 default:
    			 expression** args = new expression*[node->inumberOfChildren];
           for(int i=1;i<node->inumberOfChildren;i++)
          	 args[i] = createCouenneExpression(node->m_mChildren[i], vars, domain);
           expression* base = new exprMax(args, node->inumberOfChildren);
           delete[] args;
           return base;
        }
     case OS_NUMBER :
    	 return new exprConst(((OSnLNodeNumber*)node)->value);
     case OS_PI :
    	 assert(false);
    	 //TODO
//    	 return new exprConst(PI);
     case OS_VARIABLE : {
        OSnLNodeVariable* varnode = (OSnLNodeVariable*)node;
        if (varnode->coef == 0.)
       		return new exprConst(0.);
        if (varnode->coef == 1.)
       		return vars[varnode->idx]->clone(domain);
        if (varnode->coef == -1.)
       		return new exprOpp(vars[varnode->idx]->clone(domain));
     		return new exprMul(new exprConst(varnode->coef), vars[varnode->idx]->clone(domain));
     }
     default:
        cout << node->snodeName << " NOT IMPLEMENTED!!" << endl;
        break;
  }
	
	return NULL;
}
#endif

void OSFunction::print(ostream& out) const
{
	out << "OSFunction on variables";
	for (size_t i = 0; i < sparsity.size(); ++i)
		out << ' ' << sparsity[i];
	out << endl;
}

} /* namespace minlp */
