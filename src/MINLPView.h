// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef _MINLPVIEW_H_
#define _MINLPVIEW_H_

#include "standard.h"
#include "MINLPData.h"

/** Defines a view on a MINLP.
 * Abstract base class.
 */
class MINLPView {
protected:
	/** Defines a view on a objective or constraint of a MINLP.
	 */
	class ObjConView {
	public:
		virtual ~ObjConView() { }
	
		virtual int dim() const=0;
		virtual const string& name() const=0;
		virtual SepQcFunc::ftype functype() const=0;
		virtual Func::CurvatureType curvature() const=0;
		
		/** Evaluates the constraint.
		 * TODO: Can throw an evaluation exception.
		 * @param x The point where to evaluate the constraint.
		 * @return The value of the constraint in x.
		 */
		virtual double evaluate(const UserVector<double>& x) const=0;
		
		/** Evaluates the constraint.
		 * TODO: Can throw an evaluation exception.
		 * @param x The point where to evaluate the constraint.
		 * @return The value of the constraint in x.
		 */
		virtual double evaluate(const double* x) const { return evaluate(dvector(x, dim())); }
		
		/** Computes the gradient of the constraint.
		 * TODO: Can throw an evaluation exception.
		 * @param grad The vector to store the gradient in.
		 * @param x The point where to compute the gradient.
		 */
		virtual void gradient(UserVector<double>& grad, const UserVector<double>& x) const=0;
		
		/** Computes the gradient of the constraint.
		 * TODO: Can throw an evaluation exception.
		 * @param grad The vector to store the gradient in.
		 * @param x The point where to compute the gradient.
		 */
		virtual void gradient(double* grad, const double* x) const;
		
		/** Computes the product of the gradient of the constraint with a vector.
		 * TODO: Can throw an evaluation exception.
		 * @param x The point where to compute the gradient.
		 * @param factor The vector to multiply with.
		 * @return The product between gradient evaluated in x and the factor.
		 */
		virtual double gradientmult(const UserVector<double>& x, const UserVector<double>& factor) const;

		/** Evaluates the constraint and computes is gradient.
		 * TODO: Can throw an evaluation exception.
		 * @param grad The vector where the gradient is stored.
		 * @param x The point where to evaluate the constraint and to compute the gradient.
		 * @return The value of the constraint.
		 */
		virtual double evaluate_and_gradient(UserVector<double>& grad, const UserVector<double>& x) const { gradient(grad, x); return evaluate(x); }
		
		/** Computes the product of the hessian with a vector.
		 * TODO: Can throw an evaluation exception.
		 * @param prod The vector to store the product in.
		 * @param x The point where to evaluate the constraint.
		 * @param factor The vector to multiply the hessian with.
		 */
		virtual void hessianmult(UserVector<double>& prod, const UserVector<double>& x, const UserVector<double>& factor) const=0;
	};
	
public:
	friend ostream& operator<<(ostream& out, const MINLPView& minlpview);
	/** Defines a view on a variable of a MINLP.
	 */
	class VariableView {
	public:
		friend ostream& operator<<(ostream& out, const VariableView& varview);
	
		virtual int index() const=0;
		virtual int block_nr() const=0;
		virtual int index_in_block() const=0;
		virtual const string& name() const=0;
		virtual double lower() const=0;
		virtual double upper() const=0;
		virtual pair<double, double> bounds() const { return pair<double,double>(lower(), upper()); }
		virtual bool discrete() const=0;
	};
	
	/** Defines a view on a objective of a MINLP.
	 */
	class ObjectiveView : public virtual ObjConView {
		friend ostream& operator<<(ostream& out, const ObjectiveView& objview);
	};
	
	/** Defines a view on a constraint of a MINLP.
	 */
	class ConstraintView : public virtual ObjConView {
	public:
		friend ostream& operator<<(ostream& out, const ConstraintView& conview);
		virtual int index() const=0;
		virtual bool equality() const=0;
	};
	
	/** Defines a view on a block (of variables) of a MINLP.
	 */
	class BlockView {
	public:
		virtual int size() const=0;
		
		virtual int operator()(int index) const=0;
	};
	
	//TODO: VariableIterator, ConstraintIterator, BlockIterator...


	virtual ~MINLPView() { }

	virtual int dim() const=0;
	
	virtual int connr() const=0;
	
	virtual int nr_discr() const=0;
	
	virtual int nr_blocks() const=0;
	
	virtual Pointer<BlockView> blockPtr(int block_nr) const=0;
	
	virtual Pointer<VariableView> varPtr(int index) const=0;
	
	virtual Pointer<ConstraintView> conPtr(int index) const=0;
	
	virtual Pointer<ObjectiveView> objPtr() const=0;

	// evaluation and so on...
};

#endif //_MINLPVIEW_H_
