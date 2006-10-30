// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef _MINLP_H_
#define _MINLP_H_

#include "standard.h"
#include "MINLPView.h"
#include "MINLPData.h"

/** Default view on a MINLP.
 */
class MINLP : public MINLPView {
private:
	template <class ObjOrCon> class ObjConView : public virtual MINLPView::ObjConView {
	protected:
		const ObjOrCon& objcon;
	public:
		ObjConView(const ObjOrCon& objcon_)
		: objcon(objcon_)
		{ }
		
		int dim() const { return objcon.func->dim(); }
		const string& name() const { return objcon.name; }
		SepQcFunc::ftype functype() const { return objcon.functype; }
		Func::CurvatureType curvature() const { return objcon.curvtype; }
		
		double evaluate(const UserVector<double>& x) const { return objcon.func->eval(x); }
		
		using MINLPView::ObjConView::evaluate;
		
		void gradient(UserVector<double>& grad, const UserVector<double>& x) const { objcon.func->grad(grad, x); }

		using MINLPView::ObjConView::gradient;
		
		double evaluate_and_gradient(UserVector<double>& grad, const UserVector<double>& x) const { double val; objcon.func->valgrad(val, grad, x); return val; }
		
		void hessianmult(UserVector<double>& prod, const UserVector<double>& x, const UserVector<double>& factor) const { objcon.func->HessMult(prod, x, factor); }
	};
	
public:
	class ConstraintView : public MINLPView::ConstraintView, public MINLP::ObjConView<MINLPData::Constraint> {
	public:
		ConstraintView(const MINLPData::Constraint& con_)
		: MINLP::ObjConView<MINLPData::Constraint>(con_)
		{ }

		int dim() const { return MINLP::ObjConView<MINLPData::Constraint>::dim(); }
		int index() const { return objcon.index; }
		bool equality() const { return objcon.equality; }
		const string& name() const { return MINLP::ObjConView<MINLPData::Constraint>::name(); }
		SepQcFunc::ftype functype() const { return MINLP::ObjConView<MINLPData::Constraint>::functype(); }
		Func::CurvatureType curvature() const { return MINLP::ObjConView<MINLPData::Constraint>::curvature(); }
		
		double evaluate(const UserVector<double>& x) const { return MINLP::ObjConView<MINLPData::Constraint>::evaluate(x); }
		using MINLPView::ConstraintView::evaluate;
		
		void gradient(UserVector<double>& grad, const UserVector<double>& x) const { MINLP::ObjConView<MINLPData::Constraint>::gradient(grad, x); }
		using MINLPView::ConstraintView::gradient;

		double evaluate_and_gradient(UserVector<double>& grad, const UserVector<double>& x) const { return MINLP::ObjConView<MINLPData::Constraint>::evaluate_and_gradient(grad, x); }
		
		void hessianmult(UserVector<double>& prod, const UserVector<double>& x, const UserVector<double>& factor) const { MINLP::ObjConView<MINLPData::Constraint>::hessianmult(prod, x, factor); }
	};
	
	class ObjectiveView : public MINLPView::ObjectiveView, public MINLP::ObjConView<MINLPData::Objective> {
	public:
		ObjectiveView(const MINLPData::Objective& obj_)
		: MINLP::ObjConView<MINLPData::Objective>(obj_)
		{ }

		int dim() const { return MINLP::ObjConView<MINLPData::Objective>::dim(); }
		const string& name() const { return MINLP::ObjConView<MINLPData::Objective>::name(); }
		SepQcFunc::ftype functype() const { return MINLP::ObjConView<MINLPData::Objective>::functype(); }
		Func::CurvatureType curvature() const { return MINLP::ObjConView<MINLPData::Objective>::curvature(); }
		
		double evaluate(const UserVector<double>& x) const { return MINLP::ObjConView<MINLPData::Objective>::evaluate(x); }
		using MINLPView::ObjectiveView::evaluate;
		
		void gradient(UserVector<double>& grad, const UserVector<double>& x) const { MINLP::ObjConView<MINLPData::Objective>::gradient(grad, x); }
		using MINLPView::ObjectiveView::gradient;

		double evaluate_and_gradient(UserVector<double>& grad, const UserVector<double>& x) const { return MINLP::ObjConView<MINLPData::Objective>::evaluate_and_gradient(grad, x); }
		
		void hessianmult(UserVector<double>& prod, const UserVector<double>& x, const UserVector<double>& factor) const { MINLP::ObjConView<MINLPData::Objective>::hessianmult(prod, x, factor); }
	};		

	class VariableView : public MINLPView::VariableView {
	private:
		const MINLPData::Variable& var;
	public:
		VariableView(const MINLPData::Variable& var_)
		: var(var_)
		{ }
		
		int index() const { return var.index; }
		int block_nr() const { return var.block_nr; }
		int index_in_block() const { return var.index_in_block; }
		const string& name() const { return var.name; }
		double lower() const { return var.lower; }
		double upper() const { return var.upper; }
		bool discrete() const { return var.discrete; }
	};

	
	class BlockView : public MINLPView::BlockView {
	private:
		const MINLPData::Block& block;
	public:
		BlockView(const MINLPData::Block& block_)
		: block(block_)
		{ }
		
		int size() const { return block.size(); }
		
		int operator()(int index) const { return block[index]; }
	};
	

protected:
	const MINLPData& data;
public:
	MINLP(const MINLPData& data_);
	
	~MINLP() { };

	int dim() const { return data.var.size(); }
	
	int connr() const { return data.con.size(); }
	
	int nr_discr() const { return data.discrete_var.size(); }
	
	int nr_blocks() const { return data.block.size(); }
	
	VariableView var(int index) const { return VariableView(data.var[index]); }
	
	ObjectiveView obj() const { return ObjectiveView(data.obj); }
	
	ConstraintView con(int index) const { return ConstraintView(data.con[index]); }
	
	BlockView block(int block_nr) const { return BlockView(data.block[block_nr]); }
	
	Pointer<MINLPView::BlockView> blockPtr(int block_nr) const { return new BlockView(data.block[block_nr]); }
	
	Pointer<MINLPView::VariableView> varPtr(int index) const { return new VariableView(data.var[index]); }
	
	Pointer<MINLPView::ConstraintView> conPtr(int index) const { return new ConstraintView(data.con[index]); }
	
	Pointer<MINLPView::ObjectiveView> objPtr() const { return new ObjectiveView(data.obj); }
	
};

#endif //_MINLP_H_
