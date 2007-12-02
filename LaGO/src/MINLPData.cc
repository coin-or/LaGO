// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "MINLPData.h"


MINLPData::ObjCon::ObjCon(const Pointer<SepQcFunc>& func_, const string& name_)
: name(name_), func(func_), functype(SepQcFunc::CONSTANT), curvtype(func_->get_curvature())
{	for (int k=0; k<func->block.size(); ++k) {
		if (func->b[k] && functype==SepQcFunc::CONSTANT) functype=SepQcFunc::LINEAR;
		if (func->A[k] && (functype==SepQcFunc::CONSTANT || functype==SepQcFunc::LINEAR)) functype=SepQcFunc::QUADRATIC;
		if (func->s[k]) { functype=SepQcFunc::NONQUAD; break; }
	}
}

ostream& operator<<(ostream& out, const MINLPData::ObjCon& objcon) {
//	out << "Function: " << *objcon.func;
//	if (objcon.convex_underestimator) out << "Underestimator: " << *objcon.convex_underestimator;
//	out << "Convexification characteristica: ";
//	for (map<int,double>::const_iterator it(objcon.convexification_characteristica_lower.begin()); it!=objcon.convexification_characteristica_lower.end(); ++it)
//		out << '(' << it->first << ", " << it->second << ") ";
//	out << endl;
	out << "Blocks mapped to following constraints by reformulation: ";
	for (map<int,int>::const_iterator it(objcon.reformulation_constraints_lower.begin()); it!=objcon.reformulation_constraints_lower.end(); ++it)
		out << it->first << "->" << it->second << '\t';
	out << endl;
	return out;
}

ostream& operator<<(ostream& out, const MINLPData::Constraint& con) {
	out << (MINLPData::ObjCon&)con;
//	if (con.concave_overestimator) out << "Overestimator: " << *con.concave_overestimator;
//	out << "Convexification characteristica: ";
//	for (map<int,double>::const_iterator it(con.convexification_characteristica_upper.begin()); it!=con.convexification_characteristica_upper.end(); ++it)
//		out << '(' << it->first << ", " << it->second << ") ";
//	out << endl;
	if (con.equality) {
//		out << "ineq_index: " << con.ineq_index << endl;
		out << "Blocks of upper part mapped to following constraints by reformulation: ";
		for (map<int,int>::const_iterator it(con.reformulation_constraints_upper.begin()); it!=con.reformulation_constraints_upper.end(); ++it)
			out << it->first << "->" << it->second << '\t';
		out << endl;
	}
	return out;
}


MINLPData::MINLPData()
{
}

MINLPData::MINLPData(const MinlpProblem& prob)
{
	var.reserve(prob.dim());
	con.reserve(prob.con.size());
	name=prob.prob_name;
	block.resize(prob.block.size());
	discrete_var=prob.i_discr;
	
	for (int i=0; i<prob.dim(); ++i)
		var.push_back(Variable(i, prob.lower(i), prob.upper(i), prob.discr[i], string(prob.var_names[i] ? (const char*)prob.var_names[i] : "")));
	for (int k=0; k<prob.block.size(); ++k) {
		block[k].reserve(prob.block[k].size());
		for (int i=0; i<prob.block[k].size(); ++i) {
			int i0=prob.block[k][i];
			block[k].push_back(i0);
			var[i0].block_nr=k;
			var[i0].index_in_block=i;
		}
	}
	
	for (int c=0; c<prob.con.size(); ++c) {
		con.push_back(Constraint(c, prob.con[c], prob.con_eq[c], string(prob.con_names[c] ? (const char*)prob.con_names[c] : "")));
	}
	
	obj=Objective(prob.obj);	
}

MINLPData::~MINLPData()
{
}
