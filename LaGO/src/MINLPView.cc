// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "MINLPView.h"

ostream& operator<<(ostream& out, const MINLPView& minlpview) {
	out << "MINLP: " << minlpview.dim() << " variables in " << minlpview.nr_blocks() << " blocks, " << minlpview.nr_discr() << " discrete variables, " << minlpview.connr() << " constraints." << endl;
	out << "Variables: " << endl;
	for (int i=0; i<minlpview.dim(); ++i)
		out << *minlpview.varPtr(i) << endl;
	out << "Objective: " << *minlpview.objPtr() << endl;
	out << "Constraints: " << endl;
	for (int c=0; c<minlpview.connr(); ++c)
		out << *minlpview.conPtr(c) << endl;	
	
	return out;	
}

ostream& operator<<(ostream& out, const MINLPView::VariableView& varview) {
	out << varview.index() << ' ' << '(' << varview.block_nr() << ',' << varview.index_in_block() << ")\t " << varview.name() << '\t';
	out << (varview.discrete() ? " is discrete" : " is continuous");
	out << "\t box=[" << varview.lower() << ',' << varview.upper() << ']';	
	
	return out;
}

ostream& operator<<(ostream& out, const MINLPView::ConstraintView& conview) {
	out << conview.index() << ' ' << conview.name() << " function type: " << conview.functype() << " curvature type: " << conview.curvature() << ' ';
	out << (conview.equality() ? "is equality" : "is inequality");
//	out << (conview.equality() ? "0 = " : "0 >= ") << *conview.func();
	return out;
}

ostream& operator<<(ostream& out, const MINLPView::ObjectiveView& objview) {
	out << objview.name() << " function type: " << objview.functype() << " curvature type: " << objview.curvature() << ' ';
//	out << *objview.func();
	return out;
}


void MINLPView::ObjConView::gradient(double* grad, const double* x) const {
	dvector grad_(dim());
	gradient(grad_, dvector(x, dim()));
	grad_ >> grad;
}

double MINLPView::ObjConView::gradientmult(const UserVector<double>& x, const UserVector<double>& factor) const {
	Pointer<UserVector<double> > grad(x.getemptycopy());
	gradient(*grad, x);
	return *grad*factor;
}

