// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOMINLPData.hpp"

namespace LaGO {

MINLPData::ObjCon::ObjCon(const SmartPtr<Function>& origfuncNL_, const SmartPtr<SparseVector>& origfuncLin_, double origfuncConstant_, const string& name_)
: name(name_), origfuncNL(origfuncNL_), origfuncLin(origfuncLin_), origfuncConstant(origfuncConstant_)
{ }

MINLPData::ObjCon::~ObjCon() { }

ostream& operator<<(ostream& out, const MINLPData::ObjCon& objcon) {
	out << objcon.name << endl;
	return out;
}

ostream& operator<<(ostream& out, const MINLPData::Variable& var) {
	out << var.index << ": " << var.name << " [" << var.lower << ", " << var.upper << "] ";
	if (var.discrete) out << "discrete";
	out << endl;
	return out;	
}

ostream& operator<<(ostream& out, const MINLPData::Constraint& con) {
	out << con.index << ": [" << con.lower << ", " << con.upper << "] ";
	out << (MINLPData::ObjCon&)con;
	return out;
}

MINLPData::MINLPData()
{ }

MINLPData::~MINLPData() { }

ostream& operator<<(ostream& out, const MINLPData& data) {
	out << "MINLP " << data.name << ':' << endl;
	out << data.var.size() << " variables:" << endl;
	for (unsigned int i=0; i<data.var.size(); ++i)
		out << data.var[i];
	out << "Objective: " << data.obj;
	out << data.con.size() << " constraints:" << endl;
	for (unsigned int c=0; c<data.con.size(); ++c)
		out << data.con[c];
	out << data.start_points.size() << " start points:" << endl;
	for (unsigned int i=0; i<data.start_points.size(); ++i)
		out << data.start_points[i] << endl;
	return out;
}

} // namespace LaGO
