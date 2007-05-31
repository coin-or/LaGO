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

Curvature MINLPData::ObjCon::getCurvature() const {
	if (IsNull(origfuncNL)) return CONVEXCONCAVE; // linear constraint
	if (decompfuncNL.empty()) return UNKNOWN; // curvature type not determined yet
	Curvature curv=CONVEXCONCAVE;
	for (unsigned int k=0; k<decompfuncNL.size() && curv!=UNKNOWN; ++k)
		curv=addCurvatures(curv, decompfuncNL[k]->curvature);
	return curv;
}

void printSparseVector(ostream& out, const SparseVector& v, const vector<MINLPData::Variable>& var) {
	for (int i=0; i<v.getNumElements(); ++i) {
		double coeff=v.getElements()[i];
		if (coeff==1) out << '+';
		else if (coeff==-1) out << '-';
		else if (coeff>=0) out << '+' << coeff << '*';
		else out << coeff << '*'; 
		out << var[v.getIndices()[i]].getName();
	}
}

void MINLPData::ObjCon::print(ostream& out, const vector<MINLPData::Variable>& var) const {
	out << name << ": Curvature: " << getCurvature() << " Function: " << origfuncConstant;
	if (IsValid(origfuncLin))
		printSparseVector(out, *origfuncLin, var);
	if (IsValid(origfuncNL))
		out << '+' << *origfuncNL;
	out << endl;
		
	if (IsValid(sparsitygraph) && sparsitygraph->size())
		out << "Sparsitygraph: " << *sparsitygraph;
	 
	if (!decompfuncNL.empty()) {
		out << "Decomposed function: linear=" << decompfuncConstant;
		if (IsValid(decompfuncLin))
			printSparseVector(out, *decompfuncLin, var);
		out << endl;
		for (unsigned int k=0; k<decompfuncNL.size(); ++k)
			out << "Block " << k << ": " << *decompfuncNL[k];
		out << "Variables mapping:";
		for (unsigned int i=0; i<decompmapping.size(); ++i) {
			if (decompmapping[i].empty()) continue;
			out << ' ' << var[i].getName() << "->";
			for (unsigned int j=0; j<decompmapping[i].size(); ++j)
				out << '(' << decompmapping[i][j].first << ',' << decompmapping[i][j].second << ')';
		}			  
		out << endl;
	}	 
}

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

void MINLPData::Constraint::print(ostream& out, const vector<MINLPData::Variable>& var) const {
	out << index << ": [" << lower << ", " << upper << "] ";
	MINLPData::ObjCon::print(out, var);
}

ostream& operator<<(ostream& out, const MINLPData::Constraint& con) {
	out << con.index << ": [" << con.lower << ", " << con.upper << "] ";
	out << (MINLPData::ObjCon&)con;
	return out;
}

MINLPData::MINLPData()
{ }

MINLPData::~MINLPData() { }

void MINLPData::getBox(DenseVector& lower, DenseVector& upper, const vector<int>& indices) const {
	lower.resize(indices.size());
	upper.resize(indices.size());
	for (unsigned int i=0; i<indices.size(); ++i) {
		const Variable& var(getVariable(indices[i]));
		lower[i]=var.getLower();
		upper[i]=var.getUpper();
	}
}
void MINLPData::getBox(DenseVector& lower, DenseVector& upper) const {
	lower.resize(var.size());
	upper.resize(var.size());
	for (unsigned int i=0; i<var.size(); ++i) {
		lower[i]=var[i].getLower();
		upper[i]=var[i].getUpper();
	}
}


bool MINLPData::isConvex() const {
	if (!(obj.getCurvature() & CONVEX)) return false;
	for (unsigned int c=0; c<con.size(); ++c) {
		Curvature curv=con[c].getCurvature();
		// if there is a lower bound, then the constraint function need to be concave
		if (con[c].lower>-getInfinity() && !(curv & CONCAVE)) return false;
		// if there is an upper bound, then the constraint function need to be concave
		if (con[c].upper< getInfinity() && !(curv & CONVEX )) return false;
	}
	return true;
}

ostream& operator<<(ostream& out, const MINLPData& data) {
	out << "MINLP " << data.name << ": ";
	if (data.isConvex()) out << "is convex";
	out << endl;
	out << data.var.size() << " variables:" << endl;
	for (unsigned int i=0; i<data.var.size(); ++i)
		out << data.var[i];
	out << "Objective: ";
	data.obj.print(out, data.var);
	out << data.con.size() << " constraints:" << endl;
	for (unsigned int c=0; c<data.con.size(); ++c)
		data.con[c].print(out, data.var);
	out << data.start_points.size() << " start points:" << endl;
	for (unsigned int i=0; i<data.start_points.size(); ++i)
		out << data.start_points[i] << endl;
	return out;
}

} // namespace LaGO
