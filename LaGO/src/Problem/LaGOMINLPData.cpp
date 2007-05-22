// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOConfig.h 94 2007-05-21 13:54:40Z stefan $

#include "LaGOMINLPData.hpp"
#include "LaGOBlockFunction.hpp"

namespace LaGO {

MINLPData::ObjCon::ObjCon(const SmartPtr<Function>& origfuncNL_, const SmartPtr<SparseVector>& origfuncLin_, double origfuncConstant_, const string& name_)
: name(name_), origfuncNL(origfuncNL_), origfuncLin(origfuncLin_), origfuncConstant(origfuncConstant_)
{ }

MINLPData::ObjCon::~ObjCon() { }

ostream& operator<<(ostream& out, const MINLPData::ObjCon& objcon) {
	return out;
}

ostream& operator<<(ostream& out, const MINLPData::Constraint& con) {
	out << (MINLPData::ObjCon&)con;
	return out;
}

MINLPData::MINLPData()
{ }

MINLPData::~MINLPData() { }

ostream& operator<<(ostream& out, const MINLPData& data) {
	
	return out;
}

} // namespace LaGO
