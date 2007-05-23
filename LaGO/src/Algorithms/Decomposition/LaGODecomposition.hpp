// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGODECOMPOSITION_HPP_
#define LAGODECOMPOSITION_HPP_

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"

namespace LaGO {

class Decomposition {
private:
	MINLPData& data;
public:
	Decomposition(MINLPData& data_)
	: data(data_)
	{ }
	
	void decompose();
	
	void decompose(MINLPData::ObjCon& objcon);
	
}; // class Decomposition
	
} // namespace LaGO

#endif /*LAGODECOMPOSITION_HPP_*/
