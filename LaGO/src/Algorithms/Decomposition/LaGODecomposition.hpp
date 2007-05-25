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
	
	void findConnectedComponents(SparsityGraph& graph);
	bool setComponent(const SparsityGraph::Node& node, int comp);
public:
	Decomposition(MINLPData& data_)
	: data(data_)
	{ }
	
	void decompose();
	
	void decompose(MINLPData::ObjCon& objcon);
	
	void computeSparsityGraph(MINLPData::ObjCon& objcon, list<DenseVector>& samplepoints, const vector<int>& nonzeros);
	
	
}; // class Decomposition
	
} // namespace LaGO

#endif /*LAGODECOMPOSITION_HPP_*/
