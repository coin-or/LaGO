// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOSPARSITY_HPP_
#define LAGOSPARSITY_HPP_

#include "LaGObase.hpp"
// from Cgc (Coin Graph Classes)
#include "DynNet.h"

namespace LaGO {
	
class SparsityGraphNode {
public:
	int varindex;
	bool isquad;
	int component;

	SparsityGraphNode()
	: varindex(-1), isquad(true), component(-1)
	{ }

	SparsityGraphNode(int varindex_)
	: varindex(varindex_), isquad(true), component(-1)
	{ }

	friend ostream& operator<<(ostream& out, const SparsityGraphNode& node) {
		out << "var " << node.varindex << " in comp. " << node.component << (node.isquad ? " is quad." : " is nonquad.");
		return out;
	}
}; // SparsityGraphNode



class SparsityGraphEdge {
public:
	double hessian_entry;
	
	SparsityGraphEdge()
	: hessian_entry(0.)
	{ }

	SparsityGraphEdge(double hessian_entry_)
	: hessian_entry(hessian_entry_)
	{ }

	friend ostream& operator<<(ostream& out, const SparsityGraphEdge& edge) {
		out << edge.hessian_entry;
		return out;
	}
}; // SparsityGraphEdge

class SparsityGraph : public Cgc::DynNet<SparsityGraphNode, SparsityGraphEdge>, public ReferencedObject {
public:
	SparsityGraph(const int numNodes,const int numArcs)
	: Cgc::DynNet<SparsityGraphNode, SparsityGraphEdge>(numNodes, numArcs)
	{ }
};

} // namespace LaGO

#endif /*LAGOSPARSITY_HPP_*/
