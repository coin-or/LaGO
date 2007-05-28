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
	
	/** Finds connected components in a sparsity graph and checks whether they belong to quadratic or nonquadratic parts of the function.
	 */
	void findConnectedComponents(vector<bool>& component_isnonquad, SparsityGraph& graph);
	/** Does a depth-first-search to find all nodes in a component. 
	 * @return True if the component found is nonquadratic, false else.
	 */
	bool setComponent(const SparsityGraph::Node& node, int comp);

	void createDecomposedFunctions(MINLPData::ObjCon& objcon, DenseVector& refpoint, const vector<int>& nonzeros, const list<int>& lin_nonzeros, const vector<bool>& component_isnonquad, bool have_quadratic_component);	
public:
	Decomposition(MINLPData& data_)
	: data(data_)
	{ }
	
	void decompose();
	
	void decompose(MINLPData::ObjCon& objcon);
	
	void computeSparsityGraph(MINLPData::ObjCon& objcon, list<int>& lin_nonzeros, list<DenseVector>& samplepoints, const vector<int>& nonzeros);
}; // class Decomposition
	
} // namespace LaGO

#endif /*LAGODECOMPOSITION_HPP_*/
