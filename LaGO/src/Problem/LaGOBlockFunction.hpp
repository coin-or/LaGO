// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOBLOCKFUNCTION_HPP_
#define LAGOBLOCKFUNCTION_HPP_

#include "LaGObase.hpp"
#include "LaGOSymSparseMatrix.hpp"
#include "LaGOSparsity.hpp"
#include "LaGOCurvature.hpp"

namespace LaGO {

/** A block from a larger nonlinear function.
 */
class BlockFunction : public ReferencedObject {
public:
	/** The indices of the variables this function is defined for.
	 */
	vector<int> indices;
	
	/** The nonquadratic part.
	 */
	SmartPtr<Function> nonquad;
	
	/** The quadratic part.
	 */
	SmartPtr<SymSparseMatrix> quad;
	
	SmartPtr<SparsityGraph> sparsitygraph;

	Curvature curvature;
	
	list<DenseVector> samplepoints;

	BlockFunction()
	{ }
		 	
	BlockFunction(const vector<int>& indices_)
	: indices(indices_), curvature(UNKNOWN)
	{ }
	
	friend ostream& operator<<(ostream& out, const BlockFunction& block) {
		out << "Curvature: " << block.curvature << " Function: ";
		if (IsValid(block.quad)) out << *block.quad << endl;
		if (IsValid(block.nonquad)) out << *block.nonquad << endl;
		if (IsValid(block.sparsitygraph)) out << "SparsityGraph: " << *block.sparsitygraph;
		return out;
	} 
	
}; // class BlockFunction

} // namespace LaGO

#endif /*LAGOBLOCKFUNCTION_HPP_*/
