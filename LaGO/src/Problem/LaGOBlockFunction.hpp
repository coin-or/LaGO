// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOBLOCKFUNCTION_HPP_
#define LAGOBLOCKFUNCTION_HPP_

#include "LaGObase.hpp"
#include "LaGOSymSparseMatrix.hpp"

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
	
	BlockFunction(const vector<int>& indices_)
	: indices(indices_)
	{ }
	
}; // class BlockFunction

} // namespace LaGO

#endif /*LAGOBLOCKFUNCTION_HPP_*/
