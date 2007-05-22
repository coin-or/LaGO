// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOConfig.h 94 2007-05-21 13:54:40Z stefan $

#ifndef LAGOBLOCKFUNCTION_HPP_
#define LAGOBLOCKFUNCTION_HPP_

#include "LaGObase.hpp"

namespace LaGO {

/** A block from a larger function.
 */
class BlockFunction : public ReferencedObject {
private:
	/** The indices of the variables this function is defined for.
	 */
	vector<int> indices;
public:
	BlockFunction(const vector<int>& indices_)
	: indices(indices_)
	{ }
	
}; // class BlockFunction

} // namespace LaGO

#endif /*LAGOBLOCKFUNCTION_HPP_*/
