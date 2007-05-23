// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOConfig.h 94 2007-05-21 13:54:40Z stefan $

#ifndef LAGOSYMMATRIX_HPP_
#define LAGOSYMMATRIX_HPP_

#include "LaGObase.hpp"
#include "LaGOMatrix.hpp"

namespace LaGO {

class SymMatrix : public Matrix {
	public:
		SymMatrix(int dim)
		: Matrix(dim, dim)
		{ }	
	
}; // class SymMatrix
	
} // namespace LaGO

#endif /*LAGOSYMMATRIX_HPP_*/
