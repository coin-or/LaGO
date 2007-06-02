// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef OSISOLVER_HPP_
#define OSISOLVER_HPP_

#include "LaGObase.hpp"

#ifdef COIN_HAS_CPX
#define OsiXXXSolverInterface OsiCpxSolverInterface
#include "OsiCpxSolverInterface.hpp"
#else
#define OsiXXXSolverInterface OsiClpSolverInterface
#include "OsiClpSolverInterface.hpp"
#endif

namespace LaGO {

class OsiSolver : public OsiXXXSolverInterface, public ReferencedObject {
public:
	OsiSolver();

}; // class OsiSolver
	
} // namespace LaGO

#endif /*OSISOLVER_HPP_*/
