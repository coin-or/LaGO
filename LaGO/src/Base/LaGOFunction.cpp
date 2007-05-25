// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOFunction.hpp"

namespace LaGO {
	
ostream& operator<<(ostream& out, FunctionEvaluationError& error) {
	out << error.className() << "::" << error.methodName() << ": " << error.message();
	return out;
}


}; // namespace LaGO
