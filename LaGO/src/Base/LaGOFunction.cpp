// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOConfig.h 94 2007-05-21 13:54:40Z stefan $

#include "LaGOFunction.hpp"

namespace LaGO {
	
ostream& operator<<(ostream& out, FunctionEvaluationError& error) {
	out << error.className() << "::" << error.methodName() << ": " << error.fileName();
	return out;
}


}; // namespace LaGO
