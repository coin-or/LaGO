// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGORandomNumber.hpp"

double gennor(double a, double b) { return 0; }
double genexp(double a) { return 0; }
double genunf(double a, double b) { return 0; }

namespace LaGO {

double getRandom(double lb, double ub) {
	if (lb<=-getInfinity() && ub>=getInfinity()) return gennor(0., 1.);
	if (lb<=-getInfinity()) return ub-genexp(1.);
	if (ub>= getInfinity()) return lb+genexp(1.);
	return genunf(lb, ub);
}

} // namespace LaGO
