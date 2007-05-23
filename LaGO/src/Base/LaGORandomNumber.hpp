// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGORANDOMNUMBER_HPP_
#define LAGORANDOMNUMBER_HPP_

#include "LaGObase.hpp"

namespace LaGO {

inline int getRandom(int lb, int ub) {
   return lb+(int)((ub-lb+0.99) * CoinDrand48()); // as suggested in the rand()-manual
}

double getRandom(double lb, double ub);

} // namespace LaGO

#endif /*LAGORANDOMNUMBER_HPP_*/
