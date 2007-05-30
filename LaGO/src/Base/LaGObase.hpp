// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOBASE_HPP_
#define LAGOBASE_HPP_

#include <exception>
#include <typeinfo>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>

// to get INFINITY
//#if !defined(__GNUC__)
//#define __USE_ISOC99 1
//#endif

#include <cmath>

//STL
#include <vector>
#include <map>
#include <set>
#include <list>
#include <string>

using namespace std;

#include "LaGOConfig.h"

#include "CoinHelperFunctions.hpp"
//#include "CoinSmartPtr.hpp"

#ifdef COIN_HAS_FILIB
#define FILIB_EXTENDED
#include "interval/interval.hpp"
using filib::interval;
#endif

#include "LaGOSmartPtr.hpp"
#include "LaGODenseVector.hpp"
#include "LaGOSparseVector.hpp"
#include "LaGOFunction.hpp"
#ifdef COIN_HAS_FILIB
#include "LaGOIntervalVector.hpp"
#endif

namespace LaGO {

inline double getInfinity() { return 1E+300; }

//inline double relDistance(const double& value1, const double& value2, double norm, const double& invalid_norm=1E+14) {
//	if (value1==value2) return 0.;
//	if (norm<0) norm=-norm;
//	if (norm<invalid_norm)
//		return (value2-value1)/(1+norm);
//	norm=CoinMin(fabs(value1), fabs(value2));
//	if (norm<invalid_norm)			
//		return (value2-value1)/(1+norm);
//	return 1.;
//}
//
//inline double relDistance(const double& value1, const double& value2) {
//	return relDistance(value1, value2, CoinMin(fabs(value1), fabs(value2)));
//}
//

} // namespace LaGO

#endif /*LAGOBASE_HPP_*/
