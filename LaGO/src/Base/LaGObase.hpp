// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: base.hpp 94 2007-05-21 13:54:40Z stefan $

#ifndef LAGOBASE_HPP_
#define LAGOBASE_HPP_

#include <exception>
#include <typeinfo>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>

// to get INFINITY
#if !defined(__GNUC__)
#define __USE_ISOC99 1
#endif

#include <cmath>

//STL
#include <vector>
#include <map>
#include <set>
#include <list>
#include <string>

using namespace std;

#include "LaGOConfig.h"
//
//#include "DenseVector.hpp"
//#include "SparseVector.hpp"
//#include "Timer.hpp"

#include "CoinHelperFunctions.hpp"
//#include "CoinSmartPtr.hpp"

#include "LaGOSmartPtr.hpp"
#include "LaGODenseVector.hpp"
#include "LaGOSparseVector.hpp"
#include "LaGOFunction.hpp"


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
inline int random(int lb, int ub) {
   return lb+(int)((ub-lb+0.99) * CoinDrand48()); // as suggested in the rand()-manual
}

} // namespace LaGO

#endif /*LAGOBASE_HPP_*/
