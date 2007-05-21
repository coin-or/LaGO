// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.

// $Id: base.hpp 18 2007-04-02 16:25:41Z stefan $

#ifndef BASE_HPP_
#define BASE_HPP_

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

using namespace std;

#include "LaGOConfig.h"
//#include "Referenced.hpp"
//#include "SmartPtr.hpp"
//
//#include "DenseVector.hpp"
//#include "SparseVector.hpp"
//#include "Timer.hpp"

#include "CoinHelperFunctions.hpp"

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

#endif /*BASE_HPP_*/
