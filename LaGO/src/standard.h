// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

// standard include files
#ifndef STANDARD_H
#define STANDARD_H

#include "LaGOConfig.h"

#ifndef __INTEL_COMPILER
#include <ansidecl.h> // to get GCC_VERSION
#endif
#include <cstdlib>
#include <cstdio>
#include <cmath>
#ifdef HAVE_CFLOAT
#include <cfloat>
#endif
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
// #ifndef finite
// #define finite MY_C_FINITE
// #endif


inline double sqrt(int x) { return sqrt((double)x); }
#ifndef _WIN32
inline long pow(int a, unsigned int b) { return (long)pow((double)a, (double)b); }
inline long pow(int a, int b) { return (long)pow((double)a, (double)b); }
#endif

#include <cassert>

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#ifdef FILIB_AVAILABLE
// must be included before tools.h, because it has its own MIN and MAX routines
#define FILIB_EXTENDED
#include "interval/interval.hpp"
using filib::interval;
template <class Type>
inline bool operator<(const interval<Type>& x, const interval<Type>& y) { return x.slt(y); };
#endif

#include "tools.h"
#include "uservector.h"

#endif // STANDARD_H

