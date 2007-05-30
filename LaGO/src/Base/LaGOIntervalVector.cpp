// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOIntervalVector.hpp"

// we have to include also CoinDenseVector.cpp and explicitely define the class CoinDenseVector<interval<double> > since
// some definition of template methods are otherwise not compiled for the interval<double> type
#include "CoinDenseVector.cpp"
template class CoinDenseVector<interval<double> >;

