// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOSMARTPTR_HPP_
#define LAGOSMARTPTR_HPP_

#include "LaGOConfig.h"
// get SmartPtr from Ipopt
#include "IpSmartPtr.hpp"

namespace LaGO {

typedef Ipopt::ReferencedObject ReferencedObject;
	
template<class T>
class SmartPtr : public Ipopt::SmartPtr<T> {
public:
	SmartPtr()
	: Ipopt::SmartPtr<T>()
	{ }
	
	SmartPtr(const SmartPtr<T>& copy)
	: Ipopt::SmartPtr<T>(copy)
	{ }
	
	SmartPtr(T* ptr)
	: Ipopt::SmartPtr<T>(ptr)
	{ }
	
	~SmartPtr() { }	
	
}; // class SmartPtr

	
} // namespace LaGO

#endif /*LAGOSMARTPTR_HPP_*/
