// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: base.hpp 94 2007-05-21 13:54:40Z stefan $

#ifndef LAGOSPARSEVECTOR_HPP_
#define LAGOSPARSEVECTOR_HPP_

#include "LaGObase.hpp"
#include "CoinPackedVector.hpp"

namespace LaGO {

class SparseVectorCreator;	
	
class SparseVector : public ReferencedObject, public CoinPackedVector {
public:
	SparseVector()
	: CoinPackedVector(false)
	{}
	
	SparseVector(const SparseVectorCreator& creator);

	/** A sparsevector containing only one element.
	 */	
	SparseVector(int index, double value)
	: CoinPackedVector(1, &index, &value, false)
	{ }
	
}; // class SparseVector

class SparseVectorCreator : public map<int,double> {
public:
	SmartPtr<SparseVector> getSparseVector() const {
		return new SparseVector(*this);
	}
	
	void insert(int index, double value) {
		map<int,double>::insert(pair<int, double>(index, value));
	}
};	

	
} // namespace LaGO

#endif /*LAGOSPARSEVECTOR_HPP_*/
