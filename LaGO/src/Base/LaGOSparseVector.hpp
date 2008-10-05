// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

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
	
	SparseVector(int n_val, int* indices, double* values)
	: CoinPackedVector(n_val, indices, values)
	{ }
	
	void scale(double factor);
	
	friend ostream& operator<<(ostream& out, const SparseVector& v) {
		const int* ind=v.getIndices();
		const double* el=v.getElements();
		for (int i=v.getNumElements(); i>0; --i, ++ind, ++el)
			out << '(' << *ind << ',' << *el << ')' << ' ';
		return out; 
	}
	
}; // class SparseVector

class SparseVectorCreator : public map<int,double> {
public:
	SmartPtr<SparseVector> getSparseVector() const {
		return new SparseVector(*this);
	}
	
	void insert(int index, double value) {
		map<int,double>::insert(pair<int, double>(index, value));
	}
	
	void add(const SparseVector& v);

	void addBlockVector(const SparseVector& v, const vector<int>& indices);
};	

	
} // namespace LaGO

#endif /*LAGOSPARSEVECTOR_HPP_*/