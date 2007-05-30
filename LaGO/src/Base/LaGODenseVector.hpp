// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGODENSEVECTOR_HPP_
#define LAGODENSEVECTOR_HPP_

#include "LaGObase.hpp"
#include "LaGORandomNumber.hpp"

#include "CoinDenseVector.hpp"

namespace LaGO {

class DenseVector : public CoinDenseVector<double>, public ReferencedObject {
public:
	DenseVector()
	{ }
	
	DenseVector(int size, double init_value=0.)
	: CoinDenseVector<double>(size, init_value)
	{ }

	DenseVector(int size, double* elements)
	: CoinDenseVector<double>(size, elements)
	{ }
	
	double operator()(int index) const {
		assert(index>=0 && index<getNumElements());
		return getElements()[index];
	}
	
	void setToBlock(const DenseVector& v, const vector<int>& indices) {
		unsigned int N=indices.size();
		resize(N);
		double* elem=getElements();
		for (unsigned int i=0; i<N; ++i, ++elem)
			*elem=v(indices[i]);
	}
	
	void setElementsOfBlock(const DenseVector& v, const vector<int>& indices) {
		int N=indices.size();
		assert(v.getNumElements()==N);
		const double* v_elem=v.getElements();
		for (int i=0; i<N; ++i, ++v_elem) {
			assert(indices[i]>=0 && indices[i]<getNumElements()); 
			getElements()[indices[i]]=*v_elem;
		}
	} 
	
	void setRandom(const DenseVector& lower, const DenseVector& upper);
	
	double operator*(const DenseVector& v) const {
		assert(getNumElements()==v.getNumElements());
		double ret=0;
		const double* x_=getElements();
		const double* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) ret+=*x_ * *v_;
		return ret;
	}
	
	DenseVector& operator+=(const DenseVector& v) {
		assert(getNumElements()==v.getNumElements());
		double* x_=getElements();
		const double* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) *x_+=*v_;
		return *this;
	}
};

template<class T>
ostream& operator<<(ostream& out, const CoinDenseVector<T>& v) {
	for (int i=0; i<v.getNumElements(); ++i)
		out << v.getElements()[i] << ' ';
	return out;
}
	
	
} // namespace LaGO


#endif /*LAGODENSEVECTOR_HPP_*/
