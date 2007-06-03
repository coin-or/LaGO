// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGODENSEVECTOR_HPP_
#define LAGODENSEVECTOR_HPP_

#include "LaGObase.hpp"
#include "LaGORandomNumber.hpp"

#include "CoinDenseVector.hpp"
#include "CoinPackedVector.hpp"

namespace LaGO {

class DenseVector : public CoinDenseVector<double>, public ReferencedObject {
public:
	DenseVector()
	{ }
	
	DenseVector(int size, double init_value=0.)
	: CoinDenseVector<double>(size, init_value)
	{ }

	DenseVector(int size, const double* elements)
	: CoinDenseVector<double>(size, elements)
	{ }
	
	DenseVector(int size, const CoinPackedVector& v)
	: CoinDenseVector<double>(size)
	{ const int* ind=v.getIndices();
		const double* el=v.getElements();
		for (int i=v.getNumElements(); i>0; --i, ++ind, ++el) {
			assert(*ind>=0 && *ind<size);
			getElements()[*ind]=*el; 
		}
	}
	
	DenseVector(const DenseVector& v, const vector<int>& indices)
	{ setToBlock(v, indices);
	}
	
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

	double operator*(const CoinPackedVector& v) const {
		double ret=0;
		const double* v_el=v.getElements();
		const int* v_ind=v.getIndices();
		for (int i=v.getNumElements(); i>0; --i, ++v_el, ++v_ind) {
			assert(0<=*v_ind && *v_ind<getNumElements());
			ret+=getElements()[*v_ind] * *v_el;
		}
		return ret;
	}
	
	DenseVector& operator+=(const DenseVector& v) {
		assert(getNumElements()==v.getNumElements());
		double* x_=getElements();
		const double* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) *x_+=*v_;
		return *this;
	}
	
	DenseVector& operator+=(const CoinPackedVector& v) {
		const double* v_el=v.getElements();
		const int* v_ind=v.getIndices();
		for (int i=v.getNumElements(); i>0; --i, ++v_el, ++v_ind) {
			assert(0<=*v_ind && *v_ind<getNumElements());
			getElements()[*v_ind]+=*v_el;
		}
		return *this;
	}

	void addVector(double factor, const DenseVector& v) {
		if (factor==1.) operator+=(v);
		if (factor==0.) return;
		assert(getNumElements()==v.getNumElements());
		double* x_=getElements();
		const double* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) *x_+=factor**v_;
	}
	
	double euclidianDistance(const DenseVector& v) const {
		assert(getNumElements()==v.getNumElements());
		double ret=0.;
		const double* x_=getElements();
		const double* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) ret+=(*x_-*v_)*(*x_-*v_);
		return sqrt(ret);
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
