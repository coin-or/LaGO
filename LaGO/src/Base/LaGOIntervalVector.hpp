// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOINTERVALVECTOR_HPP_
#define LAGOINTERVALVECTOR_HPP_

#include "LaGObase.hpp"

#include "CoinDenseVector.hpp"
#include "CoinPackedVector.hpp"

template <class Type>
inline bool operator<(const interval<Type>& x, const interval<Type>& y) { return x.slt(y); }
template <class Type>
inline bool operator>(const interval<Type>& x, const interval<Type>& y) { return y.slt(x); }
template <class Type>
inline Type& operator+=(Type& x, const interval<Type>& y) { 
	throw CoinError("Type += interval<Type> is not defined", "Type& operator+=(Type&, const interval<Type>&)", "");
	return x;
}
template <class Type>
inline interval<Type> CoinAbs(const interval<Type> value) {
	throw CoinError("CoinAbs not defined for an interval", "interval<Type> CoinAbs(const interval<Type> value)", "");
  return interval<Type>(0);
}
template <class Type>
inline interval<Type> CoinMax(const interval<Type> value) {
	throw CoinError("CoinMax not defined for an interval", "interval<Type> CoinMax(const interval<Type> value)", "");
  return interval<Type>(0);
}

namespace LaGO {

inline double translateNInftyCoin2Filib(const double& val) {
	if (val==-getInfinity()) return filib::fp_traits<double>::ninfinity(); 
	return val;
}
inline double translatePInftyCoin2Filib(const double& val) {
	if (val==getInfinity()) return filib::fp_traits<double>::infinity(); 
	return val;
}
inline double translateNInftyFilib2Coin(const double& val) {
	if (val==filib::fp_traits<double>::ninfinity()) return -getInfinity(); 
	return val;
}
inline double translatePInftyFilib2Coin(const double& val) {
	if (val==filib::fp_traits<double>::infinity()) return getInfinity(); 
	return val;
}

// TODO: multiplication of infinite intervals or 0*infinite interval

//#undef FILIB_UPWARD_MULTIPLIES
//#undef FILIB_DOWNWARD_MULTIPLIES
//#undef FILIB_RESET
//#define FILIB_UPWARD_MULTIPLIES(R,A,B)		R=filib::fp_traits<double>::upward_multiplies(A,B)
//#define FILIB_DOWNWARD_MULTIPLIES(R,A,B)	R=filib::fp_traits<double>::downward_multiplies(A,B)
//#define FILIB_RESET filib::fp_traits<double>::reset()

class IntervalVector : public CoinDenseVector<interval<double> >, public ReferencedObject {
public:
	IntervalVector()
	{ }
	
	IntervalVector(int size, interval<double> init_value=interval<double>(0.))
	: CoinDenseVector<interval<double> >(size, init_value)
	{ }

	IntervalVector(const CoinDenseVector<double>& low, const CoinDenseVector<double>& up)
	: CoinDenseVector<interval<double> >(low.getNumElements())
	{ assert(up.getNumElements()==low.getNumElements());
		const double* low_=low.getElements();
		const double* up_=up.getElements();
		interval<double>* el=getElements();
		for (int i=getNumElements(); i>0; --i, ++low_, ++up_, ++el)
			*el=interval<double>(translateNInftyCoin2Filib(*low_),translatePInftyCoin2Filib(*up_));	
	}
	
	IntervalVector(const IntervalVector& v, const vector<int>& indices)
	: CoinDenseVector<interval<double> >(indices.size())
	{	interval<double>* el=getElements();
		for (int i=0; i<getNumElements(); ++i, ++el) {
			assert(indices[i]>=0 && indices[i]<v.getNumElements());
			*el=v(indices[i]);
		}
	}	
	
	~IntervalVector() { }
	
	interval<double> operator()(int index) const {
		assert(index>=0 && index<getNumElements());
		return getElements()[index];
	}
	
	void setToBlock(const IntervalVector& v, const vector<int>& indices) {
		unsigned int N=indices.size();
		resize(N);
		interval<double>* elem=getElements();
		for (unsigned int i=0; i<N; ++i, ++elem)
			*elem=v(indices[i]);
	}
	
	void setElementsOfBlock(const IntervalVector& v, const vector<int>& indices) {
		int N=indices.size();
		assert(v.getNumElements()==N);
		const interval<double>* v_elem=v.getElements();
		for (int i=0; i<N; ++i, ++v_elem) {
			assert(indices[i]>=0 && indices[i]<getNumElements()); 
			getElements()[indices[i]]=*v_elem;
		}
	} 
	
	interval<double> operator*(const IntervalVector& v) const {
		assert(getNumElements()==v.getNumElements());
		interval<double> ret(0.);
		const interval<double>* x_=getElements();
		const interval<double>* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) ret+=*x_ * *v_;
		return ret;
	}

	interval<double> operator*(const CoinDenseVector<double>& v) const {
		assert(getNumElements()==v.getNumElements());
		interval<double> ret(0.);
		const interval<double>* x_=getElements();
		const double* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_)
			if (*v_) ret+=*x_ * *v_;
		return ret;
	}

	interval<double> operator*(const CoinPackedVector& v) const {
		interval<double> ret(0.);
		const int* ind=v.getIndices();
		const double* el=v.getElements();
		for (int i=v.getNumElements(); i>0; --i, ++ind, ++el) {
			assert(*ind>=0 && *ind<getNumElements());
			if (*el) ret+=*el*getElements()[*ind];
		}
		return ret;
	}
	
	IntervalVector& operator+=(const IntervalVector& v) {
		assert(getNumElements()==v.getNumElements());
		interval<double>* x_=getElements();
		const interval<double>* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) *x_+=*v_;
		return *this;
	}
	
	void addVector(const interval<double>& factor, const IntervalVector& v) {
		if (factor==interval<double>::ZERO()) return;
		if (factor==interval<double>::ONE()) { operator+=(v); return; }
		assert(getNumElements()==v.getNumElements());
		interval<double>* x_=getElements();
		const interval<double>* v_=v.getElements();
		for (int i=getNumElements(); i>0; --i, ++x_, ++v_) *x_+=(factor* *v_);
	}
};
	
	
} // namespace LaGO


#endif /*LAGOINTERVALVECTOR_HPP_*/
