// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOMATRIX_HPP_
#define LAGOMATRIX_HPP_

#include "LaGObase.hpp"

namespace LaGO {

class Matrix : public ReferencedObject {
	friend ostream& operator<<(ostream& out, const Matrix& m) { m.print(out); return out; }
protected:
	int ncols;
	int nrows;
	
public:
	Matrix()
	: ncols(0), nrows(0)
	{ }
	
	Matrix(const int& ncols_, const int& nrows_)
	: ncols(ncols_), nrows(nrows_)
	{ }
	
	virtual ~Matrix() { }

	int getNumCols() const { return ncols; }
	int getNumRows() const { return nrows; }
	
	virtual void print(ostream& out) const=0;
	
	/** Performs the operation y = y + a * A*x
	 */
	virtual void addMultVector(DenseVector& y, const DenseVector& x, double a=1.) const {
		DenseVector z(getNumRows());
		multVector(z,x,a);
		y+=z;
	}
//	virtual void addMultVector(SparseVector& y, const SparseVector& x, double a=1.) const {
//		SparseVector z;
//		multVector(z,x,a);
//		y+=z;
//	}

	/** Performs the operation y = a * A*x
	 */
	virtual void multVector(DenseVector& y, const DenseVector& x, double a=1.) const=0;
//	virtual void multVector(SparseVector& y, const SparseVector& x, double a=1.) const=0;

	/** Performs the operation y^T*A*x.
	 */
	virtual double yAx(const DenseVector& y, const DenseVector& x) const {
		DenseVector z(getNumRows());
		multVector(z,x);
		return y*z;
	}	
//	virtual double yAx(const SparseVector& y, const SparseVector& x) const {
//		SparseVector z;
//		multVector(z,x);
//		return y*z;		
//	}
	
	/** Performs the operation x^T*A*x.
	 */
	virtual double xAx(const DenseVector& x) const {
		return yAx(x,x);
	}
//	virtual double xAx(const SparseVector& x) const {
//		return yAx(x,x);
//	}	
	
#ifdef COIN_HAS_FILIB
	virtual bool canIntervalEvaluation() const { return false; }
	
	virtual void addMultVector(IntervalVector& y, const IntervalVector& x, const interval<double>& a=interval<double>(1.)) const {
		IntervalVector z(getNumRows());
		multVector(z,x,a);
		y+=z;
	}

	virtual void multVector(IntervalVector& y, const IntervalVector& x, const interval<double>& a=interval<double>(1.)) const {
		throw CoinError("interval evaluation not possible", "multVector", "Matrix");
	}

	virtual interval<double> xAx(const IntervalVector& x) const {
		throw CoinError("interval evaluation not possible", "xAx", "Matrix");
	}
#endif
	
}; // class Matrix	
	
} // namespace LaGO

#endif /*LAGOMATRIX_HPP_*/
