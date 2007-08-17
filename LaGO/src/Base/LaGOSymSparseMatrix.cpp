// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSymSparseMatrix.hpp"
// for eigenvalue computation.. just let Ipopt handle the Lapack interface :-)
#include "IpLapack.hpp"

namespace LaGO {

void SymSparseMatrixCreator::add(double factor, const SymSparseMatrix& A) {
	assert(getDim()==A.getNumCols());
	const int* rowind=A.getRowIndices();
	const int* colind=A.getColIndices();
	const double* value=A.getValues();
	for (int i=A.getNumNonzeros(); i>0; --i) {
		if (*rowind<*colind) operator[](pair<int,int>(*colind,*rowind))+=factor* *value;
		else operator[](pair<int,int>(*rowind,*colind))+=factor* *value;
		++rowind;
		++colind;
		++value;
	}
}

/** Adds A to our matrix, and also renames the indices.
 */
void SymSparseMatrixCreator::addBlockMatrix(double factor, const SymSparseMatrix& A, const vector<int>& indices) {
	assert((int)indices.size()==A.getNumCols());
	const int* rowind=A.getRowIndices();
	const int* colind=A.getColIndices();
	const double* value=A.getValues();
	for (int i=A.getNumNonzeros(); i>0; --i) {
		int row=indices[*rowind];
		int col=indices[*colind];
		assert(row<getDim() && col<getDim());
		if (row<col) operator[](pair<int,int>(col,row))+=factor* *value;
		else operator[](pair<int,int>(row,col))+=factor* *value;
		++rowind;
		++colind;
		++value;
	}
}

void SymSparseMatrixCreator::scale(double factor) {
	if (factor==0.) clear();
	else
		for (iterator it(begin()); it!=end(); ++it)
			it->second*=factor;
}


void SymSparseMatrixCreator::cleanUpperDiagonal() {
	iterator it(begin());
	while (it!=end()) {
		if (it->first.first<it->first.second) { // row < col
			iterator next(it); ++next;
			operator[](pair<int,int>(it->first.second,it->first.first))+=it->second;
			erase(it);
			it=next;
		} else
			++it;
	}
}


SymSparseMatrix::~SymSparseMatrix() {
	delete[] value;
	delete[] rowind;
	delete[] colind;
}

void SymSparseMatrix::set(const SymSparseMatrixCreator& creator) {
	delete[] value;
	delete[] rowind;
	delete[] colind;
	nz=0;
	value=new double[creator.size()];
	rowind=new int[creator.size()];
	colind=new int[creator.size()];

	for (SymSparseMatrixCreator::const_iterator it(creator.begin()); it!=creator.end(); ++it) {
		if (it->second==0) continue;
		rowind[nz]=it->first.first;
		colind[nz]=it->first.second;
		value[nz]=it->second;
		++nz;
	}
}

void SymSparseMatrix::addMultVector(DenseVector& y, const DenseVector& x, double a) const {
	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;

	for (int i=nz; i>0; --i) {
		if (*rowind_==*colind_) y[*rowind_]+=a * *value_ * x[*colind_];
		else {
			y[*rowind_]+=a * *value_ * x[*colind_];
			y[*colind_]+=a * *value_ * x[*rowind_];
		}
		++rowind_;
		++colind_;
		++value_;
	}
}

double SymSparseMatrix::yAx(const DenseVector& y, const DenseVector& x) const {
	double ret=0;
	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;

	for (int i=nz; i>0; --i) {
		if (*rowind_==*colind_) ret+=y[*rowind_]* *value_ * x[*colind_];
		else ret+=*value_ * (y[*rowind_]*x[*colind_]+x[*rowind_]*y[*colind_]);
		++rowind_;
		++colind_;
		++value_;
	}

	return ret;
}

double SymSparseMatrix::xAx(const DenseVector& x) const {
	double ret=0;
	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;

	for (int i=nz; i>0; --i) {
		if (*rowind_==*colind_) {
			double x_=x[*rowind_];
			if (x_) ret+=x_ * *value_ * x_;
		} else {
			ret+=2* *value_ * x[*rowind_] * x[*colind_];
		}
		++rowind_;
		++colind_;
		++value_;
	}

	return ret;
}

#ifdef COIN_HAS_FILIB
// y = y + a * Ax
void SymSparseMatrix::addMultVector(IntervalVector& y, const IntervalVector& x, const interval<double>& a) const {
	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;

	if (a.isPoint()) {
		for (int i=nz; i>0; --i) {
			if (*rowind_==*colind_) y[*rowind_]+=a* *value_ * x[*colind_];
			else {
				y[*rowind_]+=a* *value_ * x[*colind_];
				y[*colind_]+=a* *value_ * x[*rowind_];
			}
			++rowind_;
			++colind_;
			++value_;
		}
	} else {
		IntervalVector tmp(x.getNumElements());
		addMultVector(tmp, x, interval<double>(1.)); // Ax;  for this call, a=1. is a point interval
		y.addVector(a, tmp);
	}
}

interval<double> SymSparseMatrix::xAx(const IntervalVector& x) const {
	double ret=0;
	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;

	for (int i=nz; i>0; --i) {
		if (*rowind_==*colind_) {
			ret+=*value_ * sqr(x[*rowind_]);
		} else {
			ret+=(2* *value_) * (x[*rowind_] * x[*colind_]);
		}
		++rowind_;
		++colind_;
		++value_;
	}

	return ret;
}

/** Computes a tight interval for a*x+b where x is an interval (and a and b are not).
 */
interval<double> envelope(double a, double b, const interval<double>& x) {
	double min=filib::fp_traits<double>::ninfinity();
	double max=filib::fp_traits<double>::infinity();
	double extreme=-b/(2*a);
	if (extreme>=x.inf() && extreme<=x.sup()) { // extreme value in interval
		if (a>0) min=-b*b/(4*a);
		else 	max=-b*b/(4*a);
	}
	double atleft= (a*x.inf()+b)*x.inf();
	double atright=(a*x.sup()+b)*x.sup();
	if (atleft<min) min=atleft;
	if (atleft>max) max=atleft;
	if (atright<min) min=atright;
	if (atright>max) max=atright;
	return interval<double>(min, max);
}

interval<double> SymSparseMatrix::xAx_bx(const IntervalVector& x, const DenseVector& b) const {
	interval<double> ret(0.);
	const interval<double>& zero(interval<double>::ZERO());

	IntervalVector coeff(x.getNumElements());
	DenseVector diag(x.getNumElements());

	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;
	for (int i=nz; i>0; --i) {
		if (*value_)
			if (*rowind_==*colind_)
				diag[*rowind_]+=*value_;
			else
				coeff[*rowind_]+=(2* *value_) * x[*colind_];
		++rowind_;
		++colind_;
		++value_;
	}

	interval<double>* coeff_=coeff.getElements();
	double* diag_=diag.getElements();
	const interval<double>* x_=x.getElements();
	const double* b_=b.getElements();
	for (int i=0; i<x.getNumElements(); ++i, ++coeff_, ++diag_, ++x_, ++b_) {
		if (*x_==zero) continue;
		if (x_->isPoint()) ret+=x_->inf()*(*coeff_+*diag_*x_->inf()+*b_);
		else {
//			clog << *x_ << '*' << *coeff_ << '=' << *x_ * *coeff_ << '\t' << mult(*x_,*coeff_) << endl;
			ret+=*x_ * *coeff_ + envelope(*diag_, *b_, *x_);
		}
	}
	return ret;
}
#endif // COIN_HAS_FILIB

bool SymSparseMatrix::computeEigenValues(DenseVector& eigval, DenseVector* eigvec) const {
	int dim=getNumCols();
	if (eigval.getNumElements()<dim) eigval.resize(dim);

	double* storage;
	if (eigvec) {
		if (eigvec->getNumElements()<dim*dim)
			eigvec->resize(dim*dim);
		storage=eigvec->getElements();
	} else {
		storage=new double[dim*dim];
		CoinZeroN(storage, dim*dim);
	}

  // we copy the content of the matrix into storage
  const double* val=value;
  const int* row=rowind;
  const int* col=colind;
  for (int i=nz; i>0; --i) {
//  	if (*row<*col)
  		storage[*col*dim+*row]=*val;
//  	else
  		storage[*row*dim+*col]=*val;
  	++row; ++col; ++val;
  }
//  clog << "matrix:";
//  for (int i=0; i<dim*dim; ++i)
//	 	clog << ' ' << storage[i];
//	clog << endl;

  int info;
  Ipopt::IpLapackDsyev(eigvec!=NULL, dim, storage, dim, eigval.getElements(), info);
//  clog << "eigenvalues:";
//  for (int i=0; i<dim; ++i)
//  	clog << ' ' << eigval[i];
//  clog << endl;

  if (!eigvec) delete[] storage;

  return (info==0);
}

bool SymSparseMatrix::computeMinMaxEigenValue(double& mineig, double& maxeig) const {
	DenseVector eigval(getNumCols());
	bool ret=computeEigenValues(eigval);
	if (!ret) return false;
	double* e=eigval.getElements();
	mineig=*e;
	maxeig=*e;
	++e;
	for (int i=getNumCols(); i>1; --i, ++e) {
		if (*e<mineig) mineig=*e;
		else if (*e>maxeig) maxeig=*e;
	}
	return true;
}

void SymSparseMatrix::print(ostream& out) const {
	if (!nz) return;
	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;

	out << "SymSparseMatrix:";
	for (int i=nz; i>0; --i)
		out << " (" << *rowind_++ << ',' << *colind_++ << ")=" << *value_++;
}

} // namespace LaGO
