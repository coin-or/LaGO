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
		if (*rowind<=*colind) operator[](pair<int,int>(*rowind,*colind))+=factor* *value;
		else operator[](pair<int,int>(*colind,*rowind))+=factor* *value;
		++rowind;
		++colind;
		++value;
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
  	if (*row<*col)
  		storage[*col*dim+*row]=*val;
  	else
  		storage[*row*dim+*col]=*val;
  	++row; ++col; ++val;
  }

  int info;
  Ipopt::IpLapackDsyev(eigvec!=NULL, dim, storage, dim, eigval.getElements(), info);

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
