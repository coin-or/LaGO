// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOConfig.h 94 2007-05-21 13:54:40Z stefan $

#include "LaGOSymSparseMatrix.hpp"

namespace LaGO {

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

void SymSparseMatrix::print(ostream& out) const {
	const int* rowind_=rowind;
	const int* colind_=colind;
	const double* value_=value;
	
	for (int i=nz; i>0; --i)
		out << '(' << *++rowind_ << ',' << *++colind_ << ")=" << *++value_ << ' ';
}	
	
} // namespace LaGO
