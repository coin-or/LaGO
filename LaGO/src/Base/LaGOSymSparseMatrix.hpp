// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOSYMSPARSEMATRIX_HPP_
#define LAGOSYMSPARSEMATRIX_HPP_

#include "LaGObase.hpp"
#include "LaGOSymMatrix.hpp"

namespace LaGO {

class SymSparseMatrix;

class SymSparseMatrixCreator : public map<pair<int,int>, double> {
private:
	int dim;

public:
	SymSparseMatrixCreator(int dim_=0)
	: dim(dim_)
	{ }

	int getDim() const { return dim; }
	void setDim(int dim_) { dim=dim_; }

	/** Inserts element into matrix.
	 * If entry exists already, the new value is added to it.
	 * If row<col (in upper triangle), the indices are swaped.
	 */
	void insert(int row, int col, double value) {
		if (row<col) operator[](pair<int,int>(col,row))+=value;
		else operator[](pair<int,int>(row,col))+=value;
	}

	void add(double factor, const SymSparseMatrix& A);
	/** Adds A to our matrix, and also renames the indices.
	 */
	void addBlockMatrix(double factor, const SymSparseMatrix& A, const vector<int>& indices);

	void scale(double factor);

	/** Adds all elements in upper diagonal part into lower diagonal part.
	 */
	void cleanUpperDiagonal();

}; // SymSparseMatrixCreator

/** A symmetric matrix in triplet format.
 * Please give only A[i][j] or A[j][i], otherwise they are sum up.
 */
class SymSparseMatrix : public SymMatrix {
private:
  /** Number of nonzero-elements.
  */
	int nz;
	/** The values of the matrix.
	*/
	double* value;
	/** The row indices of the values.
	*/
	int* rowind;
	/** The column indices of the values.
	*/
	int* colind;

public:
	SymSparseMatrix(int dim=0)
	: SymMatrix(dim), nz(0), value(NULL), rowind(NULL), colind(NULL)
	{ }

	SymSparseMatrix(const SymSparseMatrixCreator& creator)
	: SymMatrix(creator.getDim()), nz(0), value(NULL), rowind(NULL), colind(NULL)
	{ set(creator);
	}

	~SymSparseMatrix();

	void set(const SymSparseMatrixCreator& creator);

	const double* getValues() const { return value; }
	const int* getRowIndices() const { return rowind; }
	const int* getColIndices() const { return colind; }
	int getNumNonzeros() const { return nz; }
	
	/** Multiplies all matrix elements with a value.
	 */
	void scale(double factor);

	void addMultVector(DenseVector& y, const DenseVector& x, double a=1.) const;
//	void addMultVector(SparseVector& y, const SparseVector& x, double a=1.) const;

	void multVector(DenseVector& y, const DenseVector& x, double a=1.) const {
		y.clear();
		addMultVector(y,x,a);
	}
//	void multVector(SparseVector& y, const SparseVector& x, double a=1.) const;

	double yAx(const DenseVector& y, const DenseVector& x) const;
//	double yAx(const SparseVector& y, const SparseVector& x) const;

	double xAx(const DenseVector& x) const;
//	double xAx(const SparseVector& x) const;

#ifdef COIN_HAS_FILIB
	bool canIntervalEvaluation() const { return true; }

	void addMultVector(IntervalVector& y, const IntervalVector& x, const interval<double>& a=interval<double>(1.)) const;

	void multVector(IntervalVector& y, const IntervalVector& x, const interval<double>& a=interval<double>(1.)) const {
		y.clear();
		addMultVector(y,x,a);
	}

	interval<double> xAx(const IntervalVector& x) const;

	interval<double> xAx_bx(const IntervalVector& x, const DenseVector& b) const;
#endif
	/** Computes eigenvalues and eigenvectors.
	 * @param eigval Storage for the eigenvalues.
	 * @param eigvec If not null, then this is storage for the eigenvectors.
	 * @return True, if computation succeeded. False otherwise.
	 */
	bool computeEigenValues(DenseVector& eigval, DenseVector* eigvec=NULL) const;

	/** Computes the minimal and maximal eigenvalue.
	 * @return True, if computation succeeded. False otherwise.
	 */
	bool computeMinMaxEigenValue(double& mineig, double& maxeig) const;

	void print(ostream& out) const;
};


} // namespace LaGO

#endif /*LAGOSYMSPARSEMATRIX_HPP_*/
