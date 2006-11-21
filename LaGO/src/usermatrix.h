// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

// usermatrix.h

#ifndef USERMATRIX_H
#define USERMATRIX_H

// LAGO-stuff:
#include "standard.h"
#include "param.h"

#ifdef LANCZOS_AVAILABLE
// Helmberg-stuff:
#include <hel_tools/matrix.h>
#include <hel_tools/lanczpol.h>
#else
typedef int Integer;
typedef double Real;
inline double d_sign(double a,double b) { return (((b>=0)&&(a>=0))||((b<0)&&(a<0)))?a:-a; }
#endif

// TNT-stuff:
#include <tnt_cmat.h>
#include <tnt_vec.h>

/** Abstract class for a user defined quadratic matrix.
    Following methods are abstract:
    - MultV(dvector&, dvector&)
    - solve(dvector&, dvector&)
*/
class UserMatrix
#ifdef LANCZOS_AVAILABLE
: public Lanczosmatrix
#endif
{
      /** Prints out information about the matrix.
          Calls a.print(out).
          @param out The ostream to print to.
          @param a The UserMatrix to get informations about.
          @return The ostream.
          @see print(ostream&)
      */
      friend ostream& operator << (ostream& out, const UserMatrix& a)
      {  a.print(out); return out; };

   protected:
      /** The dimension of this matrix.
      */
      int dim_;

   public:
      /** (Standard-)Constructor for the dimension.
          @param n The dimension, default is 0.
      */
      UserMatrix(int n=0)
      : dim_(n)
      { }

      /** Virtual Destructor.
      */
      virtual ~UserMatrix() { }

      /** Gives the dimension.
          @return The dimension as int.
      */
      int dim() const { return dim_; };

      /** Multiplies with a UserVector<double>.
          Abstract.
          @param y The UserVector<double> to store the result in: y = *this * x.
          @param x The UserVector<double> to multiply with.
          @see MultV(double*, double*)
      */
      virtual void MultV(UserVector<double>& y,const UserVector<double>& x) const=0;

      /** Multiplies the matrix with a double*.
          @param y The double* to store the result in.
          @param x The double* to multiply with this matrix.
          @see MultV(dvector&, const dvector&)
      */
      virtual void MultV(double* y, const double* x) const {
        dvector xx(x, dim()), yy(dim());
        MultV(yy,xx);
        for (int i=0; i<dim(); i++) y[i]=yy[i];
      }

      /** Multiplies from left and from right with UserVector<double>'s.
          @param y The UserVector<double> to multiply with from left.
          @param x The UserVector<double> to multiply with from right.
          @return The double value y*this*x.
          @see xAx(const UserVector<double>&)
      */
      virtual double yAx(const UserVector<double>& y, const UserVector<double>& x) const {
        Pointer<UserVector<double> > z(y.getemptycopy());
        MultV(*z,x);
        return *z*y;
      }

      /** Multiplies the quadratic product of this matrix and a UserVector<double>.
          @param x The UserVector<double> to multiply with.
          @return The double yAx(x,x).
          @see yAx(const UserVector<double>&, const UserVector<double>&)
      */
      virtual double xAx(const UserVector<double>& x) const
      {  return yAx(x,x);  };

      /** Adds to a UserVector<double> the product of a double value and this matrix, multiplied with a UserVector<double>.
          @param y The UserVector<double> to store the result in: y + val * A * x
          @param x The UserVector<double> to multiply this matrix with.
          @param val The double to multiply with.
      */
      virtual void AddMult(UserVector<double>& y, const UserVector<double>& x, const double val) const {
        y.AddMult(val, (*this * x));
      };

      /** Multiplies this matrix with a dvector.
          @param x The dvector to multiply with.
          @return The result (this*x).
          @see MultV(dvector&, const dvector&)
      */
      dvector operator*(const dvector& x) const {
        dvector y(x.size());
        MultV(y, x);
        return y;
      }

      /** Multiplies this matrix with a double*.
          @param x The double* to multiply with.
          @return The result as a new(!!) double*.
          @see MultV(double*, double*)
      */
      double* operator*(const double* x) const {
        double* y=new double[dim()];
        MultV(y,x);
        return y;
      }

      Pointer<UserVector<double> > operator*(const Pointer<UserVector<double> > x) const {
      	Pointer<UserVector<double> > y(x->getemptycopy());
      	MultV(*y,*x);
      	return y;
      }

#ifdef LANCZOS_AVAILABLE
      /** The dimension.
          @return The dimension of this matrix.
          @see dim()
      */
      Integer lanczosdim() const { return dim(); }

      /** Multiplies (-1) times this matrix with some vectors.
          y = - *this * x.
          @param x The vectors to multiply with.
          @param y The vectors to store the results in.
          @return 0, if all went right.
          @see MultMv(double*, double*)
      */
      int lanczosmult(const Matrix& x, Matrix &y) const {
        for (int i=0; i<x.coldim(); i++)
          ((UserMatrix*)this)->MultV(y.get_store()+i*x.rowdim(), x.get_store()+i*x.rowdim());
        y*=-1;
        return 0;
      }

      /** Calculates the eigenvalues with the lanczos-algorithmus.
          Calculates the first eig_vec.size() eigenvalues, starting with the smallest one.
          @param eig_vec A vector of dvector's to store the eigenvectors in.
          @param eig_val A vector of double's to store the eigenvalues in.
          @param param Some parameters for lanczos, default is NULL.
          @return The error-code from lanczos.
          @see eig_lanczos(dvector&, double&)
      */
      int eig_lanczos(vector<dvector>& eig_vec, vector<double>& eig_val, Param *param=NULL) const;

      /** Calculates the minimum eigenvalue and eigenvector with the lanczos-algorithmus.
          @param eig_vec A dvector to store the eigenvector in.
          @param eig_val A double to store the minimum eigenvalue in.
          @param param Some parameters for lanczos, default is NULL.
          @return The error-code form lanczos.
          @see eig_lanczos(vector<dvector>&, vector<double>&)
      */
      int eig_lanczos(dvector& eig_vec, double& eig_val, Param *param=NULL) const;
#endif

#ifdef ARPACK_AVAILABLE
      /** Multiplication for use this matrix with ARPACK.
					@param x The vector to multiply with the matrix.
					@param y The vector to store the result in.
					@see MultV(double*, double*)
			*/
      void arpackmult(double* x, double* y) { MultV(y,x); }

      /** Calculates the minimum eigenvalue and eigenvector with ARPACK.
          @param eig_vec A dvector to store the eigenvector in.
          @param eig_val A double to store the minumum eigenvalue in.
          @param param Some parameters, default is NULL.
          @return The number of converged eigenvalues, should be at least 1.
      */
      int eig_arpack(dvector& eig_vec, double& eig_val, Param* param=NULL);

      /** Calculates the smallest eigenvalues and eigenvectors with ARPACK.
          @param eig_vec A vector of dvectors to store the eigenvectors in.
          @param eig_val A vector of doubles to store the eigenvalues in.
          @param param Some parameters, default is NULL.
          @return The number of converges eigenvalues, should be at least eig_vec.size().
			*/
      int eig_arpack(vector<dvector>& eig_vec, vector<double>& eig_val, Param* param=NULL);
#endif

      /** Calculates the smallest eigenvalues and eigenvectors with the QL-algorithm.
          Converts the matrix to a DenseMatrix and calls eig_ql from there.
          @param eig_vec A vector of dvectors to store the eigenvectors in.
          @param eig_val A vector of doubles to store the eigenvalues in.
          @param param Some parameters, default is NULL.
	        @return The return value from imtql2, which is 0, if all went good or the number of the eigenvalue, which wasn't founds after 30 iterations.
          @see DenseMatrix::eig_ql(vector<dvector>&, vector<double>&, Param*)
      */
			virtual int eig_ql(vector<dvector> &eig_vec, vector<double> &eig_val) const;

      /** Calculates the smallest eigenvalues of this matrix.
          If dim>=50, calls the lanczos-algorithm.
          Otherwise a DenseMatrix is created and the ql-algorithm called.
          @param eig_vec The dvectors to store the eigenvectors in.
          @param eig_val The doubles to store the eigenvalues in.
          @param param Parameters for Lanczos or ARPACK, default is NULL.
          @return Depends on the choosen eigenvalue algorithm. It's 0, if all went good.
          @see eig_ql(vector<dvector>&, vector<double>&)
          @see eig_lanczos(vector<dvector>&, vector<double>&, Param*)
          @see eig_arpack(vector<dvector>&, vector<double>&, Param*)
      */
      int eig(vector<dvector>& eig_vec, vector<double>& eig_val, Param *param=NULL);

      /** Calculates the smallest eigenvalue of this matrix.
          @param eig_vec The dvector to store the smallest eigenvector in.
          @param eig_val The double to store the smallest eigenvalue in.
          @param param Parameters for the lanczos-algorithm, default is NULL.
          @return Depends on the choosen eigenvalue algorithm. It's 0, if all went good.
          @see eig(vector<dvector>&, vector<double>&, Param*)
      */
      int eig(dvector& eig_vec, double& eig_val, Param *param=NULL);

      /** Prints out this matrix.
          Generates a DenseMatrix with the values of this matrix and prints this one.
          @param out The ostream to print to.
          @see DenseMatrix::print(ostream&)
      */
      virtual void print(ostream &out) const;
};

#ifdef FILIB_AVAILABLE
class IntervalCompliantMatrix : public UserMatrix {
	public:
		IntervalCompliantMatrix(int n=0)
		: UserMatrix(n)
		{ }

		virtual void MultV(IntervalVector& y, const IntervalVector& x) const=0;

		using UserMatrix::MultV;

		virtual interval<double> yAx(const IntervalVector& y, const IntervalVector& x) const {
			IntervalVector Ax(x.dim()); this->MultV(Ax,x);
			return y*Ax;
		}
		using UserMatrix::yAx;

		virtual interval<double> xAx(const IntervalVector& x) const { return yAx(x, x); }
		using UserMatrix::xAx;
};
#endif

/** A UserMatrix with an (int,int)-operator to read elements of the matrix.
*/
class ExtUserMatrix
#ifdef FILIB_AVAILABLE
: public IntervalCompliantMatrix
#else
: public UserMatrix
#endif
{
  public:
    /** (Default-)Constructor.
        @param n The dimension of the matrix, default is 0.
    */
    ExtUserMatrix(int n=0)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix(n)
#else
		: UserMatrix(n)
#endif
    { }

    /** Operator to read one element of the matrix.
        Abstract.
        @param row The row of the element to read.
        @param col The col of the element to read.
        @return The element at index (row, col).
    */
    virtual double operator()(int row, int col) const=0;
};

/** Wrapper class to shift a matrix.
    *this == *A + diag(shift).
*/
class ShiftMatrix : public UserMatrix {
  protected:
    /** The matrix to shift.
    */
    Pointer<UserMatrix> A;

    /** The shift-value.
    */
    dvector shift;

  public:
    /** Constructor for a UserMatrix and a shift-value.
        @param A_ The UserMatrix to shift.
        @param shift_ The shift-value.
    */
    ShiftMatrix(Pointer<UserMatrix> A_, const dvector& shift_)
    : UserMatrix(shift_.dim()), A(A_), shift(shift_)
    { }

    /** Constructor for a dimension and a shift-value.
        So, this is the matrix diag(shift_).
        @param n The dimension.
        @param shift_ The shift value.
    */
    ShiftMatrix(int n, const dvector& shift_)
    : UserMatrix(n), A(0), shift(shift_)
    { }

    /** Computes the product of this matrix with a UserVector<double>.
        @param y The UserVector<double> to store the result in: y = A*x + diag(shift)*x.
        @param x The UserVector<double> to multiply with the matrix.
    */
    void MultV(UserVector<double>& y, const UserVector<double>& x) const {
      if (A) A->MultV(y,x);
      else y=0;
      y+=x.diagmult(shift);
		}

    /** Multiplies the matrix with a double*.
        @param y The double* to store the result in: y = A*x + diag(shift)*x.
        @param x The double* to multiply with this matrix.
        @see MultV(dvector&, const dvector&)
    */
    void MultV(double* y, const double* x) const {
      if (A) {
        A->MultV(y,x);
        for (int i=0; i<dim(); i++) y[i]+=shift(i)*x[i];
      }
      else for (int i=0; i<dim(); i++) y[i]=shift(i)*x[i];
    }

    /** Multiplies from left and from right with UserVector<double>.
        @param y The UserVector<double> to multiply with from left.
        @param x The UserVector<double> to multiply with from right.
        @return The double value y*A*x + y*diag(shift)*x.
    */
    double yAx(const UserVector<double>& y, const UserVector<double>& x) const {
      return (A ? A->yAx(y,x) : 0) + y * x.diagmult(shift);
    }

    /** Prints information about this matrix.
        Prints the shift-value and the wrapped matrix.
        @param out The ostream to print to.
    */
    void print(ostream &out) const {
      out << "ShiftMatrix: dim=" << dim() << " shift=" << shift << endl;
      if  (A) out << *A;
    }
};

/** Wrapper class to multiply a matrix with -1.
    *this == - *A.
*/
class MinusMatrix
#ifdef FILIB_AVAILABLE
: public IntervalCompliantMatrix
#else
: public UserMatrix
#endif
{
  protected:
    /** The matrix to multiply with -1.
    */
    Pointer<UserMatrix> A;

  public:
    /** Constructor for a UserMatrix.
        @param A_ The UserMatrix to wrap.
    */
    MinusMatrix(Pointer<UserMatrix> A_)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix
#else
		: UserMatrix
#endif
    (A_ ? A_->dim() : 0), A(A_)
    { assert(A_ != NULL);
    }

    /** Computes the product of this matrix with a dvector.
        @param y The dvector to store the result in: y = - A*x.
        @param x The dvector to multiply with the matrix.
    */
    void MultV(UserVector<double>& y, const UserVector<double>& x) const {
      A->MultV(y,x);
      y*=-1;
    }
		
    /** Multiplies the matrix with a double*.
        @param y The double* to store the result in.
        @param x The double* to multiply with this matrix.
        @see MultV(dvector&, const dvector&)
    */
    void MultV(double* y, const double* x) const {
      A->MultV(y,x);
      for (int i=0; i<dim(); i++) y[i]*=-1;
    }

    /** Multiplies from left and from right with UserVector<double>'s.
        @param y The UserVector<double> to multiply with from left.
        @param x The UserVector<double> to multiply with from right.
        @return The double value -y*this*x.
    */
    double yAx(const UserVector<double>& y, const UserVector<double>& x) const {
      return -A->yAx(y,x);
    }

#ifdef FILIB_AVAILABLE
		virtual void MultV(IntervalVector& y, const IntervalVector& x) const {
			assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A));
			((const IntervalCompliantMatrix*)(const UserMatrix*)A)->MultV(y,x);
			y*=-1.;
		}

		virtual interval<double> yAx(const IntervalVector& y, const IntervalVector& x) const {
			IntervalVector Ax(x.dim());
			assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A));
			((const IntervalCompliantMatrix*)(const UserMatrix*)A)->MultV(Ax,x);
			return -(y*Ax);
		}
#endif

    /** Prints information about this matrix.
        Prints the wrapped matrix.
        @param out The ostream to print to.
    */
    void print(ostream& out) const {
      out << "MinusMatrix: dim=" << dim() << endl << *A;
    }
};

/** A class to represent the sum of two matrices with optional scaling.
    this == a*A + b*B, where A and B are UserMatrices ( with ((A||B) != NULL) ) and a and b are double's.
*/
class SumMatrix
#ifdef FILIB_AVAILABLE
: public IntervalCompliantMatrix
#else
: public UserMatrix
#endif
{
  public:
    Pointer<const UserMatrix> A, B;
    double a, b;

    /** Constructor for two matrices and optional multipliers.
        A_ or B_ needs to be not NULL.
        @param A_ First matrix.
        @param B_ Second matrix.
        @param a_ First multiplier.
        @param b_ Second multiplier.
    */
    SumMatrix(Pointer<const UserMatrix> A_, Pointer<const UserMatrix> B_=NULL, double a_=1., double b_=1.)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix
#else
		: UserMatrix
#endif
    (A_ ? A_->dim() : (B_ ? B_->dim() : 0)), A(A_), B(B_), a(a_), b(b_)
    { assert(A_ || B_);
    }

    SumMatrix(Pointer<UserMatrix> A_, Pointer<UserMatrix> B_=NULL, double a_=1., double b_=1.)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix
#else
		: UserMatrix
#endif
    (A_ ? A_->dim() : (B_ ? B_->dim() : 0)), A(A_), B(B_), a(a_), b(b_)
    { assert(A_ || B_);
    }

    SumMatrix(UserMatrix* A_, UserMatrix* B_=NULL, double a_=1., double b_=1.)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix
#else
		: UserMatrix
#endif
    (A_ ? A_->dim() : (B_ ? B_->dim() : 0)), A(A_), B(B_), a(a_), b(b_)
    { assert(A_ || B_);
    }

		void MultV(dvector& y, const dvector& x) const {
			if (A) {
				A->MultV(y,x);
				y*=a;
				if (B) {
				  dvector y2(y.dim());
				  B->MultV(y2,x);
					y.AddMult(b, y2);
				}
			}	else { // no A, so there is a B
			  B->MultV(y,x);
			  y*=b;
			}
    }

    void MultV(UserVector<double>& y, const UserVector<double>& x) const {
      if (A) {
        A->MultV(y,x);
        y*=a;
        if (B) {
          Pointer<UserVector<double> > y2(y.getemptycopy());
          B->MultV(*y2, x);
          y.AddMult(b, *y2);
        }
      } else {
        B->MultV(y,x);
        y*=b;
      }
    }
		using UserMatrix::MultV;

#ifdef FILIB_AVAILABLE
		virtual void MultV(IntervalVector& y, const IntervalVector& x) const {
			if (A) assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A));
			if (B) assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)B));
			if (A) {
				((const IntervalCompliantMatrix*)(const UserMatrix*)A)->MultV(y,x);
				y*=a;
				if (B) {
				  IntervalVector y2(y.dim());
				  ((const IntervalCompliantMatrix*)(const UserMatrix*)B)->MultV(y2,x);
					y.AddMult(interval<double>(b), y2);
				}
			}	else { // no A, so there is a B
			  ((const IntervalCompliantMatrix*)(const UserMatrix*)B)->MultV(y,x);
			  y*=b;
			}
		}
#endif

    void print(ostream& out) const {
      out << "SumMatrix: dim=" << dim() << " a= " << a << " b= " << b << " A= ";
      if (A) out << endl << *A; else out << " NULL ";
      out << " B= ";
      if (B) out << endl << *B; else out << " NULL " << endl;
    }

};

/** A diagonal matrix.
*/
class DiagMatrix: public ExtUserMatrix {
   public:
      /** Indicates, whether the matrix is an unitary-matrix.
      */
      bool one;

      /** A Pointer to the diagonal elements.
      */
      Pointer<UserVector<double> > diag;

      /** Constructor for a Pointer to a UserVector<double>.
          Sets the Pointer to the diagonal elements to the given Pointer.
          @param b The diagonal elements.
      */
      DiagMatrix(Pointer<UserVector<double> > b_)
      : ExtUserMatrix(b_->dim()), diag(b_), one(*b_==1)
      { }

      /** Computes the diagonal-multiplication of this matrix with a UserVector<double>.
          @param y The UserVector<double> to store the result in: y=diag.diagmult(x)
          @param x The UserVector<double> to multiply with.
      */
      void MultV(UserVector<double>& y, const UserVector<double>& x) const {
        if (one) y=x;
        else diag->diagmult(y,x);
      }


#ifdef FILIB_AVAILABLE
			void MultV(IntervalVector& y, const IntervalVector& x) const {
				if (one) y=x;
				else for (int i=0; i<dim(); i++) (y[i]=x(i))*=(*diag)(i);
			}
#endif

#ifdef FILIB_AVAILABLE
			using IntervalCompliantMatrix::MultV;
#else
			using UserMatrix::MultV;
#endif

      /** Multiplies this matrix form left and right with UserVector<double>'s.
          @param y The UserVector<double> to multiply with from left.
          @param x The UserVector<double> to multiply with from right.
          @return The result x*this*y.
      */
      double yAx(const UserVector<double>& y,const UserVector<double>& x) const {
        if (one) return y*x;
        Pointer<UserVector<double> > z(y.getemptycopy());
        diag->diagmult(*z, x);
        return *z*y;
      }

#ifdef FILIB_AVAILABLE
			using IntervalCompliantMatrix::yAx;
#else
			using UserMatrix::yAx;
#endif
      
			double operator()(int row, int col) const {
        if (row==col)
          if (one) return 1.;
          else return (*diag)(row);
        else return 0;
      }

      /** Print's out information about this matrix.
          Print's the dimension and the diagonal elements.
          @param out The ostream to print to.
      */
      void print(ostream &out) const {
        out << "DiagMatrix: dim: " << dim() << " Diag: " << *diag;
      }
};

class BlockMatrix
#ifdef FILIB_AVAILABLE
: public IntervalCompliantMatrix
#else
: public UserMatrix
#endif
{
  public:
    /** The matrices of each block.
    */
    vector<Pointer<UserMatrix> > A;

    /** The block-structure.
    */
    vector<ivector> block;

		void set_dim() {
			dim_=0;
			for (int i=0; i<block.size(); i++) dim_+=block[i].size();
		}
		
    /** (Standard-)Constructor for the dimension.
        @param n The dimension, default is 0.
    */
    BlockMatrix(int n=0)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix(n)
#else
		: UserMatrix(n)
#endif
    { };


    /** Constructor for the block-structure.
        @param block_ The block-structure.
    */
    BlockMatrix(const vector<ivector>& block_)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix(),
#else
		: UserMatrix(),
#endif
    block(block_), A(block_.size())
    { set_dim(); }

    /** Constructor for the block-structure and matrix-structure.
        @param block_ The block-structure to copy.
        @param A_ The Pointer to the matrices of each block.
    */
    BlockMatrix(const vector<ivector>& block_, const vector<Pointer<UserMatrix> >& A_)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix(),
#else
		: UserMatrix(),
#endif
    block(block_), A(A_)
    { set_dim(); }

    /** Constructor for one block with a simple block-structure.
        @param A_ The pointer to the matrix for this block.
    */
    BlockMatrix(Pointer<UserMatrix> A_)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix(A_ ? A_->dim() : 0),
#else
		: UserMatrix(A_ ? A_->dim() : 0),
#endif
    block(1), A(1)
    { block[0].resize(dim_);
      for (int i=0; i<dim_; i++) block[0][i]=i;
      A[0]=A_;
    }

    /** Copy-Constructor.
        @param b The BlockMatrix to copy.
    */
    BlockMatrix(const BlockMatrix& b)
#ifdef FILIB_AVAILABLE
		: IntervalCompliantMatrix(b.dim()),
#else
		: UserMatrix(b.dim()),
#endif
    A(b.A), block(b.block)
    { }

    /** Multiplies with a UserVector<double>.
        @param y The UserVector<double> to store the result (this*x) in.
        @param x The UserVector<double> to multiply with.
    */
    void MultV(UserVector<double>& y, const UserVector<double>& x) const;
		using UserMatrix::MultV;

    /** Multiplies from left and from right with UserVector<double>'s.
        @param y The UserVector<double> to multiply with from left.
        @param x The UserVector<double> to multiply with from right.
        @return The double value y*this*x.
        @see UserMatrix::xAx(const UserVector<double>&)
    */
    double yAx(const UserVector<double>& y,const UserVector<double>& x) const;

#ifdef FILIB_AVAILABLE
		void MultV(IntervalVector& y, const IntervalVector& x) const;

		interval<double> yAx(const IntervalVector& y, const IntervalVector& x) const;
#endif

    /** Gives a pointer to the matrix of the k-th block.
        @param k The block-nr.
        @return The matrix for the k-th block.
    */
    Pointer<UserMatrix> operator[](int k) const {
      return A[k];
    }

    /** Prints out some information about this matrix.
        Prints the dimension, the block-structure and the matrix of each defined block.
        @param out The ostream to print to.
    */
    void print(ostream &out) const {
      out << "BlockMatrix: dim: " << dim() << " Blocks: " << block.size() << endl;
      for (int i=0; i<block.size(); i++) {
        out << "block " << i << ": " << block[i];
        if (A[i]) out << *A[i];
      }
    }

};

//--------------------------------------------------------

/** A wrapper class for a dense TNT-matrix.
*/
class DenseMatrix : public ExtUserMatrix {
  private:
    /** Reduce a real symmetric matrix to a symmetric tridiagonal matrix.
        This subroutine is a translation of the algol procedure tred2,
        num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
        handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

        This subroutine reduces a real symmetric matrix to a
        symmetric tridiagonal matrix using and accumulating
        orthogonal similarity transformations.

        This version dated august 1983.

        a and z may coincide.  If distinct, a is unaltered.
        @param nm The row dimension of two-dimensional array parameters as declared in the calling program dimension statement.
        @param n The order of the matrix.
        @param a The real symmetric input matrix. Only the lower triangle of the matrix need be supplied.
        @param d At end contains the diagonal elements of the tridiagonal matrix.
        @param e At end contains the subdiagonal elements of the tridiagonal matrix in its last n-1 positions. e(1) is set to zero.
        @param z At end contains the orthogonal transformation matrix produced in the reduction.
    */
    void tred2(int nm, int n, double* a, double *d, double *e, double *z) const;

    /** Finds the eigenvalues and eigenvectors of a symmetric tridiagonal matrix by the implicit ql-method.
        This subroutine finds the eigenvalues and eigenvectors of a symmetric tridiagonal matrix by the implicit ql method.
        The eigenvectors of a full symmetric matrix can also be found if  tred2  has been used to reduce this full matrix to tridiagonal form.

        This subroutine is a translation of the algol procedure imtql2, num. math. 12, 377-383(1968) by martin and wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
        handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).

        calls pythag for dsqrt(a*a + b*b) .

        this version dated august 1983.

        @param nm must be set to the row dimension of two-dimensional array parameters as declared in the calling program dimension statement.
        @param n The order of the matrix.
        @param d The diagonal elements of the input matrix.
                 On output: The eigenvalues in ascending order. If an error exit is made, the eigenvalues are correct but unordered for indices 1,2,...,ierr-1.
        @param e The subdiagonal elements of the input matrix in its last n-1 positions. e(1) is arbitrary.
                 On output: Destroyed.
        @param z The transformation matrix produced in the reduction by tred2, if performed.
                 If the eigenvectors of the tridiagonal matrix are desired, z must contain the identity matrix.
                 On output: The orthonormal eigenvectors of the symmetric tridiagonal (or full) matrix. if an error exit is made, z contains the eigenvectors associated with the stored eigenvalues.

        @return zero for normal return, j if the j-th eigenvalue has not been determined after 30 iterations.
    */
    int imtql2(int nm, int n, double *d, double *e, double *z) const;

  protected:
    /** The matrix to wrap.
    */
    TNT::Matrix<double> A;

  public:
    /** Indicates, whether it is allowed to destroy the matrix during the eigenvalue computation.
        True save's memory and time.
    */
    bool allow_destroy;

    /** Standard-Constructor.
        @param allow_destroy_ Indicates, whether it is allowed to destroy the matrix during the eigenvalue computation, default is false.
        @see DenseMatrix(int, double, bool)
        @see DenseMatrix(UserMatrix&, bool)
        @see DenseMatrix(TNT::Matrix<double>&, bool)
    */
    DenseMatrix(bool allow_destroy_=false)
    : ExtUserMatrix(), A(), allow_destroy(allow_destroy_)
    { }

    /** Constructor for given dimension and optional initial argument.
        @param n The dimension.
        @param val A value for all elements of the matrix, default is 0.
        @param allow_destroy_ Indicates, whether it is allowed to destroy the matrix during the eigenvalue computation.
        @see DenseMatrix(bool)
        @see DenseMatrix(UserMatrix&, bool)
        @see DenseMatrix(TNT::Matrix<double>&, bool)
    */
    DenseMatrix(int n, double val=0., bool allow_destroy_=false)
    : ExtUserMatrix(n), A(n, n, val), allow_destroy(allow_destroy_)
    { }

    /** Copy-Constructor for a UserMatrix.
        Creates a dense matrix from a UserMatrix.
        It gets the collumns by calling MultV with standard-basis-vectors.
        @param A_ The UserMatrix to copy.
        @param allow_destroy_ Indicates, whether it is allowed to destroy the matrix during the eigenvalue computation.
        @see DenseMatrix(bool)
        @see DenseMatrix(int, double, bool)
        @see DenseMatrix(TNT::Matrix<double>&, bool)
    */
    DenseMatrix(const UserMatrix& A_, bool allow_destroy_=false);

    /** Copy-Constructor for a ExtUserMatrix.
        @param A_ The matrix to copy.
        @param allow_destroy_ Indicates, whether it is allowed to destroy the matrix during the eigenvalue computation.
    */
    DenseMatrix(const ExtUserMatrix& A_, bool allow_destroy_=false)
    : ExtUserMatrix(A_.dim()), A(A_.dim(), A_.dim(), 0.), allow_destroy(allow_destroy_)
    { for (int i=0; i<dim(); i++)
    		for (int j=0; j<dim(); j++)
    		  A[i][j]=A_(i,j);
    }

    /** Copy-Constructor for a TNT::Matrix<double>.
        @param A_ The matrix to copy.
        @param allow_destroy_ Indicates, whether it is allowed to destroy the matrix during the eigenvalue computation.
        @see DenseMatrix(bool)
        @see DenseMatrix(int, double, bool)
        @see DenseMatrix(UserMatrix&, bool)
    */
    DenseMatrix(const TNT::Matrix<double>& A_, bool allow_destroy_=false)
    : ExtUserMatrix(A_.num_rows()), A(A_), allow_destroy(allow_destroy_)
    { }

    /** Copy-Constructor.
        @param D The DenseMatrix to copy.
    */
    DenseMatrix(const DenseMatrix& D)
    : ExtUserMatrix(D.dim()), A(D.A), allow_destroy(D.allow_destroy)
    { }

    /** Gives the elements of this matrix as double*.
        Rowwise.
    */
    operator double*() {
      return *(double**)A;
    }

    /** Assign-Operator for a DenseMatrix.
        @param D The DenseMatrix to copy.
        @see operator=(const double)
        @return This matrix.
    */
    DenseMatrix& operator=(const DenseMatrix& D) {
      if (this != &D) {
        dim_=D.dim();
        A=D.A;
        allow_destroy=D.allow_destroy;
      }
      return *this;
    }

    /** Assign-Operator for a double.
        @param scalar A double to set all elements of the matrix to.
        @see operator=(const DenseMatrix&)
    */
    DenseMatrix& operator=(const double scalar) {
      A=scalar;
      return *this;
    }

    /** Returns one row of this matrix.
        @param i The index of the row.
        @return The i'th row as dvector.
        @see set_row(dvector&, const int)
        @see operator()(int, int)
    */
    dvector operator[](int i) const {
      return dvector(A[i], dim());
    }

    /** Set's a dvector to one row of this matrix.
        @param i The index of the row to get.
        @see operator[](int)
    */
/*    void set_row(dvector &a, const int i) {
      a=operator[](i);
    }
*/
    /** Return one element of the matrix.
        @param row The row.
        @param col The column.
        @return The element at index [row, col] (C-style).
        @see operator[](int)
    */
    double& operator()(int row, int col) {
      return A[row][col];
    }

    double operator()(int row, int col) const {
      return A[row][col];
    }

    /** Building the sum of this matrix and another one.
        @param B Second summand.
        @return The sum *this + B as DenseMatrix.
        @see operator-(const DenseMatrix&)
    */
    DenseMatrix operator+(const DenseMatrix &B) const {
      return DenseMatrix(A + B.A);
    }

    /** Building the difference of this matrix and another one.
        @param B Matrix to substract.
        @return The difference *this - B as DenseMatrix.
        @see operator+(const DenseMatrix&)
    */
    DenseMatrix operator-(const DenseMatrix &B) const {
      return DenseMatrix(A - B.A);
    }

    /** Building the product of this matrix and another one.
        @param B Matrix to multiply with.
        @return The product *this * B as DenseMatrix.
        @see UserMatrix::operator*(const dvector&)
        @see UserMatrix::operator*(const double*)
    */
    DenseMatrix operator*(const DenseMatrix &B) const {
      return DenseMatrix(A * B.A);
    }

    /** Multiplication with a UserVector<double>.
        Calls the *-operator form the TNT-matrix.
        @param y_ The UserVector<double> to store the result in.
        @param x_ The UserVector<double> to multiply with.
        @see MultV(double*, const double*)
    */
    void MultV(UserVector<double>& y_, const UserVector<double>& x_) const;

#ifdef FILIB_AVAILABLE
		void MultV(IntervalVector& y, const IntervalVector& x) const;
#endif

    /** Multiplication with a double*.
        Calls the *-operator from the TNT-matrix.
        @param y_ The double* to store the result in.
        @param x_ The double* to multiply with.
        @see MultV(dvector&, const dvector&)
    */
    void MultV(double* y_, const double* x_) const;

    void set_random(const dvector& lambda);

    /** Computes the eigenvalues and eigenvectors with the implicit ql-method.
        If allow_destroy is true, the matrix itselfe is destroyed and holds the computed eigenvectors (in transposed form) after exit.
        @param eig_vec A vector of dvectors for the eigenvectors.
        @param eig_val A vector of double for the eigenvalues.
        @return The return value from imtql2, which is 0, if all went good or the number of the eigenvalue, which wasn't founds after 30 iterations.
        @see tred2(int, int, double*, double*, double*, double*)
        @see imtql2(int, int, double*, double*, double*)
        @see eig_ql(dvector&, double&)
    */
    int eig_ql(vector<dvector> &eig_vec, vector<double> &eig_val) const;

    /** Computes the minimum eigenvalues and eigenvectors with the implicit ql-method.
        Computes all eigenvalues and look's than for the smallest one.
        @param eig_vec A dvectors for the eigenvector.
        @param eig_val A double for the eigenvalue.
        @return The return value from imtql2, which is 0, if all went good or the number of the eigenvalue, which wasn't founds after 30 iterations.
        @see eig_ql(vector<dvector>&, vector<double>&)
    */
    int eig_ql(dvector &eig_vec, double &eig_val) const;

    /** Print some information.
        Calls the output-operator from the TNT-matrix.
        @param out The ostream to print to.
    */
    void print(ostream &out) const {
      out << "DenseMatrix: " << A;
    }
};

class SparseMatrix {
  protected:
    /** A map to construct the matrix.
        Maps (col, row) to values.
    */
    map<pair<int,int>, double> values;
	
    /** Number of nonzero-elements.
        If the matrix is not defined, this is -1.
    */
  	int nz;
  	/** The values of the matrix.
  	    size=nz.
  	*/
  	double* val;
  	/** The row indices of the values.
  	    size=nz.
  	*/
  	int* row_ind;
  	/** The indices of the starting columns.
  	    size=dim()+1.
  	*/
  	int* col_ptr;

		/** The number of rows and columns.
		*/
		int rows_, cols_;	

  public:
    /** Constructor for a given dimension.
        @param rows__ The number of rows_.
		    @param cols__ The number of columns.
    */
    SparseMatrix(int rows__, int cols__)
    : rows_(rows__), cols_(cols__), val(0), row_ind(0), col_ptr(0), nz(-1)
    { }

    /** Copy-Constructor for a SparseMatrix.
        @param A_ The SparseMatrix to copy.
    */
    SparseMatrix(const SparseMatrix& A_)
		: rows_(A_.rows_), cols_(A_.cols_), nz(A_.nz),
		  val(A_.nz>=0 ? new double[A_.nz] : NULL),
		  row_ind(A_.nz>=0 ? new int[A_.nz] : NULL),
		  col_ptr(A_.nz>=0 ? new int[A_.cols_+1] : NULL)
    { if (nz>0) {
    		memcpy(val, A_.val, nz * sizeof(double));
	      memcpy(row_ind, A_.row_ind, nz * sizeof(double));
  	    memcpy(col_ptr, A_.col_ptr, (cols_+1) * sizeof(double));
  	  }
    }

		/** Copy-Constructor for a UserMatrix.
		    @param A_ The UserMatrix to copy.
        @param no_finish Indicates, whether this SparseMatrix shouldn't be finished after copying A_.
		*/
		SparseMatrix(const UserMatrix& A_, bool no_finish=false);

		/** Copy-Constructor for an ExtUserMatrix.
		    @param A_ The UserMatrix to copy.
        @param no_finish Indicates, whether this SparseMatrix shouldn't be finished after copying A_.
		*/
		SparseMatrix(const ExtUserMatrix& A_, bool no_finish=false)
		: rows_(A_.dim()), cols_(A_.dim()), val(NULL), row_ind(NULL), col_ptr(NULL), nz(-1)
		{ for (int i=0; i<rows_; i++)
				for (int j=0; j<cols_; j++)
					AddElement(i,j,A_(i,j));
			if (!no_finish) finish();
		}
		
		/** Destructor.
		    Deletes val, row_ind and col_ptr, if not NULL.
		*/
		virtual ~SparseMatrix() {
		  if (val) delete[] val;
		  if (row_ind) delete[] row_ind;
		  if (col_ptr) delete[] col_ptr;
		}

		/** Resizes the matrix, if it is not finished already.
				Elements with indices outside the dimension are removed.
		*/
		void resize(int rows__, int cols__);

		int rows() const { return rows_; }
		int cols() const { return cols_; }
		
		const int* GetRowInd() const { return row_ind; }
		const int* GetColPtr() const { return col_ptr; }
		const double* GetVal() const { return val; }
		double* GetVal() { return val; }
		
		virtual double operator()(int row, int col) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
			assert(val!=NULL);
#endif
		  for (int i=col_ptr[col]; i<col_ptr[col+1]; i++)
		    if (row_ind[i]==row) return val[i];
		  return 0.;
		}
		
		int nonzeros() const {
			if (nz>=0) return nz;
			else return values.size();
		}

    /** Adds an element to the map.
        @param row The row.
        @param col The column.
        @param v The value.
		    @param check_zero If set to true (default), and |v|<rtol, the element is not added.
        @see values
    */
    void AddElement(int row, int col, double v, bool check_zero=true) {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(val==NULL);
      assert(0<=row && row<rows_);
      assert(0<=col && col<cols_);
#endif
      if (check_zero && fabs(v)<rtol) return;
      values.insert(pair<pair<int,int>, double>(pair<int,int>(col, row), v));
    }

    void AddToElement(int row, int col, double v, bool check_zero=true) {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(val==NULL);
      assert(0<=row && row<rows_);
      assert(0<=col && col<cols_);
#endif
			if (check_zero && fabs(v)<rtol) return;
			map<pair<int,int>, double >::iterator it(values.find(pair<int,int>(col, row)));

			if (it==values.end()) values.insert(pair<pair<int,int>, double>(pair<int,int>(col, row), v));
			else it->second+=v;
    }

    /** Adds a UserMatrix to this SparseMatrix.
        Computes the elements by multiplication with unit vectors.
        This matrix needs to be unfinished to do this.
        @param A_ The UserMatrix to add.
        @return This matrix.
    */
    SparseMatrix& operator+=(const UserMatrix& A_);

    /** Adds a ExtUserMatrix to this SparseMatrix2.
        This matrix needs to be unfinished to do this.
        @param A_ The ExtUserMatrix to add.
        @return This matrix.
    */
    SparseMatrix& operator+=(const ExtUserMatrix& A_) {
#ifndef NO_SPARSEMATRIX_ASSERTS
			assert(rows_==A_.dim());
			assert(cols_==A_.dim());
#endif
    	for (int row=0; row<rows_; row++)
    		for (int col=0; col<cols_; col++)
    			AddToElement(row, col, A_(row,col));
			return *this;
    }

    SparseMatrix& operator=(const double v);

    SparseMatrix& operator*=(const double v);

		void set_block(const SparseMatrix& A, const ivector& indices);

    /** Finish the matrix.
        Constructs the arrays val, row_ind and col_ptr, using the map.
        Clears values.
    */
    void finish();

		/** Multiplication with a double*.
		    @param y The double* to store the result in.
		    @param x The double* to multipliy with.
		*/
		virtual void MultV(double* y, const double* x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(val);
#endif
		  memset(y, 0, rows_ * sizeof(double));
		  int j=0;
		  for (int col=0; col<cols_; col++, x++)
		    for (; j<col_ptr[col+1]; j++)
		      y[row_ind[j]] += *x * val[j];
		}

		/** Multiplication with a dvector.
		    @param y The dvector to store the result in.
		    @param x The dvector to multipliy with.
		*/
    virtual void MultV(dvector& y, const dvector& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(y.dim()==rows_);
      assert(x.dim()==cols_);
#endif
      MultV((Pointer<double>)y, (Pointer<double>)x);
    }

		/** Multiplication with a UserVector<double>.
		    @param y The UserVector<double> to store the result in.
		    @param x The UserVector<double> to multipliy with.
		*/
    virtual void MultV(UserVector<double>& y, const UserVector<double>& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(val);
			assert(y.dim()==rows_);
			assert(x.dim()==cols_);
#endif
      y=0;
		  int j=0;
      double xi;
		  for (int col=0; col<cols_; col++)
        if (xi=x(col))
  		    for (j=col_ptr[col]; j<col_ptr[col+1]; j++)
	  	      y[row_ind[j]] += xi * val[j];
		}

#ifdef FILIB_AVAILABLE
		virtual void MultV(IntervalVector& y, const IntervalVector& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(val);
			assert(y.dim()==rows_);
			assert(x.dim()==cols_);
#endif
		  int j=0;
			interval<double> zero(0.);
			y=zero;
		  for (int col=0; col<cols_; col++)
        if (x(col)!=zero)
  		    for (j=col_ptr[col]; j<col_ptr[col+1]; j++)
	  	      y[row_ind[j]] += val[j] * x(col);
		}

		virtual interval<double> yAx(const IntervalVector& y, const IntervalVector& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(val);
			assert(y.dim()==rows_);
			assert(x.dim()==cols_);
#endif
      interval<double> ret(0.);
			interval<double> zero(0.);

		  int j=0;
		  for (int col=0; col<cols_; col++)
        if (x(col)!=zero)
  		    for (j=col_ptr[col]; j<col_ptr[col+1]; j++)
	  	      ret+=y(row_ind[j]) * x(col) * val[j];

      return ret;
    }
    
    virtual interval<double> xAx(const IntervalVector& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
		assert(val);
		assert(rows_==cols_);
		assert(x.dim()==cols_);
#endif
		interval<double> ret(0.);
		interval<double> zero(0.);

		int j=0;
		for (int col=0; col<cols_; col++)
		if (x(col)!=zero) {
			for (j=col_ptr[col]; j<col_ptr[col+1]; j++)
				if (row_ind[j]==col) ret+=val[j]*sqr(x(col));
				else ret+=x(row_ind[j]) * x(col) * val[j];
		}
		return ret;
    }
#endif

		/** Multiplication with a SparseVector<double>.
		    @param y The SparseVector<double> to store the result in.
		    @param x The SparseVector<double> to multipliy with.
		*/
		virtual void MultV(SparseVector<double>& y, const SparseVector<double>& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
		  assert(val);
			assert(y.dim()==rows_);
			assert(x.dim()==cols_);
#endif
		  double* y0=new double[cols_]; memset(y0, 0, cols_*sizeof(double));
		  int j;
      SparseVector<double>::VectorElement* v=x.head->next;
      while (v) {
		    for (j=col_ptr[v->index]; j<col_ptr[v->index+1]; j++)
		      y0[row_ind[j]] += v->value * val[j];
		    v=v->next;
		  }
		  y.set(y0);
		  delete y0;
		}

    virtual double yAx(const UserVector<double>& y, const UserVector<double>& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
      assert(val);
			assert(y.dim()==rows_);
			assert(x.dim()==cols_);
#endif
      double ret=0;

		  int j=0;
      double xi;
		  for (int col=0; col<cols_; col++)
        if (xi=x(col))
  		    for (j=col_ptr[col]; j<col_ptr[col+1]; j++)
	  	      ret+=y(row_ind[j]) * xi * val[j];

      return ret;
    }

    virtual double yAx(const dvector& y, const dvector& x) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
			assert(y.dim()==rows_);
			assert(x.dim()==cols_);
      assert(val);
#endif
      double ret=0;
      double* y0=(Pointer<double>)y;
      double* x0=(Pointer<double>)x;

      int j=0;
      for (int col=0; col<cols_; col++, x0++)
        for (; j<col_ptr[col+1]; j++)
          ret+=y[row_ind[j]] * *x0 * val[j];

      return ret;
    }

    virtual double xAx(const UserVector<double>& x) const { return yAx(x,x); }

    /** Adds to a UserVector<double> the product of a double value and this matrix, multiplied with a UserVector<double>.
        @param y The UserVector<double> to store the result in: y + alpha * A * x
        @param x The UserVector<double> to multiply this matrix with.
        @param alpha The double to multiply with.
    */
    virtual void AddMult(UserVector<double>& y, const UserVector<double>& x, const double alpha) const {
#ifndef NO_SPARSEMATRIX_ASSERTS
			assert(y.dim()==rows_);
			assert(x.dim()==cols_);
      assert(val);
#endif
		  int j=0;
      double xi;
		  for (int col=0; col<cols_; col++)
        if (xi=alpha*x(col))
		      for (j=col_ptr[col]; j<col_ptr[col+1]; j++)
		        y[row_ind[j]] += xi * val[j];
    };

    /** Plots the edges of the sparsity pattern to a file, which can be read by gnuplot.
        For each entry (i,j) in the matrix, it prints the i and j in one line.
        @param filename The name of the file to print to.
    */
    void plot(char* filename) const;

    /** Prints the matrix.
        Prints the row, column and value for each non-zero entry in this SparseMatrix2.
    */
    virtual void print(ostream& out) const;

};

class SparseMatrix2 : public ExtUserMatrix, public SparseMatrix {
	public:
    /** Constructor for a given dimension.
        @param n The dimension.
    */
    SparseMatrix2(int n)
    : SparseMatrix(n, n), ExtUserMatrix(n)
    { }

    /** Copy-Constructor for a SparseMatrix.
        @param A_ The quadratic(!) SparseMatrix to copy.
    */
    SparseMatrix2(const SparseMatrix& A_)
		: SparseMatrix(A_), ExtUserMatrix(A_.rows())
    { assert(A_.rows()==A_.cols());
		}

		/** Copy-Constructor for a UserMatrix.
		    @param A_ The UserMatrix to copy.
        @param no_finish Indicates, whether this SparseMatrix shouldn't be finished after copying A_.
		*/
		SparseMatrix2(const UserMatrix& A_, bool no_finish=false)
		: SparseMatrix(A_, no_finish), ExtUserMatrix(A_.dim())
		{ }

		/** Copy-Constructor for an ExtUserMatrix.
		    @param A_ The UserMatrix to copy.
        @param no_finish Indicates, whether this SparseMatrix shouldn't be finished after copying A_.
		*/
		SparseMatrix2(const ExtUserMatrix& A_, bool no_finish=false)
		: SparseMatrix(A_, no_finish), ExtUserMatrix(A_.dim())
		{ }
		
		double operator()(int row, int col) const {	return SparseMatrix::operator()(row, col);	}

    SparseMatrix2& operator+=(const UserMatrix& A_) {	SparseMatrix::operator+=(A_);	return *this;	}

    SparseMatrix2& operator+=(const ExtUserMatrix& A_) { SparseMatrix::operator+=(A_); return *this; }

    SparseMatrix2& operator=(const double v) { SparseMatrix::operator=(v); return *this; }

    SparseMatrix2& operator*=(const double v) {	SparseMatrix::operator*=(v); return *this; }

		void set_block(const SparseMatrix& A, const ivector& indices) { SparseMatrix::set_block(A, indices); }

		void MultV(double* y, const double* x) const { SparseMatrix::MultV(y,x); }

    void MultV(dvector& y, const dvector& x) const { SparseMatrix::MultV(y,x); }

    void MultV(UserVector<double>& y, const UserVector<double>& x) const { SparseMatrix::MultV(y,x); }

#ifdef FILIB_AVAILABLE
    void MultV(IntervalVector& y, const IntervalVector& x) const { SparseMatrix::MultV(y,x); }
#endif

		void MultV(SparseVector<double>& y, const SparseVector<double>& x) const { SparseMatrix::MultV(y,x); }

    double yAx(const UserVector<double>& y, const UserVector<double>& x) const { return SparseMatrix::yAx(y,x); }

    double yAx(const dvector& y, const dvector& x) const { return SparseMatrix::yAx(y,x); }

#ifdef FILIB_AVAILABLE
		using IntervalCompliantMatrix::yAx;
#else		
		using UserMatrix::yAx;
#endif
    
		double xAx(const UserVector<double>& x) const { return SparseMatrix::xAx(x); }

#ifdef FILIB_AVAILABLE
		interval<double> xAx(const IntervalVector& x) const { return SparseMatrix::xAx(x); }
#endif

    virtual void AddMult(UserVector<double>& y, const UserVector<double>& x, const double alpha) const { SparseMatrix::AddMult(y,x,alpha); }
		
		void make_symmetric();
		
    /** Clears the matrix and sets the elements to random elements.
        First, the diagonal elements are set.
        @param num_el The number of elements in one half of the matrix.
        @param max The maximum absolute value of the elements.
    */
    void set_random(int num_el, double max=1.);
		
    /** Gives a symmetric random SparseMatrix2.
        First, the diagonal elements are set.
        @param n The dimension.
        @param num_el The number of elements in one half of the matrix.
        @param max The maximum absolute value of the elements.
        @return A new SparseMatrix2.
    */
    static SparseMatrix2* random(int n, int num_el, double max=1.) {
		  SparseMatrix2* A=new SparseMatrix2(n);
		  A->set_random(num_el, max);
		  return A;
		}
		
    virtual void print(ostream& out) const { SparseMatrix::print(out); }
};

#endif // USERMATRIX_H
