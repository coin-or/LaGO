// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef USERVECTOR_H
#define USERVECTOR_H

// TNT-stuff
#ifndef TNT_NO_BOUNDS_CHECK
#define TNT_BOUNDS_CHECK
#endif
#include <tnt_vec.h>

template <class Type> class DenseVector;

/** Abstract class for a UserVector.
    You need to implement the following method:
    - UserVector<Type>* getemptycopy(int=dim()) const
    - int dim() const
    - Type& operator[](int)
    - Type& operator()(int) const
    - UserVector<Type>& operator=(const UserVector<Type>&)
    - UserVector<Type>& operator=(const Type&)
    - UserVector<Type>& operator+=(const UserVector<Type>&)
    - UserVector<Type>& operator-=(const UserVector<Type>&)
    - UserVector<Type>& operator*=(const Type&)
    - UserVector<Type>& operator/=(const Type&)
    - operator const Pointer<Type>() const
*/
template <class Type> class UserVector {
  /** Outputs the components of this UserVector.
      Calls a.print(out).
      @param out An output-stream to print to.
      @param v The UserVector to print.
      @return The ostream out.
      @see print(ostream&)
  */
	friend ostream& operator << (ostream& out, const UserVector<Type>& v) {
		v.print(out);
		return out;
	}
  
	friend void operator>>(const UserVector<Type>& v, Type* out) {
		for (int i=0; i<v.dim(); ++i) out[i]=v(i);
	}

 	friend istream& operator >> (istream& in, UserVector<Type>& v) {
		for (int i=0; i<v.dim(); i++) in >> v[i];
		return in;
	}

  public:
    /** Standard-Constructor.
        @see dim()
    */
    UserVector() { }

    /** Virtual Destructor.
    */
    virtual ~UserVector() { }

    /** Gives an empty UserVector of this type.
        Abstract.
        @param n The size of the new UserVector.
        @return A Pointer to a new UserVector with size n.
    */
    virtual Pointer<UserVector<Type> > getemptycopy(int n) const=0;

    /** Gives an empty UserVector of this type and dimension.
        @return A Pointer to a new UserVector with the size of this one.
        @see getemptycopy(int)
    */
    virtual Pointer<UserVector<Type> > getemptycopy() const {
      return getemptycopy(dim());
    }

    /** Gives a copy of this UserVector.
        @return A Pointer to a new UserVector, which is a copy of this one.
    */
    virtual Pointer<UserVector<Type> > getcopy() const {
      Pointer<UserVector<Type> > ret(getemptycopy());
      *ret=*this;
      return ret;
    }

    /** Gives a Pointer to a copy of one block of this UserVector.
        @param block The block-describtion.
        @return A Pointer to a copy of the block, described in block.
    */
    virtual Pointer<UserVector<Type> > getcopy(const UserVector<int>& block) const {
      Pointer<UserVector<Type> > ret(getemptycopy(block.size()));
			for (int i=0; i<block.size(); i++) ret->SetElement(i, (*this)(block(i)));
      return ret;
    }

    /** The dimension.
        Abstract.
        @return The dimension of this UserVector.
    */
    virtual int dim() const=0;

    /** The size.
        @return The dimension.
        @see dim()
    */
    int size() const { return dim(); }

    /** Tests the new size.
        @param n The new size.
    */
    virtual void resize(const int n) {
      assert(n>=0);
    }

    /** Gives a reference to one component of this UserVector.
        Abstract.
        Use operator(int), when you want to read a component only!
        @param i The index of the element.
        @return A reference to the component at the i-th index.
        @see operator[](const int)
    */
    virtual Type& operator[](const int i)=0;

    /** Gives one element of this UserVector.
        @param i The index of the element.
        @return The element at index i, not a reference to it.
        @see operator[](int)
    */
    virtual Type operator[](const int i) const {
      return (*this)(i);
    }

    /** Gives a copy of one component of this UserVector.
        Use this one, when you want to read a component only!
        Abstract
        @param i The index of the element.
        @return A copy of the the component at the i-th index.
        @see operator[](const int)
    */
    virtual Type operator()(const int i) const=0;

    /** Gives a block form this UserVector as DenseVector.
        @param block The block structure.
        @return The block of this UserVector, which was described in block.
        @see DenseVector<Type>::DenseVector(const UserVector<Type>&, UserVector<int>&)
    */
    virtual DenseVector<Type> operator()(const UserVector<int>& block) const {
      return DenseVector<Type>(*this, block);
    }

    /** Gives one part of this UserVector as DenseVector.
        @param low The lower index.
        @param up The upper index.
        @return The elements low to up as new UserVector.
        @see DenseVector<Type>::DenseVector(const UserVector<Type>&, const int, const up)
    */
    virtual DenseVector<Type> operator()(const int low, const int up) const {
      return DenseVector<Type>(*this, low, up);
    }

    /** Sets one element of this UserVector.
        @param i The index of the element.
        @param v The value to set.
        @see operator[](const int i)
    */
    virtual void SetElement(const int i, const Type& v) {
      (*this)[i]=v;
    }

    /** Sets one block of this UserVector.
        @param v The UserVector to set.
        @param block The block-describtion.
    */
    virtual void set_block(const UserVector<Type>& v, const UserVector<int>& block) {
      for (int i=0; i<block.size(); i++) (*this)[block[i]]=v(i);
    }

    /** Cast-Operator to a const Pointer of Type.
    */
    virtual operator const Pointer<Type>() const=0;

    /** Assign-Operator for a UserVector.
        Abstract.
        @param v The UserVector to assign to this UserVector.
        @return This UserVector.
    */
    virtual UserVector<Type>& operator=(const UserVector<Type>& v)=0;

    /** Assign-Operator for a value from type Type.
        Abstract.
        @param v The value to assign to all elements of this vector.
        @return This UserVector.
    */
    virtual UserVector<Type>& operator=(const Type& v)=0;

    /** Assign-Operator for a Type*.
        @param v The Type*.
        @param n The length of v.
        @return This UserVector.
    */
    virtual void set(const Type* v, int n) {
      for (int i=0; i<n; i++) SetElement(i, v[i]);
    }

    /** Sets this UserVector to the values of v.
        Be careful with this one. v needs to have the length dim().
        @param v The Type*.
        @return This UserVector.
    */
    void set(const Type* v) { set(v, dim()); }

    /** Adds a UserVector, multiplied by a Type, to this UserVector.
        @param a The Type to multiply v with.
        @param v The UserVector to add, after multiplication with a.
        @return This UserVector
    */
    virtual UserVector<Type>& AddMult(Type a, const UserVector<Type>& v) {
      return (*this)+=a*v;
    }

    /** Adds a UserVector to this one.
        Abstract.
        @param v The UserVector to add.
        @return This UserVector.
    */
    virtual UserVector<Type>& operator+=(const UserVector<Type>& v)=0;

    /** Substracts a UserVector from this one.
        Abstract.
        @param v The UserVector to substract.
        @return This UserVector.
    */
    virtual UserVector<Type>& operator-=(const UserVector<Type>& v)=0;

    /** Multiplies all components with a value from type Type.
        Abstract.
        @param v The value to multiply with.
        @return This UserVector.
    */
    virtual UserVector<Type>& operator*=(const Type& v)=0;

    /** Divides all components through a value from type Type.
        Abstract.
        @param v The value to divide.
        @return This UserVector.
    */
    virtual UserVector<Type>& operator/=(const Type& v)=0;

    /** Compares each element of this UserVector with a value of type Type.
        @param v The value to compare with.
        @return True, if all elements are greater-equal v. False else.
    */
    virtual bool operator>=(const Type& v) const {
      for (int i=0; i<dim(); i++) if ((*this)(i)<v) return false;
      return true;
    }

    /** Compares each element of this UserVector with a value of type Type.
        @param v The value to compare with.
        @return True, if all elements are equal to v. False else.
    */
    virtual bool operator==(const Type& v) const {
      for (int i=0; i<dim(); i++) if ((*this)(i)!=v) return false;
      return true;
    }

    /** Compares this UserVector with another one.
        @param v The UserVector to compare with.
        @return True, if this UserVectors is equal to v. False else.
    */
    virtual bool operator==(const UserVector<Type>& v) const {
      assert(v.dim()==dim());
      for (int i=0; i<dim(); i++) if ((*this)(i)!=v(i)) return false;
      return true;
    }

    /** Computes the scalar product of this UserVector and another one.
        @param v The UserVector to multiply with.
        @return The scalar product *this * v.
    */
		virtual Type operator*(const UserVector<Type>& v) const {
      assert(v.dim()==dim());
      Type val=Type();
      for (int i=0; i<dim(); i++) val+=(*this)(i)*v(i);
      return val;
    }

		/** Computes the product of this UserVector and a value from type v.
        @param v The value to multiply with.
        @return A new DenseVector, which is this one multiplied by v.
    */
    virtual DenseVector<Type> operator*(const Type& v) const {
      return DenseVector<Type>(*this)*v;
    }

    /** Computes the product of a value from type Type and a UserVector.
        @param v1 The value of type Type.
        @param v2 The UserVector to multiply with.
        @return A new DenseVector, which is the product of v1 and v2.
        @see operator*(const Type&)
    */
    friend DenseVector<Type> operator*(const Type& v1, const UserVector<Type>& v2) {
      return v2 * v1;
    }

    /** Computes the sum of this UserVector and another UserVector.
        @param v The UserVector to add.
        @return A new DenseVector, which is the sum of this one and v.
        @see operator+=(const UserVector<Type>&)
    */
    virtual DenseVector<Type> operator+(const UserVector<Type>& v) const {
      return DenseVector<Type>(*this) + v;
    }

    /** Computes nothing.
        @return This UserVector as const.
    */
    virtual const UserVector<Type>& operator+() const {
      return *(const UserVector<Type>*)this;
    }

    /** Computes nothing.
        @return This UserVector.
    */
    virtual UserVector<Type>& operator+() {
      return *this;
    }

    /** Computes the difference between this UserVector and another one.
        @param v The UserVector to substract.
        @return The difference *this - v as a new DenseVector.
        @see operator-=(UserVector<Type>&)
    */
    virtual DenseVector<Type> operator-(const UserVector<Type>& v) const {
      return DenseVector<Type>(*this)-v;
    }

    /** Calculates this UserVector multiplied with -1.
        Abstract.
        @return A new DenseVector, which is this one, multiplied by -1.
    */
    virtual DenseVector<Type> operator-() const {
      return -DenseVector<Type>(*this);
    }

    /** Calculates the diagonal-multiplication.
        Multiplies each component of this UserVector with the corresponding component of v.

        This vector is not altered.
        @param y The UserVector to store the result in.
        @param v The UserVector to "multiply" with.
    */
    virtual void diagmult(UserVector<Type>& y, const UserVector<Type>& v) const {
      y=diagmult(v);
    }

    /** Calculates the diagonal-multiplication.
        Multiplies each component of this UserVector with the corresponding component of v.

        This vector is not altered.
        @param v The UserVector to "multiply" with.
        @return A new DenseVector, which is the diagonal-product of this one and v.
    */
    virtual DenseVector<Type> diagmult(const UserVector<Type>& v) const {
      return DenseVector<Type>(*this).diagmult(v);
    }

    /** Square of 2-Norm.
        @return The square of the 2-norm of this UserVector.
    */
    virtual Type sq_norm2() const {
      Type ret=0;
      for (int i=0; i<dim(); i++) ret+=(*this)(i)*(*this)(i);
      return ret;
    }

    virtual Type dist(const UserVector<Type>& v) const {
			return (Type)sqrt((*this-v).sq_norm2());
    }

    /** Calculates the mean value of the elements.
        @return The arithmetic mean value.
    */
    virtual Type mean_value() const {
      Type val=0;
      for (int i=0; i<dim(); i++) val+=(*this)(i);
			val*=(Type)(1./dim());
      return val;
    }

    /** Calculates the standard deviation.
        @return The standard deviation.
        @see mean_value()
    */
    virtual Type standard_deviation() const {
      Type mval(mean_value());
      Type val=0;
      for (int i=0; i<dim(); i++) val+=(mval-(*this)(i))*(mval-(*this)(i));
			val*=(Type)(1./dim());
      return (Type)sqrt(val);
    }

    /** Sets this vector to random values.
        @param lb The lower bound of the random values.
        @param ub The upper bound of the random values.
    */
    void set_random(const Type lb, const Type ub) {
      assert(lb<=ub);
      for (int i=0; i<dim(); i++) (*this)[i]=random(lb, ub);
    }

    /** Sets this vector to random values.
        @param low The lower bound of the random vector.
        @param up The upper bound of the random vector.
    */
    void set_random(const UserVector<Type>& low, const UserVector<Type>& up) {
    	for (int i=0; i<dim(); i++) (*this)[i]=random(low(i), up(i));
    }

    /** Sets this vector to random values with a specific sparsity.
        @param lb The lower bound of the random values.
        @param ub The upper bound of the random values.
        @param sparsity The sparsity, between 0 and 1.
    */
    void set_random(const Type lb, const Type ub, const double sparsity) {
      assert(0<=sparsity && sparsity<=1);
      assert(lb<=ub);
      (*this)=0;
      int ne=(dim()*sparsity);
      for (int i=0; i<ne; i++)
        (*this)[random((int)(dim()/ne)*i,(int)(dim()/ne)*(i+1))]=random(lb, ub);
    }

    /** Prints this UserVector.
        Prints all elements of this UserVector.
        @param out The ostream to print to.
    */
    virtual void print(ostream& out) const {
      for (int i=0; i<dim(); i++) out << (*this)(i) << " ";
      out << endl;
    }

};

/** Wrapper-class for a dense TNT-Vector.
*/
template <class Type>
class DenseVector : public UserVector<Type> {
  public:
    /** Gives a random DenseVector.
        @param n The dimension.
        @param lb The lower bound.
        @param ub The upper bound.
    */
    static Pointer<DenseVector<Type> > random(const int n, const Type lb, const Type ub) {
      Pointer<DenseVector<Type> > ret=new DenseVector<Type>(n);
      ret->set_random(lb, ub);
      return ret;
    }

    /** Gives a random DenseVector with specified sparsity.
        @param n The random DenseVector.
        @param lb The lower bound.
        @param ub The upper bound.
        @param sparsity The sparsity, between 0 and 1.
    */
    static Pointer<DenseVector<Type> > random(const int n, const Type lb, const Type ub, const double sparsity) {
      Pointer<DenseVector<Type> > ret=new DenseVector<Type>(n);
      ret->set_random(lb, ub, sparsity);
      return ret;
    }

    /** Gives a random DenseVector.
        @param low The lower bound.
        @param up The upper bound.
    */
    static Pointer<DenseVector<Type> > random(const DenseVector<Type> low, const DenseVector<Type> up) {
    	Pointer<DenseVector<Type> > ret(new DenseVector<Type>(low.dim()));
    	ret->set_random(low, up);
    	return ret;
    }

  protected:
    /** The TNT-vector to wrap.
    */
    TNT::Vector<Type> x;

  public:
    /** (Standard-)Constructor for optional dimension and initial value.
        @param n The dimension, default is 0.
        @param v Initial value, default is Type().
    */
    DenseVector(int n=0, const Type& v=Type())
    : x(n, v)
    { }

    /** Copy-Constructor for a DenseVector.
        @param v The DenseVector to copy.
    */
    DenseVector(const DenseVector<Type>& v)
    : x(v.x)
    { }

    /** Copy-Constructor for a UserVector.
        @param v The UserVector to copy.
    */
    DenseVector(const UserVector<Type>& v)
    : x(v.dim())
    { for (int i=0; i<dim(); ++i) x[i]=v(i);
    }

    /** Constructor for a part of a UserVector.
        Constructs DenseVector v[low]...v[up].
        @param v The UserVector, the elements should be taken from.
        @param low The lower index of the elements.
        @param up The upper index of the elements.
    */
    DenseVector(const UserVector<Type>& v, const int low, const int up)
    : x(up-low+1)
    { for (int i=low; i<=up; i++) x[i-low]=v(i);
    }

    /** Constructor for a block of a UserVector.
        @param v The UserVector, a block should be taken from.
        @param block The block description.
    */
    DenseVector(const UserVector<Type>& v, const UserVector<int>& block)
    : x(block.dim())
    { for (int i=0; i<block.size(); i++) x[i]=v(block(i));
    }

    /** Constructor for a vector of pointer of UserVectors and a block-structure.
        @param v The vector of pointer of UserVectors.
        @param block The block-structure.
    */
    DenseVector(const vector<DenseVector<Type>* >& v, const vector<DenseVector<int> >& block)
    : x(0)
    { int s=0; for (int i=0; i<block.size(); i++) s+=block[i].size();
      resize(s);
      for (int k=0; k<block.size(); k++) if (v[k]) set_block(*v[k], block[k]);
    }

    DenseVector(const vector<Pointer<UserVector<Type> > >& v, const vector<DenseVector<int> >& block)
    : x(0)
    { int s=0; for (int i=0; i<block.size(); i++) s+=block[i].size();
      resize(s);
      for (int k=0; k<block.size(); k++) if (v[k]) set_block(*v[k], block[k]);
    }

    DenseVector(const vector<Pointer<DenseVector<Type> > >& v, const vector<DenseVector<int> >& block)
    : x(0)
    { int s=0; for (int i=0; i<block.size(); i++) s+=block[i].size();
      resize(s);
      for (int k=0; k<block.size(); k++) if (v[k]) set_block(*v[k], block[k]);
    }

    /** Constructor for a vector of UserVectors and a block-structure.
        @param v The vector of UserVectors.
        @param block The block-structure.
    */
    DenseVector(const vector<DenseVector<Type> >& v, const vector<DenseVector<int> >& block)
    : x(0)
    { int s=0; for (int i=0; i<block.size(); i++) s+=block[i].size();
      resize(s);
      for (int k=0; k<block.size(); k++) set_block(v[k], block[k]);
    }

    /** Constructor for an array of Type's.
        @param v The array of Type's.
        @param n The length of this array.
    */
    DenseVector(const Type* v, const int n)
    : x(n, v)
    { }

    /** Copy-Constructor for a TNT-Vector.
        @param x_ The TNT-Vector to copy.
    */
    DenseVector(const TNT::Vector<Type>& x_)
    : x(x_)
    { }

    Pointer<UserVector<Type> > getemptycopy(int n) const {
      return new DenseVector<Type>(n);
    }

    Pointer<UserVector<Type> > getemptycopy() const {
      return new DenseVector<Type>(dim());
    }

    Pointer<UserVector<Type> > getcopy() const {
      return new DenseVector<Type>(*this);
    }

    Pointer<UserVector<Type> > getcopy(const UserVector<int>& block) const {
      return new DenseVector<Type>(*this, block);
    }

    int dim() const {
      return x.dim();
    }

    void resize(const int n) {
      UserVector<Type>::resize(n);
      TNT::Vector<Type> tmp(x);
      x.newsize(n);
      int i=0;
      for (; i<MIN(n, tmp.size()); i++) x[i]=tmp[i];
      for (; i<n; i++) x[i]=Type();
    }

    Type& operator[](const int i) {
      return x[i];
    }

    Type operator[](const int i) const {
      return x[i];
    }

    Type operator()(const int i) const {
      return x[i];
    }

    void SetElement(const int i, const Type& v) {
      x[i]=v;
    }

    DenseVector<Type> operator()(const UserVector<int>& block) const {
      return DenseVector<Type>(*this, block);
    }

    DenseVector<Type> operator()(const int low, const int up) const {
      return DenseVector<Type>(*this, low, up);
    }

    /** Assign-Operator for a DenseVector.
        @param v The DenseVector to assign to this one.
        @return This DenseVector.
    */
    DenseVector<Type>& operator=(const DenseVector<Type>& v) {
      x=v.x;
      return *this;
    }

    DenseVector<Type>& operator=(const UserVector<Type>& v) {
      assert(dim()==v.dim());
      for (int i=0; i<dim(); i++) x[i]=v(i);
      return *this;
    }

    DenseVector<Type>& operator=(const Type& v) {
      x=v;
      return *this;
    }

    virtual UserVector<Type>& AddMult(Type a, const UserVector<Type>& v) {
      assert(dim()==v.dim());
      for (int i=0; i<dim(); i++) x[i]+=a*v(i);
      return *this;
    }

    /** Adds a DenseVector to this one.
        @param v The DenseVector to add.
        @return This DenseVector.
    */
    DenseVector<Type>& operator+=(const DenseVector<Type>& v) {
      x=x+v.x;
      return *this;
    }

    DenseVector<Type>& operator+=(const UserVector<Type>& v) {
      assert(dim()==v.dim());
      for (int i=0; i<dim(); i++) x[i]+=v(i);
      return *this;
    }

    /** Substract a DenseVector from this one.
        @param v The DenseVector to substract.
	@return This DenseVector.
    */
    DenseVector<Type>& operator-=(const DenseVector<Type>& v) {
      x=x-v.x;
      return *this;
    }

    DenseVector<Type>& operator-=(const UserVector<Type>& v) {
      assert(dim()==v.dim());
      for (int i=0; i<dim(); i++) x[i]-=v(i);
      return *this;
    }

    DenseVector<Type>& operator*=(const Type& v) {
      for (int i=0; i<dim(); i++) x[i]*=v;
      return *this;
    }

    DenseVector<Type>& operator/=(const Type& v) {
      for (int i=0; i<dim(); i++) x[i]/=v;
      return *this;
    }

    DenseVector<Type> operator*(const Type& v) const {
      DenseVector<Type> ret(*this);
      ret*=v;
      return ret;
    }

    /** Compute the product of a Type and a DenseVector.
        @param v1 The first factor as Type.
        @param v2 The secont factor as DenseVector.
        @return A DenseVector, which is the product v1 * v2.
        @see operator*(const Type&)
    */
    friend DenseVector<Type> operator*(const Type& v1, const DenseVector<Type>& v2) {
      return v2 * v1;
    }

    /** Computes the scalar product with a DenseVector.
        @param v The DenseVector to multiply with.
        @param The scalar-product x * v.x.
    */
    Type operator*(const DenseVector<Type>& v) const {
      return TNT::dot_prod(x, v.x);
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserVector<Type>::operator*;
#endif

    /** Computes the sum of this DenseVector and another one.
        @param v The DenseVector to add with this one.
        @return The sum x + v.x as DenseVector.
    */
    DenseVector<Type> operator+(const DenseVector<Type>& v) const {
      return DenseVector<Type>(x+v.x);
    }

    DenseVector<Type> operator+(const UserVector<Type>& v) const {
      DenseVector<Type> ret(*this);
      return (ret+=v);
    }

    const UserVector<Type>& operator+() const {
      return *(const DenseVector<Type>*)this;
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserVector<Type>::operator+;
#endif

    /** Computes the difference between this DenseVector and another one.
        @param v The DenseVector to subtract.
        @return The difference x - v.x as DenseVector.
    */
    DenseVector<Type> operator-(const DenseVector<Type>& v) const {
      return DenseVector<Type>(x-v.x);
    }

    DenseVector<Type> operator-(const UserVector<Type>& v) const {
      DenseVector<Type> ret(*this);
      return (ret-=v);
    }

    DenseVector<Type> operator-() const {
      DenseVector<Type> ret(*this);
      ret*=-1;
      return ret;
    }

    /** Computes the diagonal multiplication of with another DenseVector.
        @param v The DenseVector to multiply diagonaly with.
        @return The diagonal product as DenseVector.
    */
    DenseVector<Type> diagmult(const DenseVector<Type>& v) const {
      return DenseVector<Type>(x*v.x);
    }

    DenseVector<Type> diagmult(const UserVector<Type>& v) const {
      DenseVector<Type> ret(*this);
      for (int i=0; i<dim(); i++) ret[i]*=v(i);
      return ret;
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserVector<Type>::diagmult;
#endif

    /** Returns a const Pointer to the elements of this DenseVector.
    */
    operator const Pointer<Type>() const {
      return (const Pointer<Type>)Pointer<Type>(x.begin(), false);
    }

    /** Returns a Pointer to the elements of this DenseVector.
    */
/*    operator Pointer<Type>() {
      return Pointer<Type>(x.begin(), false);
    }
*/
    /** Prints this DenseVector.
        Prints each element, seperated by space and finished with a new line.
        @param out The ostream to print to.
    */
    void print(ostream& out) const {
      for (int i=0; i<dim(); i++) out << x[i] << " ";
      out << endl;
    }

};

typedef DenseVector<double> dvector;
typedef DenseVector<int> ivector;

class SparseMatrix2;
class IntervalVector;

/** Class for a sparse vector.
    The elements are stored in a simple list sorted by the indices.
*/
template <class Type>
class SparseVector : public UserVector<Type> {
  friend class SparseMatrix;
	friend class IntervalVector;
  public:
    /** Gives a random SparseVector.
        @param n The dimension.
        @param lb The lower bound.
        @param ub The upper bound.
    */
    static Pointer<SparseVector<Type> > random(const int n, const Type lb, const Type ub) {
      Pointer<SparseVector<Type> > ret=new SparseVector<Type>(n);
      ret->set_random(lb, ub);
      return ret;
    }

    /** Gives a random SparseVector with specified sparsity.
        @param n The random DenseVector.
        @param lb The lower bound.
        @param ub The upper bound.
        @param sparsity The sparsity, between 0 and 1.
    */
    static Pointer<SparseVector<Type> > random(const int n, const Type lb, const Type ub, const double sparsity) {
      Pointer<SparseVector<Type> > ret=new SparseVector<Type>(n);
      ret->set_random(lb, ub, sparsity);
      return ret;
    }

    /** The structure of one element form my vector in my list of elements.
    */
    struct VectorElement {
      /** The index of the element.
      */
      int index;
      /** The value of the element.
      */
      Type value;
      /** A pointer to the next element.
      */
      VectorElement* next;
    };

  protected:
    /** A pointer to the head of my list of elements, which is no element itselfe.
    */
    VectorElement* head;
    /** The last element from my list.
    */
    VectorElement* last;

    /** The dimenion.
    */
    int dim_;

    /** Get an element from this vector.
        Checks first, if it's behind the last element.
        If not, starts to search for the right place at the beginning.
        @param i The index of the element to set.
        @param create_if_not_exists If true (default), a new element will be added, when no element i existed.
        @return A pointer to the of the i-th element or NULL if create_if_not_exists is false and the element didn't exist.
    */
    VectorElement* GetElement(const int i, bool create_if_not_exists=true) {
      assert(i<dim());
      if (last && (last->index<i)) {
        if (! create_if_not_exists) return NULL;
        last=last->next=new VectorElement;
        last->index=i;
        last->value=Type();
        last->next=0;
        return last;
      }
      VectorElement* p=head;
      for (VectorElement* q=p->next; q; p=q, q=q->next) {
        if (q->index==i) return q; // element exists
        if (q->index>i) { // create new one, when at right position
          if (! create_if_not_exists) return NULL;
          p->next=new VectorElement;
          p->next->index=i;
          p->next->value=Type();
          p->next->next=q;
          return p->next;
        }
      } // this is only reached, when last is 0.
      if (! create_if_not_exists) return NULL;
      last=p->next=new VectorElement;
      last->index=i;
      last->value=Type();
      last->next=0;
      return last;
    }

    /** Clears all elements, starting from the next one.
        Deletes v->next and all elements, which comes later in this list.
        Sets last to v.
        @param v The VectorElement, where to start to delete the next elements.
    */
    void clear(VectorElement* v) {
			VectorElement* next=v->next;
      for (VectorElement* p=next; next; p=next) { next=next->next; delete p; }
      last=v;
      last->next=NULL;
    }

  public:
		const VectorElement* gethead() { return head; }

    /** (Standard-)Constructor for optional dimension.
        @param n The dimension, default is 0.
    */
    SparseVector(int n=0)
    : dim_(n), head(new VectorElement), last(0)
    { assert(n>=0);
      head->next=NULL;
    }

    /** Constructor for dimension and one element.
        @param n The dimension.
        @param i The index of the element.
        @param v The value of the element.
    */
    SparseVector(int n, int i, Type v)
    : dim_(n), head(new VectorElement), last(0)
    { assert(n>=0);
      head->next=NULL;
      SetElement(i, v);
    }

    /** Constructor for dimension and two elements.
        @param n The dimension.
        @param i First index.
        @param j Second index.
        @param vi First value.
        @param vj Second value.
    */
    SparseVector(int n, int i, int j, Type vi, Type vj)
    : dim_(n), head(new VectorElement), last(0)
    { assert(n>=0);
      head->next=NULL;
      SetElement(i, vi);
      SetElement(j, vj);
    }

    /** Copy-Constructor for a SparseVector.
        @param v The SparseVector to copy.
    */
    SparseVector(const SparseVector<Type>& v)
    : dim_(v.dim()), head(new VectorElement), last(0)
    { head->next=NULL;
      for (VectorElement* p=v.head->next; p; p=p->next) {
        if (! p->value) continue;
        if (last) last=last->next=new VectorElement;
        else { last=new VectorElement; last->next=NULL; }
        if (!head->next) head->next=last;  // only in first iteration
        last->index=p->index;
        last->value=p->value;
      }
      if (last) last->next=0;
    }

    /** Copy-Constructor for a UserVector.
        @param v The UserVector to copy.
    */
    SparseVector(const UserVector<Type>& v)
    : dim_(v.dim()), head(new VectorElement), last(0)
    { head->next=NULL;
      double value;
      for (int i=0; i<v.dim(); i++) {
        if (fabs(value=v(i))<rtol) continue;
        if (last) last=last->next=new VectorElement;
        else { last=new VectorElement; last->next=NULL; }
        if (!head->next) head->next=last;
        last->index=i;
        last->value=value;
      }
      if (last) last->next=0;
    }

    /** Copy-Constructor for a part of a UserVector.
        @param v The UserVector to copy the part [low,up] from.
        @param low The lower index of the range of elements to copy.
        @param up The upper index.
    */
    SparseVector(const UserVector<Type>& v, const int low, const int up)
    : dim_(up-low+1), head(new VectorElement), last(NULL)
    { head->next=NULL;
      double value;
      for (int i=low; i<=up; i++) {
        if (fabs(value=v(i))<rtol) continue;
        if (last) last=last->next=new VectorElement;
        else { last=new VectorElement; last->next=NULL; }
        if (!head->next) head->next=last;
        last->index=i-low;
        last->value=value;
      }
      if (last) last->next=0;
    }

    /** Copy-Constructor for a UserVector and a block-description.
        @param v The UserVector to take one the block from.
        @param block The block-describtion.
    */
    SparseVector(const UserVector<Type>& v, const UserVector<int>& block)
    : dim_(block.size()), head(new VectorElement), last(NULL)
    { head->next=NULL;
    	for (int i=0; i<block.size(); i++) SetElement(i, v(block(i)));
    }

    /** Copy-Constructor for a vector of UserVectors and a vector of block-describtions.
        @param v The vector of UserVectors.
        @param block The block-describtion.
    */
    SparseVector(const vector<Pointer<UserVector<Type> > >& v, const vector<DenseVector<int> >& block)
    : dim_(0), head(new VectorElement), last(NULL)
    { head->next=NULL;
    	for (int k=0; k<block.size(); k++) dim_+=block[k].size();
			for (int k=0; k<block.size(); k++)
      	if (v[k])
					for (int i=0; i<block[k].size(); i++) SetElement(block[k](i), (*v[k])(i));
    }

    /** Copy-Constructor for an array of Types.
        @param v The array of Types.
        @param n The length of v.
    */
    SparseVector(const Type* v, const int& n)
    : dim_(n), head(new VectorElement), last(head)
    { head->next=NULL;
      for (int i=0; i<n; i++) {
        if (fabs(v[i])<rtol) continue;
        last=last->next=new VectorElement;
        if (!head->next) head->next=last;
        last->index=i;
        last->value=v[i];
      }
      if (last) last->next=NULL;
    }
    
    SparseVector(const int& n, const int& nz, const int* indices, const double* elements)
    : dim_(n), head(new VectorElement), last(head)
    { head->next=NULL;
      for (int i=0; i<nz; ++i) {
      	if (i) assert(indices[i-1]<indices[i]);
        last=last->next=new VectorElement;
        if (!head->next) head->next=last;
        last->index=indices[i];
        last->value=elements[i];
      }
      if (last) last->next=NULL;
    }

    Pointer<UserVector<Type> > getemptycopy(int n) const {
      return new SparseVector<Type>(n);
    }

	Pointer<UserVector<Type> > getemptycopy() const {
		return new SparseVector<Type>(dim());
	}

    Pointer<UserVector<Type> > getcopy() const {
      return new SparseVector<Type>(*this);
    }

    Pointer<UserVector<Type> > getcopy(const UserVector<int>& block) const {
      return new SparseVector<Type>(*this, block);
    }

    /** Destructor.
        Deletes the list by calling clear(head).
        @see clear(VectorElement*)
    */
    ~SparseVector() {
       clear(head);
       delete head;
    }

    /** Gives the dimension.
        @return The dimension dim_.
        @see resize(int)
    */
    int dim() const { return dim_; }

    /** Resize this vector.
        Deletes all elements with index greater equal n.
        @param n The new dimension.
        @see dim()
    */
    void resize(const int n) {
      UserVector<Type>::resize(n);
      dim_=n;
      if ((last && last->index<n) || (!head->next)) return;
      for (VectorElement* p=head; p; p=p->next)
        if (p->next && p->next->index>=n) {
          clear(p);
          last=p;
          return;
        }
    }

    /** Sets one element of this vector.
        If the element isn't existing, it's created.
        If the element exists and set to zero and test_for_zero is true, it's deleted.
        @param i The index of the element to set.
        @param v The value to set.
        @param test_for_zero If true, a new element is only created, if it's not zero.
    */
    void SetElement(int i, const Type& v, bool test_for_zero) {
      if (test_for_zero && fabs(v)<rtol) {
        DelElement(i);
        return;
      }
      GetElement(i)->value=v;
    }

    void SetElement(const int i, const Type& v) { SetElement(i, v, true); }

    /** Deletes one element of this vector.
        @param i The index of the element to delete.
    */
    void DelElement(int i) {
      for (VectorElement* p=head->next, *q=head; p && p->index<=i; q=p, p=p->next)
        if (p->index==i) {
          q->next=p->next;
          delete p;
          if (last==p) last=q->next;
          return;
        }
    }

    Type& operator[](const int i) {
      return GetElement(i)->value;
    }

    Type operator[](const int i) const {
      return (*this)(i);
    }

    Type operator()(const int i) const {
      assert(i<dim());
      if (last && (last->index<i)) return Type();
      for (VectorElement* p=head->next; p && p->index<=i; p=p->next)
        if (p->index==i) return p->value;
      return Type();
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserVector<Type>::operator();
#endif

    // could be improved
    UserVector<Type>& operator=(const UserVector<Type>& v) {
      assert(v.dim()==dim());
      for (int i=0; i<v.dim(); i++) SetElement(i, v(i));
      return *this;
    }

    /** Assign-Operator for a SparseVector.
        @param v The SparseVector to assign.
        @return This SparseVector.
    */
    SparseVector<Type>& operator=(const SparseVector<Type>& v) {
      assert(v.dim()==dim());
			clear(head);
      for (VectorElement* p=v.head->next; p; p=p->next) {
        if (! p->value) continue;
        last=last->next=new VectorElement; // last = head after clear
        last->index=p->index;
        last->value=p->value;
      }
			if (last==head) last=NULL;
			else last->next=NULL;
/*      VectorElement* p=head;
      VectorElement* q=v.head->next;
      for (; p->next && q; )
        if (p->next->index == q->index) {
          p->next->value=q->value;
          p=p->next; q=q->next;
        }
        else if (p->next->index < q->index) {
          VectorElement* r=p->next;
          p=p->next=r->next;
          delete r;
        }
        else {
          VectorElement* r=p->next;
          p=p->next=new VectorElement;
          p->index=q->index;
          p->value=q->value;
          p->next=r;
          q=q->next;
        }
      for (; q; q=q->next) {
        p=p->next=new VectorElement;
        p->index=q->index;
        p->value=q->value;
      }
      clear(p);
*/      return *this;
    }

    UserVector<Type>& operator=(const Type& v) {
      if (fabs(v)<rtol) {
        clear(head); last=NULL;
        return *this;
      }
      for (int i=0; i<dim(); i++) SetElement(i, v, false);
      return *this;
    }

    void set(const Type* v, int n) {
      clear(head);
      VectorElement* p=head;
      for (int i=0; i<n; i++, v++)
        if (fabs(*v)>rtol) {
          p=p->next=new VectorElement;
          p->index=i;
          p->value=*v;
        }
      p->next=NULL;
      last=p;
    }

    void set(const Type* v) { set(v, dim()); }

    UserVector<Type>& AddMult(Type a, const UserVector<Type>& v) {
      assert(dim()==v.dim());
      int i=0;
      VectorElement* p=head;
      if (p->next)
        for (; p->next->index < i && i < dim(); i++)
          if (fabs(v(i))>0) {
            VectorElement* r=p->next;
            p=p->next=new VectorElement;
            p->index=i;
            p->value=a*v(i);
            p->next=r;
          }
      for (; p->next && i<v.dim(); i++) {
        if (p->next->index == i) {
          p->next->value+=a*v(i);
          p=p->next;
        }
        else if (fabs(v(i))>0) {  // here is p->next->index > i
          VectorElement* r=p->next;
          p=p->next=new VectorElement;
          p->index=i;
          p->value=a*v(i);
          p->next=r;
        }
      }
      for (; i<v.dim(); i++)
        if (fabs(v(i))>0) {
          p=p->next=new VectorElement;
          p->index=i;
          p->value=a*v(i);
        }
      p->next=NULL;
			if (p!=head) last=p;
      return *this;
    }

    UserVector<Type>& operator+=(const UserVector<Type>& v) {
      assert(dim()==v.dim());
      int i=0;
      VectorElement* p=head;
      if (p->next)
        for (; p->next->index < i && i < dim(); i++)
          if (fabs(v(i))>0) {
            VectorElement* r=p->next;
            p=p->next=new VectorElement;
            p->index=i;
            p->value=v(i);
            p->next=r;
          }
      for (; p->next && i<v.dim(); i++) {
        if (p->next->index == i) {
          p->next->value+=v(i);
          p=p->next;
        }
        else if (fabs(v(i))>0) {  // here is p->next->index > i
          VectorElement* r=p->next;
          p=p->next=new VectorElement;
          p->index=i;
          p->value=v(i);
          p->next=r;
        }
      }
      for (; i<v.dim(); i++)
        if (fabs(v(i))>0) {
          p=p->next=new VectorElement;
          p->index=i;
          p->value=v(i);
        }
      p->next=NULL;
      last=p;
      return *this;
    }

    SparseVector<Type>& AddMult(Type a, const SparseVector<Type>& v) {
      assert(dim()==v.dim());
      VectorElement* p=head;
      VectorElement* q=v.head->next;
      while(p->next && q)
        if (p->next->index == q->index) {
          p->next->value+=a*q->value;
          p=p->next; q=q->next;
        }
        else
          if (p->next->index < q->index) p=p->next;
        else {
          VectorElement* r=p->next;
          p=p->next=new VectorElement;
          p->index=q->index;
          p->value=a*q->value;
          p->next=r;
          q=q->next;
        }
      for(; q; q=q->next) {
        p=p->next=new VectorElement;
        p->index=q->index;
        p->value=a*q->value;
        p->next=NULL;
      }
      if (!p->next) last=p;
      return *this;
    }

    /** Adds a SparseVector to this one.
        @param v The SparseVector to add.
        @return This SparseVector.
    */
    SparseVector<Type>& operator+=(const SparseVector<Type>& v) {
      assert(dim()==v.dim());
      VectorElement* p=head;
      VectorElement* q=v.head->next;
      while(p->next && q)
        if (p->next->index == q->index) {
          p->next->value+=q->value;
          p=p->next; q=q->next;
        }
        else
          if (p->next->index < q->index) p=p->next;
        else {
          VectorElement* r=p->next;
          p=p->next=new VectorElement;
          p->index=q->index;
          p->value=q->value;
          p->next=r;
          q=q->next;
        }
      for(; q; q=q->next) {
        p=p->next=new VectorElement;
        p->index=q->index;
        p->value=q->value;
        p->next=NULL;
      }
      if (!p->next) last=p;
      return *this;
    }

    UserVector<Type>& operator-=(const UserVector<Type>& v) {
      assert(dim()==v.dim());
      int i=0;
      VectorElement* p=head;
      if (p->next)
        for (; p->next->index < i && i < dim(); i++)
          if (fabs(v(i))>0) {
            VectorElement* r=p->next;
            p=p->next=new VectorElement;
            p->index=i;
            p->value=-v(i);
            p->next=r;
          }
      for (; p->next && i<v.dim(); i++) {
        if (p->next->index == i) {
          p->next->value-=v(i);
          p=p->next;
        }
        else if (fabs(v(i))>0) {  // here is p->next->index > i
          VectorElement* r=p->next;
          p=p->next=new VectorElement;
          p->index=i;
          p->value=-v(i);
          p->next=r;
        }
      }
      for (; i<v.dim(); i++)
        if (fabs(v(i))>0) {
          p=p->next=new VectorElement;
          p->index=i;
          p->value=-v(i);
        }
      p->next=NULL;
      last=p;
      return *this;
    }

    /** Substracts a SparseVector to this one.
        @param v The SparseVector to substract.
    */
    SparseVector<Type>& operator-=(const SparseVector<Type>& v) {
      assert(dim()==v.dim());
      VectorElement* p=head;
      VectorElement* q=v.head->next;
      while(p->next && q)
        if (p->next->index == q->index) {
          p->next->value-=q->value;
          p=p->next; q=q->next;
        }
        else
          if (p->next->index < q->index) p=p->next;
        else {
          VectorElement* r=p->next;
          p=p->next=new VectorElement;
          p->index=q->index;
          p->value=-q->value;
          p->next=r;
          q=q->next;
        }
      for(; q; q=q->next) {
        p=p->next=new VectorElement;
        p->index=q->index;
        p->value=-q->value;
        p->next=NULL;
      }
      if (!p->next) last=p;
      return *this;
    }

    UserVector<Type>& operator*=(const Type& v) {
      if (v==0.) {
        clear(head);
        last=NULL;
        return *this;
      }
      for (VectorElement* p=head->next; p; p=p->next) p->value*=v;
      return *this;
    }

    UserVector<Type>& operator/=(const Type& v) {
      for (VectorElement* p=head->next; p; p=p->next) p->value/=v;
      return *this;
    }

    DenseVector<Type> operator*(const Type& v) const {
      DenseVector<Type> ret(dim());
      for (VectorElement* p=head->next; p; p=p->next) ret[p->index]=v*p->value;
      return ret;
    }

    bool operator>=(const Type& v) const {
      if (v<-rtol && last && last->index<dim()-1) return false;
      for (VectorElement* p=head->next; p; p=p->next) if (p->value+rtol<v) return false;
      return true;
    }

    bool operator==(const Type& v) const {
      if (fabs(v)>rtol && last && last->index<dim()-1) return false;
      for (VectorElement* p=head->next; p; p=p->next)
        if (fabs(p->value-v)>rtol) return false;
      return true;
    }

    bool operator==(const UserVector<Type>& v) const {
      assert(v.dim()==dim());
      VectorElement* p=head->next;
      int i=0;
      if (p)
        while (p->index>i && i<dim())
          if (fabs(v(i++))>rtol) return false;
      for (; p && i<dim(); i++)
        if (p->index==i)
          if (fabs(p->value-v(i))>=rtol) return false;
          else p=p->next;
        else
          if (fabs(v(i))>=rtol) return false;
      while (i<dim()) if (fabs(v(i++))>=rtol) return false;
      return true;
    }

    Type operator*(const UserVector<Type>& v) const {
      assert(v.dim()==dim());
      Type val=0;
      for (VectorElement* p=head->next; p; p=p->next)
        val+=p->value*v(p->index);
      return val;
    }

    /** Computes the scalar product for this and another SparseVector.
        @param v The sparse vector to multiply with.
        @return The scalar product of this vector and the other one.
    */
    Type operator*(const SparseVector<Type>& v) const {
      assert(v.dim()==dim());
      Type val=0;
      for (VectorElement* p=head->next, *q=v.head; p && q; )
        if (p->index==q->index) {
          val+=p->value * q->value;
          p=p->next; q=q->next;
        }
        else if (p->index < q->index) p=p->next;
        else q=q->next;
      return val;
    }

    /** Computes the product of a Type and a SparseVector.
        @param v1 The first factor as Type.
        @param v2 The second factor as SparseVector.
        @return The product v1*v2.
    */
    friend SparseVector<Type> operator*(const Type& v1, const SparseVector<Type>& v2) {
      SparseVector<Type> ret(v2);
      for (VectorElement* p=ret.head->next; p; p=p->next)
				p->value*=v1;
			return ret;
    }

    /** Computes the sum of this SparseVector and another one.
        @param v The SparseVector to add.
        @return A new SparseVector, which is this one plus the given one.
    */
    SparseVector<Type> operator+(const SparseVector<Type>& v) const {
      assign(v.dim()==dim());
      VectorElement* p1=head->next;
      VectorElement* p2=v.head->next;
      SparseVector<Type> ret(dim());
      ret.last=ret.head;
      while(p1 && p2) {
        ret.last=ret.last->next=new VectorElement;
        if (p1->index==p2->index) {
          ret.last->index=p1->index;
          ret.last->value=p1->value+p2->value;
          p1=p1->next; p2=p2->next;
        }
        else if (p1->index<p2->index) {
          ret.last->index=p1->index;
          ret.last->value=p1->value;
          p1=p1->next;
        }
        else {
          ret.last->index=p2->index;
          ret.last->value=p2->value;
          p2=p2->value;
        }
      }
      for (p1+=p2; p1; p1=p1->next) {  // at least one of them is NULL, so continue with the one, which is not NULL
        ret.last=ret.last->next=new VectorElement;
        ret.last->index=p1->index;
        ret.last->value=p1->value;
      }
      return ret;
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserVector<Type>::operator+;
#endif

    /** Computes the difference of this SparseVector and another one.
        @param v The SparseVector to substract.
        @return The difference as SparseVector.
    */
    SparseVector<Type> operator-(const SparseVector<Type>& v) const {
      assign(v.dim()==dim());
      VectorElement* p1=head->next;
      VectorElement* p2=v.head->next;
      SparseVector<Type> ret(dim());
      if (p1 || p2) ret.last=ret.head;
      else return ret;
      while(p1 && p2) {
        ret.last=ret.last->next=new VectorElement;
        if (p1->index==p2->index) {
          ret.last->index=p1->index;
          ret.last->value=p1->value-p2->value;
          p1=p1->next; p2=p2->next;
        }
        else if (p1->index<p2->index) {
          ret.last->index=p1->index;
          ret.last->value=p1->value;
          p1=p1->next;
        }
        else {
          ret.last->index=p2->index;
          ret.last->value=-p2->value;
          p2=p2->value;
        }
      }
      for (; p1; p1=p1->next) {
        ret.last=ret.last->next=new VectorElement;
        ret.last->index=p1->index;
        ret.last->value=p1->value;
      }
      for (; p2; p2=p2->next) {
        ret.last=ret.last->next=new VectorElement;
        ret.last->index=p2->index;
        ret.last->value=-p2->value;
      }
      return ret;
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserVector<Type>::operator-;
#endif

    /** Computes the diagonal multiplication of this SparseVector and another one.
        @param v The SparseVector to multiply with.
        @return The diagonal product.
    */
    void diagmult(SparseVector<Type> y, const SparseVector<Type>& v) const {
      assert(dim()==v.dim());
      assert(dim()==y.dim());
      VectorElement* p1=head->next;
      VectorElement* p2=v.head->next;
      y=0;
      if (p1 && p2) y.last=y.head;
      else return;
      while (p1 && p2)
        if (p1->index==p2->index) {
          y.last=y.last->next=new VectorElement;
          y.last->index=p1->index;
          y.last->value=p1->value * p2->value;
          p1=p1->next; p2=p2->next;
        }
        else if (p1->index < p2->index)
          p1=p1->next;
        else p2=p2->next;
    }

    void diagmult(UserVector<Type>& y, const UserVector<Type>& v) const {
      assert(dim()==v.dim());
      assert(dim()==y.dim());
      y=0;
      for (VectorElement* p=head->next; p; p=p->next)
        y[p->index]=p->value * v(p->index);
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserVector<Type>::diagmult;
#endif

    void set_block(const UserVector<Type>& v, const UserVector<int>& block) {
      for (int i=0; i<block.size(); i++) SetElement(block(i), v(i));
    }

    operator const Pointer<Type>() const {
      Pointer<Type> ret(new Type[dim()]);
      VectorElement* p=head->next;
      for (int i=0; i<dim(); i++)
        if (p && p->index==i) {
          ret[i]=p->value;
          p=p->next;
        } else ret[i]=Type();
      return ret;
    }

    Type sq_norm2() const {
      Type val=0;
      for (VectorElement* p=head->next; p; p=p->next) val+=p->value*p->value;
      return val;
    }

    void print(ostream& out) const {
      for (VectorElement* p=head->next; p; p=p->next) out << p->index << ": " << p->value << " ";
      out << endl;
    }

};


#ifdef FILIB_AVAILABLE
class IntervalVector : public DenseVector<interval<double> > {
	public:
		IntervalVector(int n=0, const double& v=0.)
		: DenseVector<interval<double> >(n, v)
		{ }

		IntervalVector(const IntervalVector& v)
		: DenseVector<interval<double> >(v)
		{ }

    IntervalVector(const UserVector<double>& v)
    : DenseVector<interval<double> >(v.dim())
    { for (int i=0; i<dim(); i++) x[i]=v(i);
    }

		IntervalVector(const UserVector<double>& v1, const UserVector<double>& v2)
    : DenseVector<interval<double> >(v1.dim())
    { for (int i=0; i<dim(); i++) x[i]=interval<double>(v1(i), v2(i));
    }

    IntervalVector(const IntervalVector& v, const int low, const int up)
    : DenseVector<interval<double> >(v, low, up)
    { }

    IntervalVector(const IntervalVector& v, const UserVector<int>& block)
    : DenseVector<interval<double> >(v, block)
    { }

    IntervalVector(const vector<IntervalVector>& v, const vector<DenseVector<int> >& block)
    {	int s=0; for (int i=0; i<block.size(); i++) s+=block[i].size();
      resize(s);
      for (int k=0; k<block.size(); k++) set_block(v[k], block[k]);
		}

		IntervalVector& operator=(const IntervalVector& v) {
			DenseVector<interval<double> >::operator=(v);
			return *this;
		}

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using DenseVector<interval<double> >::operator=;
#endif

		interval<double> operator*(const UserVector<double>& v) const {
      assert(v.dim()==dim());
      interval<double> val;
      for (int i=0; i<dim(); i++) val+=(*this)(i)*v(i);
      return val;
    }

		interval<double> operator*(const SparseVector<double>& v) const {
      assert(v.dim()==dim());
      interval<double> val;
      for (SparseVector<double>::VectorElement* p=v.head->next; p; p=p->next)
        val+=(*this)(p->index)*p->value;
      return val;
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using DenseVector<interval<double> >::operator*;
#endif

		void AddMult(const double& a, const UserVector<double>& v) {
			for (int i=0; i<v.dim(); i++) (*this)[i]+=a*v(i);
		}

		void AddMult(const double& a, const SparseVector<double>& v) {
			for (SparseVector<double>::VectorElement* p=v.head->next; p; p=p->next)
    		(*this)[p->index]+=a*p->value;
		}

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using DenseVector<interval<double> >::AddMult;
#endif

    IntervalVector diagmult(const UserVector<double>& v) const {
      IntervalVector ret(*this);
      for (int i=0; i<dim(); i++) ret[i]*=v(i);
      return ret;
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using DenseVector<interval<double> >::diagmult;
#endif
};
#endif

/** Class to project a point onto a box.
*/
class Project {
	public:
		/** Projects a point to a box.
		    @param x The point to project.
		    @param lower The lower bound of the box.
		    @param upper The upper bound of the box.
		    @return The projected point.
		    @see project(dvector&, const dvector&, const dvector&, const dvector&)
		*/
		static dvector project(const dvector& x, const dvector& lower, const dvector& upper) {
			dvector xp(x.dim());
			project(xp, x, lower, upper);
			return xp;
		}
			
		/** Projects a point to a box.
		    @param xp The dvector to store the projected point in.
		    @param x The point to project.
		    @param lower The lower bound of the box.
		    @param upper The upper bound of the box.
		    @see project(const dvector&, const dvector&, const dvector&)
		*/
		static void project(dvector& xp, const dvector& x, const dvector& lower, const dvector& upper) {
	  	for (int j=0; j<x.dim(); j++)
			  if (x(j) < lower(j)) xp[j]=lower(j);
			  else if (x(j) > upper(j)) xp[j]=upper(j);
			  else xp[j]=x(j);
		}
		
};

/** Class to round the discrete variables of a point to the bounds.
*/
class Round {
	public:
		/** Rounds the values of the SamplePoints, which are determined by indices, to the bounds.
		    @param x The point to round.
		    @param indices The indices of the elements, which should be rounded.
		    @param lower The lower bound of the box.
		    @param upper The upper bound of the box.
		    @return The rounded point.
		    @see round(dvector&, const dvector&, vector<int>&, const dvector&, const dvector&)
		*/
		static dvector round(const dvector& x, vector<int>& indices, const dvector& lower, const dvector& upper) {
			dvector xr(x.dim());
			round(xr, x, indices, lower, upper);
			return xr;
		}
		
		/** Rounds the values of the SamplePoints, which are determined by indices, to the bounds.
		    @param xr The dvector to store the rounded point in.
		    @param x The point to round.
		    @param indices The indices of the elements, which should be rounded.
		    @param lower The lower bound of the box.
		    @param upper The upper bound of the box.
		    @see round(const dvector&, vector<int>&, const dvector&, const dvector&)
		*/
		static void round(dvector& xr, const dvector& x, vector<int>& indices, const dvector& lower, const dvector& upper) {
			xr=x;
	  	for (vector<int>::iterator i(indices.begin()); i!=indices.end(); i++)
			  xr[*i] = x(*i)<=0.5*(lower(*i)+upper(*i)) ? lower(*i) : upper(*i);
		}
};


#endif // USERVECTOR_H
