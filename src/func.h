// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef FUNC_H
#define FUNC_H

#include "standard.h"
#include "usermatrix.h"

class Func;
class Decomposition;

/** A class to represent the sparsity structure of a function.
    Depending on the preprocessing progress, not all members have to been set.
		The coefficient information is not trustable yet. Due to some design-drawbacks, it will become outdated quite fast.
*/
class SparsityInfo {
	friend class Decomposition;
	private:
		void add_nonlinear(int index) {
			nonlinear->insert(pair<int, NonlinearVariable>(index, NonlinearVariable()));
		}

		void add_quadratic(int index, double coeff_lin, double coeff_quad);

		void add_linear(int index, double coeff) {
			linear->insert(pair<int, LinearVariable>(index, LinearVariable(coeff)));
		}

		void add_sparsity_pattern(int index1, int index2, double coeff);

	public:
		/** To represent a linear variable.
		*/
		class LinearVariable {
			public:
				double coeff;
				LinearVariable(double coeff_)
				: coeff(coeff_)
				{ }

				friend ostream& operator<<(ostream& out, const LinearVariable& var) {
					out << var.coeff;
					return out;
				}
		};

		/** To represent a nonlinear variable.
		*/
		class NonlinearVariable {
			public:
		};

		/** To represent a quadratic variable.
		*/
		class QuadraticVariable {
			public:
				/** The coefficient before this variable in the linear and quadratic term: coeff_lin * x + coeff_quad * x^2
				*/
				double coeff_lin, coeff_quad;

				QuadraticVariable(double coeff_lin_, double coeff_quad_)
				: coeff_lin(coeff_lin_), coeff_quad(coeff_quad_)
				{ }

				friend ostream& operator<<(ostream& out, const QuadraticVariable& var) {
					out << var.coeff_lin << ',' << var.coeff_quad;
					return out;
				}
		};

		/** To represent a connection (entry in hessian) between two nonlinear variables.
		*/
		class NonlinearConnection {
			public:
				/** If this relates to a connection between quadratic variable, this is the coefficient for the corresponding quadratic term.
				*/
				double coeff;

				NonlinearConnection(double coeff_)
				: coeff(coeff_)
				{ }

				friend ostream& operator<<(ostream& out, const NonlinearConnection& nc) {
					out << nc.coeff;
					return out;
				}
		};

		/** The indices of the linear variables.
		*/
		Pointer<map<int, SparsityInfo::LinearVariable> > linear;
		/** The indices of the nonlinear variables.
		*/
		Pointer<map<int, SparsityInfo::NonlinearVariable> > nonlinear;

		/** The indices of the quadratic variables.
		    If set, then this set is a subset of the set of nonlinear variables.
		*/
		Pointer<map<int, SparsityInfo::QuadraticVariable> > quadratic;

		Pointer<map<pair<int, int>, SparsityInfo::NonlinearConnection> > sparsity_pattern;

		SparsityInfo(int init_level=0);

		SparsityInfo(const SparsityInfo& si);

		/** Joins two sparsity infos into one new.
		    Used when adding two functions, for example.
		*/
		SparsityInfo(const SparsityInfo& si1, const SparsityInfo& si2);

		friend ostream& operator<<(ostream& out, const SparsityInfo& si);

		bool compute_sparsity_pattern(const Func& f, const vector<dvector>& sample_set);

//		void compute_sparsity_pattern(const UserMatrix& A);

		void add(const UserVector<double>& b);
		void add(const UserVector<double>& b, const ivector& block);
		void add(const UserMatrix& A);
		void add(const UserMatrix& A, const ivector& block);
		void add(const SparsityInfo& si);
		void add(const SparsityInfo& si, const ivector& block);

		int nr_var(bool linear_=true, bool nonlinear_=true) {
			return (linear_ ? linear->size() : 0)+(nonlinear_ ? nonlinear->size() : 0);
		}
};

/** Abstract class for an iterator, which iterates over the variables of a function.
    Might iterate only over linear or only quadratic or only nonlinear/nonquadratic... variables.
*/
class VariableIterator_Type {
	protected:
		/** Indicates, whether we want to iterate over linear variable.
		*/
		bool linear;
		/** Indicates, whether we want to iterate over nonlinear variables.
		*/
		bool nonlinear;
		/** Indicates, whether we want to iterate over quadratic variables.
		    If nonlinear is set to true, the iterator will also iterate over quadratic variables.
		*/
		bool quadratic;

	public:
		typedef enum { LINEAR=0, NONLINEAR=1, QUADRATIC=3, END=4 } VarType;

		VariableIterator_Type(bool linear_=true, bool nonlinear_=true, bool quadratic_=false)
		: linear(linear_), nonlinear(nonlinear_), quadratic(quadratic_)
		{ }

		virtual ~VariableIterator_Type() { }

		/** Returns the index of the current variable.
		*/
		virtual int operator()() const=0;

		/** Returns the coefficient in the linear part of this variable, if available.
		    @return If the variable is linear, it returns it's coefficient. If the variable x is quadratic, it returns the coefficient before x.
		*/
		virtual double coeff_lin() const=0;

		/** Returns the coefficient in the quadratic part of this variable, if available.
		    @return If the variable x is quadratic, it returns the coefficient before x^2.
		*/
		virtual double coeff_quad() const=0;

		/** The type of the current variable.
		*/
		virtual VarType type() const=0;

		/** Moves to next variable.
		*/
		virtual void operator++()=0;

		/** Returns true, if we are not at the end of the variables.
		*/
		virtual operator bool() const=0;
};

class VariableIterator : public VariableIterator_Type {
	private:
		const SparsityInfo& sparsity;
		map<int, SparsityInfo::LinearVariable>::iterator it_linear;
		map<int, SparsityInfo::NonlinearVariable>::iterator it_nonlinear;
		map<int, SparsityInfo::QuadraticVariable>::iterator it_quadratic;
		/** Over which variable type am I just iterating?
		*/
		VarType whereami;

		void init();
	public:
		VariableIterator(const Func& f_, bool linear_=true, bool nonlinear_=true, bool quadratic_=false);
		VariableIterator(const SparsityInfo& sparsity_, bool linear_=true, bool nonlinear_=true, bool quadratic_=false);

		virtual int operator()() const;

		virtual double coeff_lin() const;
		virtual double coeff_quad() const;

		virtual VarType type() const { return whereami; }

		virtual void operator++();

		virtual operator bool() const { return (whereami!=END); }
};

/** The abstract base class for a function.
    You need to implement the following methods:
    - double eval(const dvector&)
    - void grad(dvector&, const dvector&)
    - void HessMult(dvector&, const dvector&, const dvector&)
		- void set_curvature(CurvatureType)
		- CurvatureType get_curvature
*/
class Func {
		/** Prints some status information.
		    Calls a.print(out).
		    @param out The ostream to print to.
		    @param a The Func to print.
		    @return The ostream.
		    @see print(ostream &)
		*/
		friend ostream& operator << (ostream& out, const Func& a)
		{  a.print(out); return out; };

		friend class VariableIterator;
		friend class MinusFunc;
		friend class SumFunc;

	protected:
		/** The dimension.
		*/
		int dim_;

		Pointer<SparsityInfo> sparsity;

		virtual SparsityInfo& get_sparsity() { assert(sparsity); return *sparsity; }

	public:
		/** A Pointer to an output-stream for function-relevant output.
		    Set in the constructor, default is out_out_p.
		    @see out_func
		*/
		Pointer<ostream> out_func_p;

      /** The dereferenced Pointer<ostream> out_func_p.
          If out_func_p==NULL, no output appears.
          Only use this in the form "out_func << ...", not as argument !
          @see Func::out_func_p
      */
#define out_func if (out_func_p) (*out_func_p)
      /** A Pointer to an ostream to print function-relevant logging-information.
          Set in the constructor, default is out_log_p.
          @see out_func_log
      */
      Pointer<ostream> out_func_log_p;
      /** The dereferenced Pointer<ostream> out_func_log_p.
          If out_func_log_p==NULL, no output appears.
          Only use this in the form "out_func_log << ...", not as argument !
          @see Func::out_func_log_p
      */
#define out_func_log if (out_func_log_p) (*out_func_log_p)

      // methods:
      /** (Standard-)Constructor for optional dimension.
          @param n The dimension, default is 0.
          @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
          @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
      */
      Func(int n=0, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : dim_(n), out_func_p(out_func_p_), out_func_log_p(out_func_log_p_)
      { };

      /** Virtual Destructor.
      */
      virtual ~Func() { }

      /** Gives the dimension.
          @return The dimension.
          @see dim_
      */
      int dim() const { return dim_; };

			/** To check, whether the sparsity information is set.
					@return True, if sparsity information is available. False, else.
			*/
			virtual bool sparsity_available() const { return sparsity; }

			/** Gives the sparsity information with read-permissions.
			*/
			virtual const SparsityInfo& get_sparsity() const { assert(sparsity); return *sparsity; }

			virtual bool compute_sparsity_pattern(const vector<dvector>& sample_set) { return get_sparsity().compute_sparsity_pattern(*this, sample_set); }

      /** Computes the value of this function for a UserVector<double>.
          Abstract.
          @param x The UserVector<double> to compute the value for.
          @return The value this(x) as double.
      */
      virtual double eval(const UserVector<double>& x) const=0;

      /** Computes the (sub)gradient for a UserVector<double>.
          Abstract.
          @param g The UserVector<double> to store the result in.
          @param x The UserVector<double> to compute the gradient for.
          @see grad(const dvector&)
          @see valgrad(double&, UserVector<double>&, const UserVector<double>&)
          @see add_grad(UserVector<double>&, UserVector<double>&)
          @see add_valgrad(double&, UserVector<double>&, UserVector<double>&)
      */
      virtual void grad(UserVector<double>& g, const UserVector<double>& x) const=0;

      /** Computes the (sub)gradient for a UserVector<double>.
          @param x The UserVector<double> to compute.
          @return The gradient as dvector.
          @see grad(UserVector<double>&, UserVector<double>&)
          @see add_grad(UserVector<double>&, const UserVector<double>&)
          @see add_valgrad(UserVector<double>&, UserVector<double>&, const UserVector<double>&)
      */
      virtual dvector grad(const dvector& x) const {
        dvector g(x.dim());
        grad(g, x);
        return g;
      }

      /** Evaluates the product of the Hessian and a UserVector<double>.
          Abstract.
          @param y The UserVector<double> to store the result in: @f$ y = \nabla^2 f(x) \cdot z @f$.
          @param x The UserVector<double> to evaluate the Hessian for.
          @param z The UserVector<double> to multiply with the Hessian.
          @see HessMult(const UserVector<double>&, const UserVector<double>&)
      */
      virtual void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const=0;

      /** Evaluates the product of the Hessian and a dvector.
          @param x The dvector to evaluate the Hessian for.
          @param z The dvector to multiply with the Hessian.
          @return The result as a dvector: @f$ \nabla^2 f(x) \cdot z @f$.
          @see HessMult(UserVector<double>&, const UserVector<double>&, const UserVector<double>&)
      */
      virtual dvector HessMult(const dvector& x, const dvector& z) const {
        dvector y(x.size());
        HessMult(y, x, z);
        return y;
      }

      /** Computes the (sub)gradient and the value of this function for a UserVector<double>.
          @param val The double to store the value this(x).
          @param g The UserVector<double> to store the gradient.
          @param x The UserVector<double> to compute the gradient and value for.
          @return 0, if all went right.
          @see grad(UserVector<double>&, const UserVector<double>&)
          @see add_grad(UserVector<double>&, const UserVector<double>&)
          @see add_valgrad(UserVector<double>&, UserVector<double>&, const UserVector<double>&)
      */
      virtual int valgrad(double& val, UserVector<double>& g, const UserVector<double>& x) const {
        val=eval(x);
        grad(g,x);
        return 0;
      };

      /** Add's the (sub)gradient of this function to a UserVector<double>.
          @param g The UserVector<double> to add the gradient to.
          @param x The UserVector<double> to compute the gradient for.
          @see grad(UserVector<double>&, const UserVector<double>&)
          @see valgrad(double&, dvector&, const dvector&)
          @see add_valgrad(double&, dvector&, const dvector&)
      */
      virtual void add_grad(UserVector<double>& g, const UserVector<double>& x) const {
        Pointer<UserVector<double> > gr(g.getemptycopy());
        grad(*gr, x);
        g+=*gr;
      };

      /** Add's the (sub)gradient and the value of this function to a UserVector<double> and a double.
          @param val The double to add the value to.
          @param g The UserVector<double> to add the gradient to.
          @param x The UserVector<double> to compute the value and the gradient for.
          @see eval(const UserVector<double>&)
          @see add_grad(UserVector<double>&, const UserVector<double>&)
      */
      virtual void add_valgrad(double& val, UserVector<double>& g, const UserVector<double>& x) const {
        val+=eval(x);
        add_grad(g,x);
      };

      /** Gives the minimum eigenvector of the Hessian in point x.
          @param x The dvector to compute the Hessian for.
          @param param Optional parameters for eigenvalue computation.
          @return The minimum eigenvector of the Hessian.
      */
      virtual double min_eig_hess(const UserVector<double>& x, Param* param=NULL) const;

      /** Gets the convex shift parameter.
          @param lower Lower bound of the box.
          @param upper Upper bound of the box.
          @param lev The level.
          @param offset Will be added to the shift, if not convex.
          @return Convex shift parameters or 0, if function is convex.
      */
//      virtual dvector get_convex_shift(const UserVector<double>& lower, const UserVector<double>& upper, int lev, double offset) const;

#ifdef FILIB_AVAILABLE
			/** Indicates, whether the function and it's gradient can be evaluated over an interval.
			*/
			virtual bool is_interval_compliant() const { return false; }

			virtual interval<double> eval(const IntervalVector& x) const {
				IntervalVector g(x.dim());
				interval<double> val;
				valgrad(val, g, x);
				return val;
			};

			virtual void grad(IntervalVector& g, const IntervalVector& x) const {
				interval<double> val;
				valgrad(val, g, x);
			};

			virtual int valgrad(interval<double>& val, IntervalVector& g, const IntervalVector& x) const {
				out_err << "Interval evaluation not implemented for this class." << endl;
				exit(-1);
				return 1;
			}
#endif
			typedef enum { CONVEX=1, CONCAVE=2, LINEAR=CONVEX | CONCAVE, UNKNOWN=4, INDEFINITE=12 } CurvatureType;
			
			friend ostream& operator<<(ostream& out, const CurvatureType& ct);

			virtual void set_curvature(CurvatureType ct)=0;
			virtual CurvatureType get_curvature() const=0;

			static CurvatureType add_curvatures(double a1, CurvatureType ct1, double a2, CurvatureType ct2);
			static CurvatureType mult_curvature(double a, CurvatureType ct);

      /** Print's out the dimension.
          @param out The ostream to print to.
          @see dim()
      */
      virtual void print(ostream &out) const
      {  out << "Func: dim=" << dim() << endl;  };
};

// -----------------------------------------------------------------

/** Represantation of a Hessian Matrix for a function.
    So it's like a wrapper-class for Func::HessMult(...)
    @see Func::HessMult(dvector&, const dvector&, const dvector&)
*/
class HessMatrix : public UserMatrix {
  protected:
    /** The UserVector<double> to compute the Hessian for.
    */
    Pointer<const UserVector<double> > x0;

    /** The function to compute the Hessian for.
    */
    Pointer<const Func> f;

  public:
    /** Constructor for a function and a dvector.
        @param f_ The function, which Hessian should be computed.
        @param x0_ The UserVector<double> to compute the Hessian of f_ for.
    */
    HessMatrix(Pointer<const Func> f_, Pointer<const UserVector<double> > x0_)
    : UserMatrix(f_->dim()), f(f_), x0(x0_)
    { }

    HessMatrix(const Func& f_, const UserVector<double>& x0_)
    : UserMatrix(f_.dim()), f(&f_, false), x0(&x0_, false)
    { }

    /** Computes the product of the Hessian and a UserVector<double>.
        @param y The UserVector<double> to store the result in: @f$ y = \nabla^2 f(x0) \cdot x @f$.
        @param x The UserVector<double> to multiply with the Hessian.
        @see Func::HessMult(dvector&, const dvector&, const dvector&)
    */
    void MultV(UserVector<double> &y, const UserVector<double> &x) const {
      f->HessMult(y, *x0, x);
    }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using UserMatrix::MultV;
#endif

    /** Prints some information about this matrix.
        Prints the dvector x0 and the function.
        @param out The ostream to print to.
    */
    void print(ostream &out) const {
      out << "HessMatrix for UserVector<double>: " << *x0 << "function: " << *f;
    }

};

// ----------------------------------------------------------------------

/** A wrapper class to multiply a function with -1.
    *this == - *f .
*/
class MinusFunc : public Func {
  protected:
    /** The function to wrap.
    */
    Pointer<Func> f;

		SparsityInfo& get_sparsity() { return sparsity ? *sparsity : f->get_sparsity(); }

  public:
    /** Constructor for a normal func.
        @param f_ The function to wrap.
        @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
        @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
    */
    MinusFunc(Pointer<Func> f_, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
    : Func(f_ ? f_->dim() : 0, out_func_p_, out_func_log_p_), f(f_)
    { assert(f_ != NULL);
    }

		virtual bool sparsity_available() const { return sparsity || f->sparsity_available(); }

		virtual const SparsityInfo& get_sparsity() const { return sparsity ? *sparsity : ((const Func*)(Func*)f)->get_sparsity(); }

		virtual bool compute_sparsity_pattern(const vector<dvector>& sample_set);

    /** Evaluates the function for a UserVector<double>.
        @param x The UserVector<double> to evaluate.
        @return The value: -f(x)
    */
    double eval(const UserVector<double>& x) const {
      return - f->eval(x);
    }

    /** Evaluates the gradient for a UserVector<double>.
        @param g The UserVector<double> to store the result in: g = - f'(x).
        @param x The UserVector<double> to evaluate.
    */
    void grad(UserVector<double>& g, const UserVector<double>& x) const {
      f->grad(g,x);
      g*=-1;
    }
		
		dvector grad(const dvector& x) const {
			return -f->grad(x);
		}

    /** Computes the product of the Hessian in one point and a UserVector<double>.
        @param y The UserVector<double> to store the result in.
        @param x The UserVector<double> to compute the Hessian for.
        @param z The UserVector<double> to multiply the Hession with.
    */
    void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
      f->HessMult(y, x, z);
      y*=-1;
    }
		
		dvector HessMult(const dvector& x, const dvector& z) const {
			return -f->HessMult(x,z);
		}

#ifdef FILIB_AVAILABLE
		bool is_interval_compliant() const { return f->is_interval_compliant(); }

		interval<double> eval(const IntervalVector& x) const {
			return -f->eval(x);
		};

		void grad(IntervalVector& g, const IntervalVector& x) const {
			f->grad(g,x);
			g*=-1;
		};

		int valgrad(interval<double>& val, IntervalVector& g, const IntervalVector& x) const {
			int ret=f->valgrad(val, g, x);
			val*=-1;
			g*=-1;
			return ret;
		}

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Func::valgrad;
#endif

#endif

		virtual void set_curvature(CurvatureType ct) {
			switch(ct) {
				case Func::CONVEX: f->set_curvature(Func::CONCAVE); break;
				case Func::CONCAVE: f->set_curvature(Func::CONVEX); break;
				default: f->set_curvature(ct); // UNKNOWN, UNDEFINED or LINEAR
			}
		}

		virtual CurvatureType get_curvature() const {
			Func::CurvatureType ct(f->get_curvature());
			switch(ct) {
				case Func::CONVEX: return Func::CONCAVE;
				case Func::CONCAVE: return Func::CONVEX;
			}
			return ct;
		}

    /** Print's the shift-value and the shifted function.
        @param out The ostream to print to.
    */
    void print(ostream& out) const {
      out << "MinusFunc: function: " << *f;
    }

};

// ----------------------------------------------------------------------

/** A function to represent the sum of two functions with optional multiplicators.
    *this == a * f(x) + b * g(x), where a and b are double's and f and g are Func's.
    f or g can also be NULL, but not both.
*/
class SumFunc: public Func {
  private:
  	/** First function in sum.
		*/
  	Pointer<Func> f;
  	/** Second function in sum.
  	*/
  	Pointer<Func> g;
  	/** Multiplicators for the functions.
  	*/
  	double a,b;

		Func::CurvatureType curv_type;

		SparsityInfo& get_sparsity() {
			assert(sparsity_available());
			if (sparsity) return *sparsity;
			return f->get_sparsity(); // this gives me only a const-member
		}

  public:
   	/** Constructor for two functions and optional multiplicators.
   		  Asserts, that f_ or g_ is not NULL.
   		  @param f_ The first function.
   		  @param g_ The second function.
   		  @param a_ The multiplicator for the first function.
   		  @param b_ The multiplicator for the second function.
   	*/
  	SumFunc(Pointer<Func> f_, Pointer<Func> g_, double a_=1., double b_=1., Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
    : Func(f_ ? f_->dim() : (g_ ? g_->dim() : 0), out_func_p_, out_func_log_p_), f(f_), g(g_), a(a_), b(b_),
		  curv_type(add_curvatures(f_ ? a_ : 0., f_ ? f_->get_curvature() : Func::LINEAR, g_ ? b_ : 0., g_ ? g_->get_curvature() : Func::LINEAR))
    { assert(f || g);
			if (!f) { f=g; a=b; g=NULL; } // move to front
			if (f && g && f->sparsity_available() && g->sparsity_available())
				sparsity=new SparsityInfo(((const Func*)(Func*)f)->get_sparsity(), ((const Func*)(Func*)g)->get_sparsity());
		};

		virtual bool sparsity_available() const {
			if (sparsity) return true;
			if (g) return false; // f and g defined, but no sparsity set
			return f->sparsity_available();
		}

		virtual const SparsityInfo& get_sparsity() const {
			assert(sparsity_available());
			if (sparsity) return *sparsity;
			return f->get_sparsity();
		}

		virtual bool compute_sparsity_pattern(const vector<dvector>& sample_set);

		/** Evaluates the function.
			  @param x The vector to evaluate the function for.
			  @return The value of the function a * f(x) + b * g(x).
			  @see valgrad(double&, UserVector<double>&, UserVector<double>&)
		*/
		double eval(const UserVector<double>& x) const {
			return a * (f ? f->eval(x) : 0) + b * (g ? g->eval(x) : 0);
		}

		/** Computes the gradient.
			  @param y The vector to store the gradient in.
			  @param x The vector to compute the gradient for.
			  @see valgrad(double&, UserVector<double>&, UserVector<double>&)
		*/
		void grad(UserVector<double>& y, const UserVector<double>& x) const;
		
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Func::grad;
#endif

		/** Computes the value and the gradient.
			  @param val A double to store the value in.
			  @param y A UserVector to store the gradient in.
			  @param x The UserVector to compute.
			  @see eval(UserVector<double>&)
			  @see grad(UserVector<double>&, UserVector<double>&)
		*/
		int valgrad(double& val, UserVector<double>& y, const UserVector<double>& x) const;
		
		/** Computes the product of the Hessian and a UserVector.
        @param y The UserVector<double> to store the result in.
        @param x The UserVector<double> to compute the Hessian for.
        @param z The UserVector<double> to multiply the Hession with.
		*/
		void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const;
		
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Func::HessMult;
#endif

#ifdef FILIB_AVAILABLE
		bool is_interval_compliant() const { return (f ? f->is_interval_compliant() : true) && (g ? g->is_interval_compliant() : true); }

		interval<double> eval(const IntervalVector& x) const {
			return a*(f ? f->eval(x) : 0)+ b*(g ? g->eval(x) : 0);
		};

		void grad(IntervalVector& y, const IntervalVector& x) const;

		int valgrad(interval<double>& val, IntervalVector& y, const IntervalVector& x) const;
#endif

		Func::CurvatureType get_curvature() const { return curv_type; }
		void set_curvature(Func::CurvatureType ct) { curv_type=ct; }

		/** Prints some information.
				Prints a, b, f and g.
		    @param out The ostream to print to.
		*/
		void print(ostream& out) const {
			out << "SumFunc a*f(x)+b*g(x): a=" << a << " b=" << b;
			if (f) out << endl << *f; else out << " f=NULL";
			if (g) out << endl << *g; else out << " g=NULL" << endl;
		}
		
};


// ----------------------------------------------------------------------

/** Abstract base class for a separable function.
    @see Func
*/
class SepFunc : public Func {
   public:
      /** The block's of indices.
      */
      vector<ivector> block;

      /** Constructor for one block.
          @param n The size of the first block.
          @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
          @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
      */
      SepFunc(int n, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : Func(n, out_func_p_, out_func_log_p_), block(1)
      { block[0].resize(n);
        for (int i=0; i<n; i++) block[0][i]=i;
      };

      /** Constructor for a block-structure.
          @param block_ The block-structure.
          @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
          @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
          @see set_dim()
      */
      SepFunc(const vector<ivector>& block_, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : Func(0, out_func_p_, out_func_log_p_), block(block_)
      { set_dim();
      };

			/** Copy constructor.
			*/
			SepFunc(const SepFunc& f)
			: Func(f.dim(), f.out_func_p, f.out_func_log_p), block(f.block)
			{ }

      /** Sets the dimension to the correct value.
          Sets dim_ to sum of the sizes of each block.
          @see dim()
      */
      void set_dim() {
        dim_=0;
        for (int i=0; i<block.size(); i++) dim_+=block[i].size();
      }

      /** Print's out information about this function.
          Print's the dimension and the block-structure.
          @param out The ostream to print to.
          @see block
      */
      virtual void print(ostream &out) const {
        out << "SepFunc: dim=" << dim() << endl;
        for (int i=0; i<block.size(); i++)
          out << "block " << i << ": " << block[i];
      };
};

// ----------------------------------------------------------------------

/** A class for a separable QC-function.
    @f$f(x)=\sum_i (x_i*A_i*x_i+2b_i*x_i+s_i(x_i)) + c@f$.
*/
class SepQcFunc: public SepFunc {
	friend class Decomposition;
	private:
		vector<Pointer<SparsityInfo> > sparsity_block;
	public:
	class VariableIterator : public VariableIterator_Type {
		private:
			const SepQcFunc& f;
			map<int, SparsityInfo::LinearVariable>::iterator it_linear;
			map<int, SparsityInfo::NonlinearVariable>::iterator it_nonlinear;
			map<int, SparsityInfo::QuadraticVariable>::iterator it_quadratic;
			/** Over which variable type am I just iterating?
			*/
			VarType whereami;
			int blocknr;
			void find_start();

		public:
			VariableIterator(const SepQcFunc& f_, bool linear_=true, bool nonlinear_=true, bool quadratic_=false);

			virtual int bnr() const { return blocknr; }
			virtual int bindex() const;
			virtual int operator()() const { return f.block[blocknr][bindex()]; }

			virtual double coeff_lin() const;
			virtual double coeff_quad() const;

			virtual VarType type() const { return whereami; }

			virtual void operator++();

			virtual operator bool() const { return (whereami!=END); }
	};
	/** Types of functions.
	*/
	typedef enum { CONSTANT, LINEAR, QUADRATIC, NONQUAD } ftype;
	friend ostream& operator<<(ostream& out, const ftype& ft);

      /** The matrices of each block.
      */
      vector<Pointer<UserMatrix> > A;

      /** The vector's of each block.
      */
      vector<Pointer<UserVector<double> > > b;

      /** The nonlinear, nonquadratic part's of each block.
      */
      vector<Pointer<Func> > s;

			/** The curvature types of each block.
			*/
			vector<Func::CurvatureType> curv_type;

      /** The constant part.
      */
      double c;

      /** Constructor for one block.
          @param n The size of the first block.
          @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
          @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
      */
      SepQcFunc(int n, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : SepFunc(n, out_func_p_, out_func_log_p_), A(1), b(1), c(0), s(1), curv_type(1, Func::UNKNOWN), sparsity_block(1)
      { }

      /** Constructor for a block-structure.
          @param block_ The block-structure to copy.
          @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
          @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
      */
      SepQcFunc(const vector<ivector>& block_, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : SepFunc(block_, out_func_p_, out_func_log_p_), A(block_.size()), b(block_.size()), c(0), s(block_.size()), curv_type(block_.size(), Func::UNKNOWN), sparsity_block(block_.size())
      { }

      SepQcFunc(Pointer<UserMatrix> A_, Pointer<UserVector<double> > b_=NULL, Pointer<Func> s_=NULL, double c_=0, Pointer<SparsityInfo> sparsity_=NULL, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : SepFunc(A_ ? A_->dim() : b_ ? b_->dim() : s_->dim(), out_func_p_, out_func_log_p_), A(1, A_), b(1, b_), s(1, s_), c(c_),
			  curv_type(1, A_ ? Func::UNKNOWN : (s_ ? s_->get_curvature() : Func::LINEAR)), sparsity_block(1, sparsity_)
      { sparsity=sparsity_; }

			/** Copy Constructor with (optional) minus.
          @param f The SepQcFunc to copy.
          @param minus_func_ Indicates, if the function should be multiplied by -1.
      */
      SepQcFunc(const SepQcFunc &f, bool minus_func);

			/** Computes a Taylor approximation.
			*/
			SepQcFunc(const SepQcFunc& f, const UserVector<double>& point, const int degree=2);

			virtual bool sparsity_available(int k) const {
				if (sparsity_block[k]) return true;
				return false;
			}

			virtual bool sparsity_available() const { return sparsity; }

			virtual SparsityInfo& get_sparsity(int k) {
				return *sparsity_block[k];
			}

			virtual const SparsityInfo& get_sparsity(int k) const {
				return *sparsity_block[k];
			}

			virtual SparsityInfo& get_sparsity() { assert(sparsity); return *sparsity; }

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using Func::get_sparsity;
#endif

			virtual void set_sparsity(int k, Pointer<SparsityInfo> si) {
				sparsity_block[k]=si;
			}
			virtual void set_sparsity();

			bool compute_sparsity(int block_nr, const vector<dvector>& sample_set, bool replace_if_quadratic=false);

			/** The type of the a block of this function.
          @param k The block number.
          @return NONQUAD, if s[k] is set, QUADRATIC if only A[k] is set, LINEAR if only b[k] is set, CONSTANT else.
      */
      virtual ftype type(int k) {
        if (s[k]) return NONQUAD;
        if (A[k]) return QUADRATIC;
				if (b[k]) return LINEAR;
        return CONSTANT;
      }

			/** Increasing the dimension by adding a variable.
					If you add a variable to an existing block, A and s must be NULL for this block.
			    @param index Index of variable.
			    @param bnum The block, where we add the variable to.
			*/
			void add_var(int index, int bnum);

      /** Multiplies this function with a double and adds another function, multiplied with another double.
          So this function will be a*this + b*f.
          @param f The function to add.
          @param a_ The double to multiply this function with, default is 1.
          @param b_ The double to multiply the given function with, default is 1.
      */
      void addmult(const SepQcFunc& f, double a_=1., double b_=1.);

      /** Computes the value for a UserVector<double>.
          @param x The UserVector<double> to compute the value for.
          @return The value as double.
      */
      double eval(const UserVector<double>& x) const;

      /** Computes the value of one block for a UserVector<double>, without the constant part.
          @param x The UserVector<double> of k'th block to evaluate.
          @param k The block number.
          @return The value as double.
      */
      double eval(const UserVector<double>& x, int k) const {
        return (A[k] ? A[k]->xAx(x) : 0) + (b[k] ? 2* (*b[k]*x) : 0) + (s[k] ? s[k]->eval(x) : 0);
      }

      /** Computes the gradient of this function.
          @param g The UserVector<double> to store the result in.
          @param x The UserVector<double> to compute the gradient for.
      */
      void grad(UserVector<double>& g, const UserVector<double>& x) const;

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using Func::grad;
#endif

      /** Computes the gradient of one block for a UserVector<double>.
			    @param g The UserVector<double> to store the gradient in.
          @param x The UserVector<double> of the k'th block to compute.
          @param bnum The block-number.
      */
			void grad(UserVector<double>& g, const UserVector<double>& x, int bnum) const;

			dvector grad(const UserVector<double>& x, int bnum) const {
				dvector g(x.dim()); grad(g,x,bnum); return g;
			}

      /** Computes the value and the gradient for a UserVector<double>.
          @param val The double to store the value in.
          @param g The UserVector<double> to store the gradient in.
          @param x The UserVector<double> to compute the value and the gradient for.
          @return The sum of the return codes over the blocks.
      */
      int valgrad(double& val, UserVector<double>& g, const UserVector<double>& x) const;

      /** Computes the product of the Hessian and a UserVector<double>.
          @param y The UserVector<double> to store the result in.
          @param x The UserVector<double> to evaluate the Hessian for.
          @param z The UserVector<double> to multiply with the Hessian.
      */
      void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const;

      /** Compute the product of one block of the Hessian and a UserVector<double>.
          @param y The UserVector<double> to store the result in.
          @param x The UserVector<double> to evaluate the k-th block of the Hessian for.
          @param z The UserVector<double> to multiply with the Hessian.
          @param k The block number.
      */
      void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z, int k) const;

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using Func::HessMult;
#endif

#ifdef FILIB_AVAILABLE
			bool is_interval_compliant() const;

			interval<double> eval(const IntervalVector& x) const;

			interval<double> eval(const IntervalVector& x, int k) const;

			void grad(IntervalVector& g, const IntervalVector& x) const;

			void grad(IntervalVector& g, const IntervalVector& x, int k) const;

			int valgrad(interval<double>& val, IntervalVector& y, const IntervalVector& x) const;
#endif

			Func::CurvatureType get_curvature(int k) const { return curv_type[k]; }
			Func::CurvatureType get_curvature() const;
			void set_curvature(int k, Func::CurvatureType ct) { curv_type[k]=ct; }
			void set_curvature(Func::CurvatureType ct);

      /** Prints information about this function.
          Prints the dimension, constant part, block-structure and linear, quadratic and nonquadratic parts.
          @param out The ostream to print to.
      */
      void print(ostream &out) const;

      /** Prints this function in a nicer format than the other print.
          @param out The ostream to print to.
			    @param var_names The variable names.
      */
			void print(ostream& out, vector<Pointer<char> > var_names) const;

};

#endif // FUNC_H
