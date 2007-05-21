// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

#ifndef PROBLEM_H
#define PROBLEM_H

#include "standard.h"
#include "usermatrix.h"
#include "func.h"

/** Class for defining a mixed-integer nonlinear minimization problem.
    If i_discr is empty the problem is continuous.
    @class MinlpProblem
    @param constraint scaling eps
		%options $\geq 0$
		%default 0
    The value $\epsilon$ to which constraints are going to be scaled. Not applied currently.
    @param constraint scaling default
		%options W*gradient | condition Jac | adjust by box ends | equalized
		%default no scaling
    How to scale the constraints.
    @param constraint scaling max patterns
    %options integer $\geq 0$
		Maximum number of different patterns to scale specific type of constraints.
    @param constraint scaling pattern <$n$>
		%options W*gradient | condition Jac | adjust by box ends | equalized | double $\geq 0$
    How to scale the constraints that match the patters.
    Example: for all the constraints that contain the letter `E'
    \texttt{constraint scaling pattern 2 : E 1.}
    \texttt{constraint scaling pattern 2 : E number or equalized or condition Jac or none}
*/
class MinlpProblem {
   /** The output-operator.
       Calls print(out).
       @param out The ostream to print to.
       @param a The MinlpProblem to print.
       @see print(ostream&)
   */
   friend ostream& operator << (ostream &out, const MinlpProblem& a) {
     a.print(out); return out;
   };

   double compute_scale(char* option, int c, double eps, UserVector<double>& x);

   protected:
      /** The dimension.
      */
      int dim_;

   public:
      /** A Pointer to an ostream to print problem-relevant information to.
          @see out_problem
      */
      Pointer<ostream> out_problem_p;
      /** A Pointer to an ostream to print problem-relevant logging-information to.
          @see out_problem_log
      */
      Pointer<ostream> out_problem_log_p;
      /** The dereferenced Pointer<ostream> out_problem_p.
          If out_problem_p==NULL, no output appears.
          Only use this in the form "out_problem << ...", not as argument !
          @see MinlpProblem::out_problem_p
      */
#define out_problem if (out_problem_p) (*out_problem_p)
      /** The dereferenced Pointer<ostream> out_problem_log_p.
          If out_problem_log_p==NULL, no output appears.
          Only use this in the form "out_problem_log << ...", not as argument !
          @see MinlpProblem::out_problem_log_p
      */
#define out_problem_log if (out_problem_log_p) (*out_problem_log_p)

			/** Types of optimization problems.
			*/
			typedef enum { QQP, MINLP } ptype;

      /** The type of the problem.
          QIP, QQP or MINLP
      */
      ptype problem_type;

      /** A name for the problem.
          @see MinlpProblem(ptype, char*, ostream*, ostream*)
      */
      Pointer<char> prob_name;

      /** The block structure of the variables.
          @see add_var(int, int, bool, double, double, char*)
      */
      vector<ivector> block;

      /** The names of the variables.
          @see add_var(int, int, bool, double, double, char*)
      */
      vector<Pointer<char> > var_names;

      // definition of the optimization problem
      /** The objective function.
          @see add_obj(SepQcFunc*)
      */
      Pointer<SepQcFunc> obj;

      /** The constraints functions.
          @see add_con(SepQcFunc&, bool, char*)
      */
      vector<Pointer<SepQcFunc> > con;

      /** The names of the constraints.
          @see add_con(SepQcFunc&, bool, char*)
      */
      vector<Pointer<char> > con_names;

      /** Indicates, which constraints are equations and which are inequation.
          @see add_con(SepQcFunc&, bool, char*)
      */
      vector<bool> con_eq;

      /** The indices of the discrete variables.
          @see add_var(int, int, bool, double, double, char*)
      */
      vector<int> i_discr;

      /** The indices of the continuous variables.
				  @see add_var(int, int, bool, double, double, char*)
      */
      vector<int> i_cont;

      /** Indicates, whether variable is discrete (true) or continuous (false).
          @see add_var(int, int, bool, double, double, char*)
      */
      vector<bool> discr;

      /** The lower bounds of the variables.
          @see add_var(int, int, bool, double, double, char*)
      */
      dvector lower;

      /** The upper bounds of the variables.
          @see add_var(int, int, bool, double, double, char*)
      */
      dvector upper;

      /** Used for fixing the discrete variables.
      */
      dvector primal_point;

      /** (Standard-)Constructor for an optimization problem.
          Variables, objective function and constraints should then be add with the appropriate methods.
          @param p_type The type of the optimization problem. By default, it's set to MINLP.
          @param prob_name A name for the problem.
          @param out_problem_p_ A Pointer to an ostream to print the problem-releated information to, default is out_out_p.
          @param out_problem_log_p_ A Pointer to an ostream to print the problem-related logging-information to, default is out_log_p.
          @see add_var(int, int, bool, double, double, char*)
          @see add_obj(SepQcFunc*)
          @see add_con(SepQcFunc&, bool, char*)
      */
      MinlpProblem(ptype p_type=MINLP, char* prob_name_=(char*)NULL, Pointer<ostream> out_problem_p_=out_out_p, Pointer<ostream> out_problem_log_p_=out_log_p)
      : problem_type(p_type), prob_name(prob_name_ ? strdup(prob_name_) : NULL), dim_(0), out_problem_p(out_problem_p_), out_problem_log_p(out_problem_log_p_)
      { };

      /** Constructor for an objective function, the lower and upper bounds and a boolean, which indicates, whether all variables are discrete or continuous.
          Set's the dimension of the problem to the dimension of the objective function.
          @param obj_ The objective function.
          @param lower_ The lower bounds of the variables.
          @param upper_ The upper bounds of the variables.
          @param discr_ If true, each variable is set to be discrete. If false, each variable is set to be continuous. Default is false.
          @param out_problem_p_ A pointer to an ostream to print the problem-releated information to, default is out_out_p.
          @param out_problem_log_p_ A pointer to an ostream to print the problem-related logging-information to, default is out_log_p.
      */
      MinlpProblem(Pointer<SepQcFunc> obj_, const dvector &lower_, const dvector &upper_, bool discr_=false, Pointer<ostream> out_problem_p_=out_out_p, Pointer<ostream> out_problem_log_p_=out_log_p)
      : dim_(obj_->dim()), prob_name(0), problem_type(MINLP), out_problem_p(out_problem_p_), out_problem_log_p(out_problem_log_p_), block(obj_->block), discr(obj_->dim(), discr_), lower(lower_), upper(upper_), var_names(obj_->dim()), primal_point(obj_->dim())
      { if (discr_) for (int i=0; i<dim(); i++) i_discr.push_back(i);
				else for (int i=0; i<dim(); i++) i_cont.push_back(i);
        add_obj(obj_);
      }

      /** The copy constructor.
          @param p The MinlpProblem to copy.
      */
      MinlpProblem(const MinlpProblem &p)
      : dim_(p.dim()), block(p.block), var_names(p.var_names),
        obj(p.obj), con(p.con), con_names(p.con_names), con_eq(p.con_eq),
        i_discr(p.i_discr),i_cont(p.i_cont), discr(p.discr),
        lower(p.lower), upper(p.upper), primal_point(p.primal_point),
        problem_type(p.problem_type), prob_name(p.prob_name),
        out_problem_p(p.out_problem_p), out_problem_log_p(p.out_problem_log_p)
      { }

			/** Copy constructor for the variable structure of one block from another MinlpProblem.
			    This does not copy the constraints.
					@param p The problem to take one block from.
					@param k The number of the block to take.
			*/
			MinlpProblem(const MinlpProblem& p, int k);

      /** Virtual Destructor.
      */
      virtual ~MinlpProblem() { }

      /** Get's the dimension of the problem.
          @return The dimension of the problem.
      */
      int dim() const { return dim_; };

      /** Add's a variable.
          Increases the dimension of the problem.

          Add's the index i_ to block bnum.

          If discr_ is true, adds the index to i_discr.
          If discr_ is false, adds the index to i_cont.

          Add's discr_ to disrc.

          Add's lower_ to lower.
          Add's upper_ to upper.
          @param i_ The index of the new variable.
          @param bnum The number of the block, where this variable belongs to.
          @param discr_ Indicates, whether variable is discrete (true) or not (false).
          @param lower_ The lower bound of this variable.
          @param upper_ The upper bound of this variable.
          @param name The name of the variable as pointer of char, default is NULL.
          @see add_obj(SepQcFunc*)
          @see add_con(SepQcFunc&, bool, char*)
      */
      virtual void add_var(int i_, int bnum, bool discr_, double lower_, double upper_, char* name=NULL);

      /** Adds a constraint to the problem.
          Adds f to con, init's the Ablock of f, copys the name and add eq to con_eq.
          @param f The function.
          @param eq Indicates, whether it's an equality-(true) or inequality-(false)-function. Default-Value is true.
          @param name The name of the constraint as pointer of char. Default is NULL.
          @see add_var(int, int, bool, double, double, char*)
          @see add_obj(SepQcFunc&)
          @see del_con(int)
      */
      virtual void add_con(Pointer<SepQcFunc> f, bool eq=true, char* name=NULL);

      /** Deletes one constraint from the problem.
          @param connr The number of the constraint to delete.
          @see add_con(SepQcFunc&, bool, char*)
      */
      virtual void del_con(int connr);

      /** Add's (or replace) the objective function.
          @param f The objective function, which should be used.
          @see add_var(int, int, bool, double, double, char*)
          @see add_con(SepQcFunc&, bool, char*)
      */
      virtual void add_obj(Pointer<SepQcFunc> f) { obj=f; };

      /** Does a second taylor approximation for this problem.
          Changes the type of this problem to QQP.
          @param point The point, where the approximation should be done.
      */
      void taylor_approx(UserVector<double>& point, const int degree=2);

      /** Checks the point for feasibility according to the constraints, not the lower and upper bounds.
          @param x The UserVector<double> to check.
          @param tol The tolerance, default is rtol.
          @return The number of not sattisfied constraints.
      */
      int feasible(const UserVector<double>& x, double tol, ostream* out=NULL);
      
      /** Prints the names of those constraints that are most violated by a given point.
       * Each constraint is scaled by max(1,||gradient of constraint at x||_2).  
       * @param x The point to check
       * @param out Where to print to.
       * @param nr How many violated constraints to print.
       * @param tol A feasibility tolerance. 
       */
      void print_most_violated_constraints(const UserVector<double>& x, ostream& out, int nr=5, double tol=1E-4);

      /** Scales the constraints, using the primal_point.
          @param param Optional parameters, default is NULL.
          @return The number of scaled constraints.
          @see scale(UserVector<double>&)
      */
      int scale(Pointer<Param> param=NULL) { return scale(primal_point, param); }

      /** Scales the constraints.
          @param param Optional parameters, default is NULL.
          @param x The UserVector<double> to compute for.
          @return The number of scaled constraints.
      */
      int scale(UserVector<double>& x, Pointer<Param> param=NULL);

			void get_sparsity(vector<Pointer<SparsityInfo> >& si) const;
			Pointer<SparsityInfo> get_sparsity(int k) const;

      /** Prints some information about the problem.
          Prints out:
          - The problem type, dimension and name.
          - Information about the objective function.
          - The block structure.
          - The names of the block-variables.
          - The indices of the discrete and continuous variables.
          - Information about the constraint-functions.
          - The lower-bounds of the dual variables.

          @param out The ostream to print to.
          @see obj
          @see block
          @see con
          @see i_discr
          @see i_cont
          @see discr
          @see con_eq
          @see lower
          @see upper
          @see problem_type
          @see prob_name
      */
      virtual void print(ostream& out) const;

			/** Prints the problem in gams format. Working only for MIQQPs!
			*/
			void print_as_gams(ostream& out) const;
};


/** A Penalty Function for a MinlpProblem.
    The function has the form @f$ \Psi (x) = \sum_{i\in I_i} \delta_i \max \{0,f_i(x)\} + \sum_{i\in I_e} \delta_i |f_i(x)|. @f$
*/
class MinlpPenaltyFunc: public Func {
  private:
    /** The MinlpProblem.
    */
    Pointer<MinlpProblem> minlp;

  public:
    /** The vector of delta-values.
    */
    dvector delta;

    /** Constructor for a MinlpProblem.
        Sets the delta-values to 1000.
        @param minlp_ The MinlpProblem, which PenaltyFunc should be represented.
        @param penalty_param The delta-values for penalty of each constraints as parameters, default is "penaltyfunc.res".
    */
    MinlpPenaltyFunc(Pointer<MinlpProblem> minlp_, Pointer<Param> penalty_param=NULL)
    : Func(minlp_->dim()), minlp(minlp_), delta(minlp_->con.size(), 1000.)
    { if (penalty_param) {
    		penalty_param->read();
	      for (int i=0; i<minlp->con.size(); i++)
  	      delta[i]=penalty_param->get_d(minlp->con_names[i], 1000.);
      }
    }

    /** Evaluates the penalty function for a UserVector<double>.
        @param x The UserVector<double> to compute the value for.
        @return The value of this function.
        @see MinlpPenaltyFunc
        @see valgrad(double&, UserVector<double>&, const UserVector<double>&)
    */
    double eval(const UserVector<double>& x) const;
		using Func::eval;

    /** Computes the value and the gradient of this penalty function.
        The gradient of the absolut of an equality constraint is handled as zero, if the gradient of the equality constraint is zero.
        @param val A double to store the value in.
        @param y A UserVector<double> to store the gradient in.
        @param x The UserVector<double> to compute.
        @return 0
        @see eval(const UserVector<double>&)
        @see grad(UserVector<double>&, const UserVector<double>&)
    */
    int valgrad(double& val, UserVector<double>& y, const UserVector<double>& x) const;
		using Func::valgrad;

    /** Computes the gradient of this penalty function.
        Calls valgrad(double&, UserVector<double>&, const UserVector<double>&).
        @param y The UserVector<double> to store the gradient in.
        @param x The UserVector<double> to compute.
        @see valgrad(double&, UserVector<double>&, const UserVector<double>&)
    */
    void grad(UserVector<double>& y, const UserVector<double>& x) const;
		using Func::grad;

    /** Computes the product of the hessian of this matrix and a UserVector<double>.
        @param y The UserVector<double> to store the result in.
        @param z The UserVector<double> to compute the hessian for.
        @param x The UserVector<double> to multiply with the hessian.
    */
    void HessMult(UserVector<double>& y, const UserVector<double>& z, const UserVector<double>& x) const;
		using Func::HessMult;

		virtual void set_curvature(CurvatureType ct) { out_err << "MinlpPenaltyFunc::set_curvature() not implemented. Aborting." << endl; exit (-1); };
		virtual CurvatureType get_curvature() const  { out_err << "MinlpPenaltyFunc::get_curvature() not implemented. Aborting." << endl; exit (-1); return Func::UNKNOWN; };

    /** Prints some information about this function.
        Prints the delta vector.
        @param out The ostream to print to.
        @see delta
    */
    void print(ostream& out) const {
      out << "MinlpPenaltyFunc: dim=" << dim() << endl;
      out << "delta-vector: " << delta;
    }

};


/** Container for a MIP problem.
    To represent a problem of the form min c^T x+d s.t. \lb b<=Ax<=\ub b, x\in [\lb x, \ub x], x_i\in Z, i\in J.
		Optimized for fast access by a MIPSolver. Adding or removing columns or constraints can be very slow.
*/
class MipProblem {
	public:
		typedef enum { CONTINUOUS, BINARY, INTEGER } VarType;

		friend ostream& operator<<(ostream& out, const VarType& type) {
			switch (type) {
				case CONTINUOUS: out << "Continuous"; break;
				case BINARY: out << "Binary"; break;
				case INTEGER: out << "Integer"; break;
				default: out << "Unknown !!";
			}
			return out;
		}
		friend ostream& operator<<(ostream& out, const MipProblem& mip);

	private:
		int nr_col, nr_row;

	// Variables
		dvector col_lower, col_upper;
		vector<Pointer<char> > col_names;
		vector<VarType> col_types;

	// Objective function
		Pointer<UserVector<double> > obj;
		double obj_const;

	// Constraints
		SparseMatrix A;
		dvector row_lower, row_upper;
		vector<Pointer<char> > row_names;

	public:
		/** Constructs an empty MIP.
		    nr_col continuous variables are constructed, which are all unbounded.
				nr_row constraints of the form 0<=0*x<=0 are constructed.
				The objective function is zero.
				@param nr_col The number of columns (variables).
				@param nr_row The number of rows (constraints).
		*/
		MipProblem(int nr_col=0, int nr_row=0);

		/** Constructs an MIP for a given column structure with no rows and no objective.
		    @param lower The lower bounds of the columns.
				@param upper The upper bounds of the columns.
				@param names The names for the columns, if available. Default is a vector of size 0.
		*/
		MipProblem(const dvector& lower, const dvector& upper, const vector<int>& i_discr=vector<int>(0), const vector<Pointer<char> >& names=vector<Pointer<char> >(0));

		/** Transforms a MINLP to a MIP, if the MINLP is a MIP in fact.
		*/
		MipProblem(const MinlpProblem& minlp);

		/** Removes all rows from the MIP.
		*/
		void reset();

		int dim() const { return nr_col; }
		int rows() const { return nr_row; }

		const SparseMatrix& getMatrix() const { return A; }
		const dvector& getRowLower() const { return row_lower; }
		const dvector& getRowUpper() const { return row_upper; }
		const dvector& getColLower() const { return col_lower; }
		const dvector& getColUpper() const { return col_upper; }
		const VarType& getColType(int i) const { assert(i<dim() && i>=0); return col_types[i]; }
		const UserVector<double>& getObj() const { assert(obj); return *obj; }
		double getObjConst() const { return obj_const; }

		void setColBounds(int index, double low, double up);
		void setColType(int index, VarType type);
		void setColName(int index, Pointer<char> name);
		void setCol(int index, double low, double up, VarType type=CONTINUOUS, Pointer<char> name=NULL);

		void setObj(Pointer<UserVector<double> > obj_, double obj_const_=0.);

		void setRow(int index, const UserVector<double>& row, double low, double up, Pointer<char> name=NULL);
		void setRow(int index, const UserVector<double>& row, const ivector& indices, double low, double up, Pointer<char> name=NULL);
		void setRowBounds(int index, double low, double up);

		void resize_col(int nr_col_);
		void resize_row(int nr_row_);

		void finish();
};

#endif // PROBLEM_H
