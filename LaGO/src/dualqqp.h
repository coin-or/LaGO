// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

#ifndef DUALQQP_H
#define DUALQQP_H

#include "standard.h"
#include "usermatrix.h"
#include "func.h"
#include "problem.h"
#include "param.h"

/** A abstract base class for a dual function of a MinlpReform.
    The dimension should be the one form the dual problem.
    You can get the dimension of the primal problem with primal_dim().

    In addition to the methods from Func, you need to implement:
    - dvector get_dual_point() to get the last dual point
    - double get_dual_val() to get the last dual value
    - dvector get_orig_point() to get the last point from the original problem
*/
class DualFunc: public Func {
   public:
      /** The MINLP problem.
      */
      MinlpProblem &qqp;

      /** The lower bounds of the dual variables.
          Should be set in a derived class.
      */
      vector<double> dual_bounds;

      /** The parameters of the Minlp-solver.
      */
      Param &param;

      /** Constructor for a MinlpReform.
          The dual dimension is set to 0 at this moment and should be set by a derived class.
          @param qqp_ The MinlpReform, this function belongs to.
          @param param_ Some parameters.
          @param out_func_p_ Pointer to an ostream to print function-related stuff to.
          @param out_func_log_p_ Pointer to an ostream to print function-related logging-stuff to.
      */
      DualFunc(MinlpProblem &qqp_, Param &param_, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : Func(0, out_func_p_, out_func_log_p_), qqp(qqp_), param(param_)
      { }

      /** The dimension of the MINLP-problem, not the dual one.
          @return The dimension of the primal problem.
      */
      int primal_dim() const { return qqp.dim(); };

      /** Gives the last dual point.
          @return The last dual point, which was used to evaluate this function.
      */
      virtual dvector get_dual_point() const=0;

      /** Gives the last dual value.
          @return The last dual value, which was calculated with this function.
      */
      virtual double get_dual_val() const=0;

      /** Gives the last solution point of the original problem.
          In some dual functions, you can get more than one point.
          @param i The number of this point to get. 0 should give you the best one.
          @return The last solution point of the given problem.
					@see get_orig_point(int)
      */
      virtual dvector get_orig_point(int i=0) const=0;

      /** The number of points, you can get with get_orig_point().
          @return 1.
      */
      virtual int nr_of_orig_points() const { return 1; }

			virtual void set_curvature(CurvatureType ct) { out_err << "DualFunc::set_curvature() not implemented. Aborting." << endl; exit (-1); };
			virtual CurvatureType get_curvature() const  { out_err << "DualFunc::get_curvature() not implemented. Aborting." << endl; exit (-1); return Func::UNKNOWN; };
};

class QqpMatrix;
class QqpExtMatrix;

/** A base class for evaluating the dual function.
    @f$ D(\lambda )=
    \sum_{k=1}^{l-1}(\rho_k^2\cdot \min\{eig_{min}(\tilde A_k(\lambda)),0\}
    +\sum_{j=r_l}^{r_{l+1}}\min\{\underline x_{j}b_{j}(\lambda ),
    \bar x_{j}b_{j}(\lambda )\})+\tilde c(\lambda ) @f$

    assignment of the dual variables:

    0:n-1 quadratic box-constraints of original variables

    n:n+block_size-1 quadratic box-constraints of the new variables

    block_size+n:block_size+n+con_size-1 remaining constraints
*/
class QqpDualFunc: public DualFunc {
    friend class QqpMatrix;
    friend class QqpExtMatrix;

		protected:
			/** The different block-types.
    			QUAD: the hessian of this block is non-zero and the linear part is zero.

			    EXTQUAD: the hessian of this block and the linear part is non-zero.

			    LIN: the hessian is zero.
			*/
			typedef enum { LIN, QUAD, EXTQUAD } b_type;

      /** Indicates for each block, which type the block of the lagrange problem has.
          QUAD, EXTQUAD or LIN.
          Initialized in init_ext.
          @see init_ext()
          @see init()
      */
      vector<b_type> block_type;

      /** The indices of the new quadratic box constraints.
      */
      ivector i_ext;

      /** The number of extended blocks.
      */
      int num_ext;

      /** The index of the start of the quadratric box constraints for the k-th block.
          Size of i_q: number of nonlinear-blocks.
      */
      ivector i_q;

      /** The number of quadratic variables.
      */
      int n_q;

      /** The mid point of the box = 0.5*(lower+upper).
          Initialized in init, should be deleted.
          @see init()
      */
      dvector mid_point;

      /** The half-length of the edges of the box.
          Initialized in init as BlockMatrix of DiagMatrix-blocks.
          @see init()
      */
      BlockMatrix W;

      /** The radius of the ball-constraints
      */
      dvector radius;

      /** The type of the ball-constraints.
          False if equality constraint and true if inequality constraint.
          Initialized in set_ball_con.
          @see set_ball_con()
          @see init()
      */
      vector<bool> ball_con;

      /** The type of scaling (currently not used).
          True if scaling and false else.
          Initialized in set_ball_con.
          @see set_ball_con()
          @see init()
      */
//      vector<bool> diag_scale;

      /** The modified b-values for the constraints.
          Size set and initialized in set_bc.
          @see set_bc()
      */
      vector<dvector> b_con;

      /** The modified b-value for the objective-function.
          Set in set_bc.
          @see set_bc()
      */
      dvector b_obj;

      /** The modified c-values for the constraints.
          Size set and initialized in set_bc.
          @see set_bc()
      */
      dvector c_con;

      /** The modified c-value for the objective-function.
          Set in set_bc.
          @see set_bc()
      */
      double c_obj;

      /** The block matrices of the modified Lagrangian
          Set in init.
          @see init()
      */
      vector<Pointer<UserMatrix> > A_lag;

      /** The linear part of the modified Lagrangian
          Set in init.
          @see set_dual_val()
      */
      mutable vector<dvector> b_lag;

      /** The constant part of the modified Lagrangian
          Set in init.
          @see set_dual_val()
      */
      mutable double c_lag;

      /** The Lagrange solution.
          @see set_subgrad()
      */
      mutable dvector x_lag;

      /** current value of the dual function.
      */
      mutable double dual_val;

      /** The subgradient.
      */
      mutable dvector subgrad;

      /** minimal eigenvalues of A_lag[k], k=0,..,block_size-1.
          eig_val[k] are the first eig_val[k].size() eigenvalues of the k'th block, starting with the smallest. */
      mutable vector<vector<double> > eig_val;

      /** eigenvectors of A_lag[k], k=0,..,block_size-1.
          eig_vec[k] are the first eig_vec[k].size() eigenvectors of the k'th block, starting with the one, belong to the smallest eigenvalue. */
      mutable vector<vector<dvector> > eig_vec;

      /** The number of quadratic-box-constraints.
          @return The number of quadratic-box-constraints.
      */
      int num_qbox() const { return n_q+num_ext; }

      /** The number of blocks.
          @return The number of blocks.
      */
      int num_blocks() const { return qqp.block.size(); };

      /** Check's, if partial Lagrange problem is extended
          Set ext[k] to true, if the bounds are not symmetric or no linear part exists.

          Or you say: Check, if modified b_k's are zero.
          @see ext
      */
      void init_ext();

      /** Initialization of the dual function.
          Allocating memory, setting counters to zero, computing u,W...

          Check's, whether the Lagrange problem is extended by calling init_ext().

          Set's the modified c,b-values by calling set_bc().
          Set's the ball-constraint type by calling set_ball_con(lower, upper).
          @see init_ext()
          @see set_bc()
          @see set_ball_con()
      */
      void init();

      /** Initialize the ball constraints.
          Initialize the ball constraints for each block with the corresponding part of lower and upper.
          Set's the type of each ball-constraint to true, if it contains one not-discrete variable.
          @see ball_con
      */
      void set_ball_con();

      /** Compute's modified b,c values for a SepQcFunc.
          Sums mid_point[k]*q.A[k]*mid_point[k] + q.b[k]*mid_point[k] over all blocks.
          Adds the constant part of q: q.c
          @param c_ The modified c-value.
          @param b_ The modified b-value.
          @param q The SepQcfunc to compute.
      */
      void get_bc(double &c_,dvector &b_,SepQcFunc &q);

      /** Initialize the modified b,c-values.
          Computes for each constraint and the objective the modified b,c-values by calling get_bc.
          @see get_bc(double&, dvector&, SepQcFunc&)
          @see b_con
          @see b_obj
          @see c_con
          @see c_obj
      */
      void set_bc();

      /** Compute the value of the dual function.
          The function value is stored in dual_val.
          @see dual_val
          @return Sum of return codes from the eigenvalue computation.
      */
      int set_dual_val() const;

      /** Evaluates a modified constraint.
          @param con_nr The number of the constraint to evaluate.
          @return The value of the constraint.
      */
      double eval_mod_con(int con_nr) const;

      /** Evaluates the modified objective.
          @return The value of the objective.
      */
      double eval_mod_obj() const;

      /** Compute a subgradient of the dual function.
          The subgradient is stored in subgrad.
          @see eval_mod_con(int)
          @see subgrad
      */
      void set_subgrad() const;

    public:
      /** number of matrix vector multiplications */
      int num_matmult;
      /** time spent by the eigenvalue calculation  */
      mutable double eig_time;

      /** The dual variables of the extended (!) dual function.
          This is the dual point, which was last used to evaluate the function.
          @see get_dual_point()
      */
      mutable dvector dual_point;

      /** Constructor for a MinlpReform and a shift.
          @param qqp_ The MinlpReform, the DualFunc consists of.
          @param param_ Parameters.
          @param out_func_p_ A pointer to an output stream for function output.
          @param out_func_log_p_ A pointer to an output stream for function logging output.
      */
      QqpDualFunc(MinlpProblem &qqp_, Param &param_, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
      : DualFunc(qqp_, param_, out_func_p_, out_func_log_p_), num_matmult(0), eig_time(0),
        eig_vec(qqp_.block.size()), eig_val(qqp_.block.size()),
        radius(qqp_.block.size()), ball_con(qqp_.block.size()),
        block_type(qqp_.block.size()), i_ext(qqp_.block.size()), i_q(qqp_.block.size()),
        c_con(qqp_.con.size()), b_obj(qqp_.dim()), b_con(qqp_.con.size()),
        A_lag(qqp.block.size()), b_lag(qqp.block.size()),
        mid_point(qqp_.dim()), W(qqp_.block)
      { init();
      };

      /** Evaluates the function value and computes the subgradient of the dual function.
          @param val A double to store the value in.
          @param g A dvector to store the subgradient in.
          @param z The dual value to compute the value and subgradient for.
          @return Return code from set_dual_val().
          @see set_dual_val()
          @see set_subgrad()
      */
      virtual int valgrad(double& val, UserVector<double>& g, const UserVector<double>& z) const;

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using DualFunc::valgrad;
#endif

      /** Evaluates this function for a dvector.
          Sets dual_point to given point and calls set_dual_val().
          This method doesn't look at the return value from set_dual_val() !
          @param z The dvector to compute the value for.
          @return The value of this function.
          @see valgrad(double&, dvector&, dvector&)
      */
      double eval(const UserVector<double>& z) const {
      	assert(z.dim() == dim());
			  dual_point=z;
			  set_dual_val();
			  return dual_val;
			}

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using DualFunc::eval;
#endif
/*
      double eval(const UserVector<double>& z) const {
        dvector z0(z);
        return eval(z0);
      }
*/
      /** Computes the gradient of this function.
          Calls valgrad().
          This method doesn't look at the return value from set_dual_val(), so use valgrad() instead !
          @param g The dvector to store the gradient in.
          @param x The dvector to compute the gradient for.
      */
/*      void grad(dvector& g, const dvector& x) const {
        double val;
        valgrad(val, g, x);
      }
*/
      void grad(UserVector<double>& g, const UserVector<double>& x) const {
				double val;
				valgrad(val, g, x);
      }
			
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using Func::grad;
#endif

      /** Does nothing.
          Should do: Compute the product of the Hessian and a dvector.
          @param y The dvector to store the result in.
          @param x The dvector to evaluate the Hessian for.
          @param z The dvector to multiply with the Hessian.
      */
//      virtual void HessMult(dvector &y, const dvector &x, const dvector &z) const { };

      /** Does nothing.
      */
      virtual void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const { };
			
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using Func::HessMult;
#endif

      /** Set the Lagrange multipliers of the quadratic box constraints.
      */
//      void set_qbox_dual_point(dvector &x,dvector &g);

      /** The last dual point of the original qqp.
          @return The dual point of the not-extended(!) QQP, used in the last iteration.
      */
      dvector get_dual_point() const;

      /** The dual point of the original qqp for an extended dual point.
          If the extended dual point (the one, you can use for eval()) is not equal to dual_point, eval(ext_dual_point) is called.
          Then get_dual_point() is called.
          @param ext_dual_point A dual point from the extended qqp.
          @see get_dual_point()
			*/
      dvector get_dual_point(dvector& ext_dual_point) const {
        if (!(ext_dual_point==dual_point)) eval(ext_dual_point);
        return get_dual_point();
      }

      /** The last dual value.
          @return The dual value, used in the last iteration.
      */
      double get_dual_val() const { return dual_val; }

      /** The dual value for an extended dual point.
          If the extended dual point (the one, you can use for eval()) is not equal to dual_point, eval(ext_dual_point) is called.
          Then get_dual_val() is called.
          @param ext_dual_point A dual point from the extended qqp.
          @see get_dual_val()
			*/
      double get_dual_val(dvector& ext_dual_point) const {
        if (!(ext_dual_point==dual_point)) eval(ext_dual_point);
        return get_dual_val();
			}

      /** Gives the solution point of the Lagrangian function.
          @param use_eig_vec The eigenvector to use to compute the lag-point.
          @return The solution point of the Lagrangian function, used in the last iteration.
      */
      dvector get_orig_point(int use_eig_vec=0) const;

      /** The dual value for an extended dual point.
          If the extended dual point (the one, you can use for eval()) is not equal to dual_point, eval(ext_dual_point) is called.
          Then get_orig_point() is called.
          @param ext_dual_point A dual point from the extended qqp.
          @see get_orig_point()
			*/
      dvector get_orig_point(dvector& ext_dual_point, int use_eig_vec=0) const {
        if (!(ext_dual_point==dual_point)) eval(ext_dual_point);
        return get_orig_point(use_eig_vec);
      }

      /** The number of eigenvectors, which are computed in each evaluation.
          @return The size of eig_vec.
      */
      int nr_of_orig_points() const { return eig_vec[0].size(); }

      /** Print's a lot of information about this function.
          @param out The ostream to print to.
      */
      void print(ostream &out) const;
};

/** A class for evaluating the dual function as in QqpDualFunc where
    the Lagrange multipliers of the quadratic box-constraints
    are eliminated by @f$\nabla L(x^\ast;\lambda)=0@f$
*/
/*
class QqpDualFunc2: public QqpDualFunc {
   public:
      /** Constructor for a MinlpReform and parameters.
          @param qqp_ The MinlpReform.
          @param param_ Parameters.
          @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
          @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
      */
/*      QqpDualFunc2(MinlpReform &qqp_,Param &param_, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p, Pointer<ostream> out_dualval_p_=NULL)
      : QqpDualFunc(qqp_,param_, out_func_p_, out_func_log_p_, out_dualval_p_)
      {  };

      /** Gives the dimension.
          Gives the dim() of upper-class less the primal dim.
      */
/*      int dim() const {
        return QqpDualFunc::dim()-primal_dim();
      };

      /** Computes the value and gradient of this function for a dvector.
          @param val The double to store the double in.
          @param g The dvector to store the gradient in.
          @param z The dvector to compute.
      */
/*      int valgrad(double &val,dvector &g,const dvector &z);
};
*/
// reduced dual function
/*
class QqpDualFunc3: public QqpDualFunc
   {
   int dual_dim; // dimension of the full dual space

   vector<int> i_dual; // indice of actual dual variables
   dvector red_dual_bounds; // lower bounds of the reduced dual func.
   bool full_dual_space; // true if the full dual space is used
   void set_full_dual_space()
   {  full_dual_space=true;dim_=dual_dim;  };
   void set_reduced_dual_space(vector<int> &iv)
   {  i_dual=iv; dim_=iv.size(); full_dual_space=false;
   red_dual_bounds.resize(dim());
   for(int i=0;i<dim();i++)
   red_dual_bounds[i]=qqp.dual_bounds[i_dual[i]];  };
   void valgrad(double &val,dvector &g,dvector &z);
   };
*/

/** One block of the block matrix of the modified QQP-Lagrangian.
*/
class QqpMatrix: public UserMatrix {
   protected:
      /** The block number. */
      int b_num;
   public:
      /** The QQP dual function. */
      QqpDualFunc &d;

      // methods
      /** Constructor for a QqpDualFunc and a block number.
          @param d_ The dual function, this matrix is based on.
          @param b_num_ The number of the block of the QqpDualFunc, this matrix should evaluate.
      */
      QqpMatrix(QqpDualFunc &d_,int b_num_)
      : UserMatrix(d_.qqp.block[b_num_].size()), d(d_), b_num(b_num_)
      { };

      /** Multiplies the matrix of the Lagrangian with a dvector.
          @param y The dvector to store the result in.
          @param x The dvector to multiply with this matrix.
      */
      void MultV(dvector &y,const dvector &x) const;

      void MultV(UserVector<double>& y, const UserVector<double>& x) const {
        dvector x0(x);
        dvector y0(y.size());
        MultV(y0,x0);
        y=y0;
      }
			
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using UserMatrix::MultV;
#endif

};

/** One extended block of the block matrix of the modified QQP-Lagrangian.
*/
class QqpExtMatrix : public UserMatrix {
   protected:
      /** The block number. */
      int b_num;
   public:
      /** The QQP dual function. */
      QqpDualFunc &d;

      /** Constructor for a QqpDualFunc and a block number.
          @param d_ The dual function, this matrix is based on.
          @param b_num_ The number of the block of the QqpDualFunc, this matrix should evaluate.
      */
      QqpExtMatrix(QqpDualFunc &d_,int b_num_)
      : UserMatrix(d_.qqp.block[b_num_].size()+1), d(d_), b_num(b_num_)
      { };

      /** Multiplies the matrix of the Lagrangian with a dvector.
          @param y The dvector to store the result in.
          @param x The dvector to multiply with this matrix.
      */
      void MultV(dvector &y, const dvector &x) const;

      void MultV(UserVector<double>& y, const UserVector<double>& x) const {
        dvector x0(x);
        dvector y0(y.size());
        MultV(y0,x0);
        y=y0;
      }
			
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
			using UserMatrix::MultV;
#endif
};

#endif // DUALQQP_H
