// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

#ifndef OPT_H
#define OPT_H

#include "standard.h"
#include "usermatrix.h"
#include "func.h"
#include "param.h"
#include "dualqqp.h"

/** Abstract base class for solving an optimization problem.
    You need to implement the following method:
    - int solve()
*/
class Solver {
      /** Output-Operator.
          Calls a.print(out).
          @param out The ostream to print to.
          @param a The Solver to print.
          @return The ostream.
      */
      friend ostream& operator << (ostream& out, const Solver& a)
      {  a.print(out); return out; };

   protected:
      /** The number of variables.
          @see dim()
      */
      int dim_;

      /** The time, needed to find the solution.
          @see time()
      */
      double time_;

      /** The optimal value.
          @see opt_val()
      */
      double opt_val_;

      /** The number of iterations.
          @see iter()
      */
      int iter_;
   public:
      /** A Pointer to an ostream to print solver-relevant information to.
          @see out_solver
      */
      Pointer<ostream> out_solver_p;
      /** A Pointer to an ostream to print solver-relevant logging-information to.
          @see out_solver_log
      */
      Pointer<ostream> out_solver_log_p;
      /** The dereferenced Pointer<ostream> out_solver_p.
          If out_solver_p==NULL, no output appears.
          Only use this in the form "out_solver << ...", not as argument !
          @see Solver::out_solver_p
      */
#define out_solver if (out_solver_p) (*out_solver_p)
      /** The dereferenced Pointer<ostream> out_solver_log_p.
          If out_solver_log_p==NULL, no output appears.
          Only use this in the form "out_solver_log << ...", not as argument !
          @see Solver::out_solver_log_p
      */
#define out_solver_log if (out_solver_log_p) (*out_solver_log_p)

      /** A tolerance value.
          @see set_tol(double)
      */
      double tol;

      /** The maximal number of iterations.
          @see set_itermax(int)
      */
      int iter_max;

      /** The solution point.
      */
      dvector sol_point;

      /** Constructor for the dimension.
          Sets sol_point to zero, time_ and tol to -1, iter_max to INF, optimal and no_sol to false.
          @param n The dimension of the problem.
          @param out_solver_p_ A Pointer to an ostream, to print solver-output to, default is out_out_p.
          @param out_solver_log_p_ A Pointer to an ostream, to print solver-related logging-output to, default is out_log_p.
      */
      Solver(int n, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p)
      : dim_(n), sol_point(n), iter_(0), time_(-1), tol(-1.), opt_val_(INFINITY), iter_max(INF),
        out_solver_p(out_solver_p_), out_solver_log_p(out_solver_log_p_)
      {  };

      /** Virtual Destructor.
      */
      virtual ~Solver() { };

      /** Gives the dimension.
          @return The dimension of the problem.
      */
      int dim() const { return dim_;  };

      /** Gives the optimal value.
          @return The optimal value.
      */
      double opt_val() { return opt_val_;  };

      /** Gives the time, needed to solve the problem.
      */
      double time() { return time_;  };

      /** Gives the number of iterations, needed to solve the problem.
      */
      int iter() { return iter_;  };

      /** Solves the problem.
          Abstract.
          @return A status code: 0, if all went right.
          @see solve(dvector&)
      */
      virtual int solve()=0;

      /** Solves the problem for a starting point.
          Sets sol_point to x.
          Calls solve().
          @param x The dvector to start the solver with.
          @see solve()
      */
      virtual int solve(dvector &x) { sol_point=x; return solve(); };

      /** Print's some information about the problem.
          Print's the dimension, the tolerance, the solution time, the number of iterations and if an optimal solution was found.
          @param out The ostream to print to.
      */
      virtual void print(ostream &out) const {
				out << "dim: " << dim() << " tol: " << tol << endl;
				out << "solution time: " << time_ << endl;
				out << "iterations: " << iter_ << endl;
			};
};


/** Abstract base class for an optimizer for solving a dual problem of the form: max { obj(x) | x_i>=lower_bound_i }.
*/
class DualSolver : public Solver {
  private:
    /** Is the threshold control enabled ?
    */
    bool threshold_cntrl;

    /** Is the convergence rate control enabled ?
    */
    bool conv_rate_cntrl;
    /** If conv_rate_cntrl is set and improvement in last minor_iter iterations is less than stopping_rho * rel_imp1, check() breakes the solving process.
    */
    double stopping_rho;
    /** The number of minor iterations for the convergence rate control.
    */
    int minor_iter;
    /** The value in the last major iteration.
    */
    double last_major_val;
    /** The value in the last iteration.
    */
    double last_val;
    /** The value in the first iteration.
    */
    double first_major_val;
    /** First relative improvement.
    */
    double max_rel_improvement;
    /** Counter for the number of iterations with improvements (serious steps).
    */
    int improve_iter;

  protected:
    /** Parameters.
    */
    Param& param;

    /** Dual function.
    */
    DualFunc& obj;

    /** Logs the present points.
        If it's time for a log entry, the dual_value, dual_point and the points from the original problem are stored.
        @see log_frequency
        @see dual_vals
        @see dual_points
        @see orig_points
		*/
    void do_log();

    /** Checks the actual iteration.
        Calls do_log().
        @param val The last value of the dual function.
				@return 0, if the solver should continue
				@return 100, if threshold is reached
				@return 10, if convergence rate is too low
				@see do_log()
    */
    int check(double val);

  public:
    /** Last relative improvement, if conv_rate_cntrl is enabled.
    */
    double last_rel_improvement;

    /** Threshold for the dual value.
    */
    double threshold;

    /** Lower bounds of the dual points.
    */
    dvector lower_bound;

    /** To log the values of the dual function.
    */
    Pointer<vector<double> > dual_vals;
    /** To log the points of the dual problem.
    */
    Pointer<vector<dvector> > dual_points;
    /** To log the points of the original problem.
    */
    Pointer<vector<vector<dvector> > > orig_points;
    /** The log frequency.
        Indicates, that every log_frequency'th point should be stored.
        0 (default) for no logging.
    */
    int log_frequency;
    /** Should we add the primal points to the sample set.
        Default is true.
    */
    bool store_primals;

    /** Constructor for a function and lower bounds.
        Sets threshold to INFINITY.
        @param f The function to maximize.
        @param lower_bound_ The lower bounds of the dual points.
        @param threshold_ Threshold.
        @param out_solver_p_ A Pointer to an ostream to print solver-relevant information to.
        @param out_solver_log_p_ A Pointer to an ostream to print solver-relevant logging-information to.
    */
    DualSolver(DualFunc& f, vector<double> lower_bound_, Param& param_, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p)
    : Solver(f.dim(), out_solver_p_, out_solver_log_p_), threshold(INFINITY), lower_bound(f.dim()), obj(f),
      log_frequency(0), store_primals(true), param(param_), last_major_val(-INFINITY), last_val(-INFINITY), improve_iter(-1)
    { for (int i=0; i<dim(); i++) lower_bound[i]=lower_bound_[i];
      param.read();
      threshold_cntrl=param.get_i("control threshold", 1);
      conv_rate_cntrl=param.get_i("control convergence rate", 1);
			stopping_rho=param.get_d("stopping rho", 0.1);
			minor_iter=param_.get_i("minor iterations", 5);
		}

		virtual ~DualSolver() { }

		double rel_improvement1() {
		  return (last_major_val-first_major_val)/fabs(first_major_val);
		}

		double rel_improvement2() {
		  return last_rel_improvement;
		}

		double rel_improvement_max() {
		  return max_rel_improvement;
		}

		int serious_steps() {
		  return improve_iter;
		}

};

/** Abstract class for a local optimizer.
    A normal Solver with an additional method to get the Lagrangian multipliers.
*/
class LocOpt : public Solver {
	protected:
		Timer timer;
	public:
		static bool nlp_solver_available();
		/** Gives a local optimizer, if available.
		*/
		static Pointer<LocOpt> get_solver(const Pointer<MinlpProblem> prob, Pointer<Param> param, char* param_prefix=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		/** Gives a LP solver, if available.
		*/
		static Pointer<LocOpt> get_lp_solver(const Pointer<MinlpProblem> prob, Pointer<Param> param, char* param_prefix=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		/** Gives a local optimizer for the original problem, if available.
	      The given problem should be the original one, only the variable bounds can be different.
		*/
		static Pointer<LocOpt> get_solver_origprob(const Pointer<MinlpProblem> prob, Pointer<Param> param, char* param_prefix=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		LocOpt(int n, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p)
		: Solver(n, out_solver_p_, out_solver_log_p_)
		{ }

		virtual ~LocOpt() { }

		virtual double time() { return timer; }

		virtual void reinit() { };

		virtual dvector get_lag_multipliers()=0;
};

class SimpleCut;


/** A general solver for mixed-integer linear programs (MIPs).
*/
class MIPSolver {
	public:
		/** The solution statuses.
		    - SOLVED problem is solved and optimum found.
				- FEASIBLE Feasible point found, but not proven optimal.
				- UNBOUNDED Problem is unbounded.
				- INFEASIBLE Problem is infeasible.
				- ITERATIONLIMITEXCEEDED The iteration limit is exceeded.
				- ABORTED Solve process aborted for some other reason.
		*/
		typedef enum { SOLVED, FEASIBLE, UNBOUNDED, INFEASIBLE, ITERATIONLIMITEXCEEDED, ABORTED, UNKNOWN } SolutionStatus;

		friend ostream& operator<<(ostream& out, const SolutionStatus status) {
			switch(status) {
				case SOLVED: out << "Solved"; break;
				case FEASIBLE: out << "Feasible"; break;
				case UNBOUNDED: out << "Unbounded"; break;
				case INFEASIBLE: out << "Infeasible"; break;
				case ITERATIONLIMITEXCEEDED: out << "Iterationlimit exceeded"; break;
				case ABORTED: out << "Aborted"; break;
				case UNKNOWN: out << "Unknown"; break;
				default: out << "SolutionStatus not known!!";
			}
			return out;
		}

		/** Used to identify a row in the LP.
		    Specialized by MIPSolver's.
				Shouldn't be copied.
		*/
		class RowItem {
			public:
				RowItem() { }
				virtual ~RowItem() { }
		};

		/** Used to identify a column in the LP.
		    Specialized by MIPSolver's.
				Shouldn't be copied.
		*/
		class ColItem {
			public:
				ColItem() { }
				virtual ~ColItem() { }
		};

		static Pointer<MIPSolver> get_solver(const MipProblem& mip, Pointer<Param> param);

		virtual ~MIPSolver() { }

		virtual int nr_col()=0;
		virtual int nr_row()=0;

		virtual void set_tol(double tol)=0;
		virtual void set_maxiter(int maxiter)=0;

		/** Next solve will be a cold start. All added rows and columns will be deleted.
		    Modified columns and rows are NOT set to their original values!
		*/
		virtual void reset()=0;

		/** Apply warm start if possible.
		*/
		virtual SolutionStatus solve()=0;

		/** Cold start from primal starting point.
		*/
		virtual SolutionStatus solve(const UserVector<double>& x) { return solve(); };

		virtual SolutionStatus solveMIP()=0;

		/** Checks the problem for feasiblity.
		    This method can never return the SolutionStatus SOLVED.
		*/
		virtual SolutionStatus feasible()=0;

		/** The last point, we found.
		    If last SolutionStatus was SOLVED, the solution. If it was FEASIBLE, a feasible point.
				Sets x to the primal values of the first x.dim() columns.
		*/
		virtual void get_primal(UserVector<double>& x)=0;
		virtual dvector get_primal() { dvector x(nr_col()); get_primal(x); return x; }
		/** Gives the primal value for a specific column.
		*/
		virtual double get_primal(const ColItem& colitem)=0;

		virtual int get_colindex(const ColItem& colitem)=0;

		/** Gives the dual variables for the rows.
				Sets mu to the dual values of the first mu.dim() rows.
		    @param mu To store the duals.
		*/
		virtual void get_dual(UserVector<double>& mu)=0;
		/** Gives the dual for a specific row.
		*/
		virtual double get_dual(const RowItem& rowitem)=0;

		/** Gives the reduced costs for the columns.
		    Sets rc to the reduced costs of the first rc.dim() columns.
		    @param rc To store the reduced costs.
		*/
		virtual void get_reducedcosts(UserVector<double>& rc)=0;
		/** Gives the reduced costs for a specific column.
		*/
		virtual double get_reducedcosts(const ColItem& colitem)=0;
		
		/** Gives the row activities for the rows.
		    @param rowact To store the row activities of the first rowact.dim() rows.
		*/
		virtual void get_rowactivity(UserVector<double>& rowact)=0;
		/** Gives the row activities for a specific row.
		*/
		virtual double get_rowactivity(const RowItem& rowitem)=0;
		

		/** The last optimal value.
		*/
		virtual double get_optval()=0;

		virtual int get_iter()=0;

		virtual void set_obj(const UserVector<double>& obj, double obj_const=0.)=0;
		virtual void modify_obj(int i, double coeff)=0;

		/** Adds a constraint to the problem.
		    @param rowitem A RowItem, where we can store the row identification.
		    @param row The row.
        @param lb The lower bounds.
				@param ub The upper bounds.
		*/
		virtual const RowItem* add_row(const UserVector<double>& row, double low, double up)=0;
		virtual const RowItem* add_row(const UserVector<double>& row, const ivector& indices, double low, double up) {
			dvector bigrow(nr_col()); bigrow.set_block(row, indices);
			return add_row(bigrow, low, up);
		};
		virtual void add_rows(const vector<Pointer<UserVector<double> > >& rows, const dvector& low, const dvector& up) {
			for (int i=0; i<rows.size(); i++) add_row(*rows[i], low[i], up[i]);
		}
		virtual void add_rows(const vector<pair<dvector, ivector> >& rows, const dvector& low, const dvector& up) {
			for (int i=0; i<rows.size(); i++) add_row(rows[i].first, rows[i].second, low[i], up[i]);
		} 
		virtual void add_rows(list<const RowItem*>& rowitems, const vector<pair<dvector, ivector> >& rows, const dvector& low, const dvector& up) {
			for (int i=0; i<rows.size(); i++) rowitems.push_back(add_row(rows[i].first, rows[i].second, low[i], up[i]));
		}
		virtual void delete_row(const RowItem& rowitem)=0;
		virtual void delete_rows(const list<const RowItem*>& rowitems) {
			for (list<const RowItem*>::const_iterator it(rowitems.begin()); it!=rowitems.end(); ++it) delete_row(**it);
		}
		virtual void modify_row(const RowItem& rowitem, double low, double up)=0;
		virtual void modify_row(int index, double low, double up)=0;

		virtual const ColItem* add_col(double low, double up, MipProblem::VarType=MipProblem::CONTINUOUS)=0;
		virtual const ColItem* add_col(const UserVector<double>& col, double obj_coeff, double low, double up, MipProblem::VarType=MipProblem::CONTINUOUS)=0;
		/** Adds several columns.
		    @param colitems A list, where we can add the ColItem's to.
				@param low Lower bounds for new columns. -INFINITY values in low are replaced by the -\infty representation of the MIPSolver.
				@param up Upper bounds for new columns. INFINITY values in up are replaced by the \infty representation of the MIPSolver.
		*/
		virtual void add_cols(list<const ColItem*>& colitems, const dvector& low, const dvector& up)=0;
		virtual void add_cols(list<const ColItem*>& colitems, vector<Pointer<UserVector<double> > >& cols, const vector<double>& obj_coeff, const dvector& low, const dvector& up)=0;
		virtual void delete_col(const ColItem& colitem)=0;
		virtual void delete_cols(const list<const ColItem*>& colitems) {
			for (list<const ColItem*>::const_iterator it(colitems.begin()); it!=colitems.end(); ++it) delete_col(**it);
		}
		virtual void modify_col(const ColItem& colitem, double low, double up, MipProblem::VarType vartype)=0;
		virtual void modify_col(const ColItem& colitem, const UserVector<double>& col, double obj_coeff, double low, double up, MipProblem::VarType vartype)=0;
		/** If you are sure, which index your variable is.
		*/
		virtual void modify_col(int index, double low, double up, MipProblem::VarType type)=0;
		
		virtual double get_collow(int index)=0;
		virtual double get_colup(int index)=0;
		
		virtual int generate_cuts(list<Pointer<SimpleCut> >& rowcuts)=0;
};

/** Wrapper class to use a MIPSolver to solve LP's, which are given as MINLP's.
*/
class LPSolver : public LocOpt {
	private:
		Pointer<MIPSolver> mipsolver;
		const Pointer<MinlpProblem> prob;
//	Pointer<Param> param;
//		char* param_prefix;

		void initmip();

		// called from solve or solve(x)
		int solved(MIPSolver::SolutionStatus);

	public:
		LPSolver(const Pointer<MinlpProblem> prob_/*, Pointer<Param> param_, char* param_prefix_=NULL*/);

		dvector get_lag_multipliers();

		void reinit() { initmip(); }

		int solve() { return solved(mipsolver->solve()); }
		int solve(dvector& start) { return solved(mipsolver->solve(start)); }
};

/** Abstract class for a solver, which minimize a function over a box.
    Trys to solve @f$min_{x\in(lower, upper)} f(x)@f$.
*/
class BoxMinimizer : public Solver {
	protected:
	  /** The function to minimize.
  	*/
	  const Func& f;

	  /** The box constraints.
  	*/
	  Pointer<dvector> lower, upper;

	public:
	  BoxMinimizer(const Func& f_, Pointer<dvector> lower_, Pointer<dvector> upper_)
	  : Solver(f_.dim()), f(f_), lower(lower_), upper(upper_)
	  { }

	  virtual void set_box(Pointer<dvector> lower_, Pointer<dvector> upper_) {
	  	lower=lower_;
	  	upper=upper_;
	  }


};

class BoxLocOpt : public BoxMinimizer {
	private:
		Pointer<LocOpt> locopt;
		Pointer<MinlpProblem> prob;
		Pointer<Param> param;

	public:
		BoxLocOpt(Func &f_, Pointer<dvector> lower_, Pointer<dvector> upper_, Pointer<Param> param_=NULL, char* param_prefix=NULL)
		: BoxMinimizer(f_, lower_, upper_), param(param_),
		  prob(new MinlpProblem(new SepQcFunc(NULL, NULL, Pointer<Func>(&f_, false)), lower_ ? *lower_ : dvector(f_.dim()), upper_ ? *upper_ : dvector(f_.dim())))
		{ locopt=LocOpt::get_solver(prob, param, param_prefix, NULL, NULL);
		}

		BoxLocOpt(const SepQcFunc& f_, Pointer<dvector> lower_, Pointer<dvector> upper_, Pointer<Param> param_=NULL, char* param_prefix=NULL)
		: BoxMinimizer(f_, lower_, upper_), param(param_),
			prob(new MinlpProblem(Pointer<SepQcFunc>((SepQcFunc*)&f_, false), lower_ ? *lower_ : dvector(f_.dim()), upper_ ? *upper_ : dvector(f_.dim())))
		{ locopt=LocOpt::get_solver(prob, param, param_prefix, NULL, NULL);
		}

		BoxLocOpt(SepQcFunc& f_, Pointer<dvector> lower_, Pointer<dvector> upper_, Pointer<Param> param_=NULL, char* param_prefix=NULL)
		: BoxMinimizer(f_, lower_, upper_), param(param_),
			prob(new MinlpProblem(Pointer<SepQcFunc>(&f_, false), lower_ ? *lower_ : dvector(f_.dim()), upper_ ? *upper_ : dvector(f_.dim())))
		{ locopt=LocOpt::get_solver(prob, param, param_prefix, NULL, NULL);
		}

		void set_box(Pointer<dvector> lower_, Pointer<dvector> upper_) {
			if (lower==lower_ && upper==upper_) return; // same bounds, so no change needed
      BoxMinimizer::set_box(lower_, upper_);
      prob->lower=*lower;
      prob->upper=*upper;
      locopt->iter_max=iter_max;
      locopt->reinit();
		}

		int solve() {
			locopt->iter_max=iter_max;
			int ret=locopt->solve();
			sol_point=locopt->sol_point;
			opt_val_=locopt->opt_val();
			return ret;
		}

		int solve(dvector& x) {
			locopt->iter_max=iter_max;
			int ret=locopt->solve(x);
			sol_point=locopt->sol_point;
			opt_val_=locopt->opt_val();
			return ret;
		}

};

#endif // OPT_H
