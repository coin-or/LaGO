// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef LAGHEU_H
#define LAGHEU_H

#include "standard.h"
#include "problem.h"
#include "relaxopt.h"
#include "linrelax.h"
#include "bcp.h"

/** Base class for Lagrangian Heuristics.
    @class LagHeu
		@param LagHeu penalty parameter
		%options double $\geq 0$
		%default 10
		Determines penalty parameter $\delta$ in LagHeu-LP and Project-LP.
*/
class LagHeu : public RelaxationSolver {
	private:
		/** LP to project a point onto the feasible set of (R) with respect to the 1-norm.
		*/
		Pointer<MIPSolver> project_LP;
		/** c^T x from linear relaxation.
		*/
		Pointer<SparseVector<double> > project_LP_obj;

		list<const MIPSolver::RowItem*> project_LP_cuts;

		/** The index of the start of the x-s=\hat x constraints.
		*/
		int x_con_start;

		Pointer<MIPSolver> lagheu_LP;
		int z_start;

		bool project_and_round_rek(dvector& x, map<int, double>& unfixed, int& switch_count);

		int switch_count_limit;

	protected:
		Pointer<MinlpNode> node;

		vector<list<const MIPSolver::ColItem*> > lagheu_LP_columns;

		/** Initialize the project LP.
		    project LP: Minimizing objective function and distance from a given point x, s.t. the coupling constraints are feasible.
				The variable bounds are initialized to [-\infty, \infty] and will be by update_project_LP.
		*/
		void make_project_LP();
		/** Sets variable bounds according to current node in the project LP.
		*/
		void update_project_LP();

		/** Sets reference point x in the project LP.
		*/
		void update_project_LP_x(const dvector& x);

		/** Solves the current project LP.
		    @param sol To store the solution.
				@return Return code of LP solver.
		*/
		MIPSolver::SolutionStatus solve_project_LP(dvector& sol);
		MIPSolver::SolutionStatus solve_project_LP(double& val);

		void project_LP_fix_variable(int i, double val);
		void project_LP_unfix_variable(int i);

		/** solves an LP to get a combination of RMP-points that
		    violates the coupling constraints as less as possible.
		    @param z The solution of the LP.
		    @return The solver status.
		*/
		void make_lagheu_LP(double delta);

		MIPSolver::SolutionStatus solve_lagheu_LP();
		MIPSolver::SolutionStatus solve_lagheu_LP(dvector &z, dvector& x);
		MIPSolver::SolutionStatus solve_lagheu_LP(vector<dvector> &z, dvector& x);

		void lagheu_LP_fixblock(int k, const MIPSolver::ColItem* colitem);
		void lagheu_LP_unfixblock(int k);

		/** Computes the violation of the coupling constraints.
		*/
		double couple_con_violation(const dvector& x);

		/** Penalty parameter in project LP
		*/
		double delta;

		/** successive projection and rounding of a given point to find a point, which is close to x, minimizes the objective function, satisfies the binary-constraints and is feasible for the linear relaxation.
		    @param x The point to try. The result is written to x.
				@return True, if successfull, false else.
		*/
		bool project_and_round(dvector& x);

	public:
		LagHeu(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob=false, double closeval_tol_=0., Pointer<dvector> diam_=NULL, Pointer<Param> param_=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p)
		: RelaxationSolver(orig_prob_, is_gams_prob, closeval_tol_, diam_, param_, out_solver_p_, out_solver_log_p_),
		  delta(param_->get_d("LagHeu penalty parameter", 10.)), switch_count_limit(2*orig_prob_->i_discr.size())
		{	set_linear_relax(linrelax_);
		}

		void set_linear_relax(Pointer<LinearRelax> linrelax_) {
			RelaxationSolver::set_linear_relax(linrelax_);
			make_project_LP();
		}

		virtual int solve(Pointer<MinlpNode> node_)=0;

		int solve() {
			assert(node);
			return solve(node);
		}

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif
};

class LagHeu1 : public LagHeu {

	public:
		LagHeu1(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob=false, double closeval_tol_=0., Pointer<dvector> diam_=NULL, Pointer<Param> param_=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p)
		: LagHeu(orig_prob_, linrelax_, is_gams_prob, closeval_tol_, diam_, param_, out_solver_p_, out_solver_log_p_)
		{ }

		int solve(Pointer<MinlpNode> node_);

#if (!defined(__GNUC__)) // icc
		using LagHeu::solve;
#endif
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif
};

/** Lagrangian Heuristic using Simulated Annealing.
    @class LagHeu_SimAnnealing
		@param Simulated Annealing iter max
		%options integer $\geq 0$
		%default MAX(20, number of blocks)
		%level 1
		The number of major iterations of the Simulated Annealing walk. So this is the number of local optimizations, we perform.
		@param Simulated Annealing weights
		%options distance, violation
		%default violation
		Which weigths to use for the states in the sim. annealing walk.
		violation defines the weight as the objective value plus the violation of the coupling constraints.
		distance defines the weight as a minimum distance to a point, which fullfills the coupling constraints and minimizes the objective.
*/
class LagHeu_SimAnnealing : public LagHeu {
	private:
		/** RMP points, stored in a better accessible structure.
		*/
		vector<vector<list<ExtremePoint>::iterator> > Z;

		/** Indices of vectors in Z, which have more than one RMP points.
		*/
		vector<int> nonsingles;

		/** Maximal number of minor iterations.
		*/
		int minor_iter_max;

		enum { VIOLATION, DISTANCE } weight_type;

		void get_Wz(dvector& Wz, const ivector& z);

		double weight(const ivector& z, const dvector& Wz);

		/** Performs a walk.
				@param z Input: Starting point. Output: Last point.
				@param weight Input: Weight of starting point. Output: Weight of last point.
				@param T temperature.
				@param l Iterations.
		*/
		void walk(ivector& z, dvector& Wz, double& weight_z, double T, int max_iter);

	public:
		LagHeu_SimAnnealing(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob=false, double closeval_tol_=0., Pointer<dvector> diam_=NULL, Pointer<Param> param_=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		int solve(Pointer<MinlpNode> node_);

#if (!defined(__GNUC__)) // icc
		using LagHeu::solve;
#endif
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif
};

/** Another lagrange heuristic.
    @class LagHeu2
    @param LagHeu2 max locopt
		%options integer $\geq 0$
		%default 100
		%level 1
		Number of best candidates to keep for local optimization.
    @param LagHeu2 max combinations
		%options integer $\geq 0$
		%default 10000
		%level 1
		Upper bound on the number of combinations of extreme points to explore.
*/
class LagHeu2 : public LagHeu {
	private:
		dvector x;

		/** { ((key, k), { (z*_kj, w_kj) }_j }_k
		*/
		map<pair<double, int>, multimap<double, list<ExtremePoint>::iterator> > W;

		Pointer<MinlpPenaltyFunc> penalty;

		map<double, dvector> candidates;
		int max_candidates;

		int max_ns;

		double ns_done;
		double ns_all;
		int lastprint;

		void search(map<pair<double, int>, multimap<double, list<ExtremePoint>::iterator> >::reverse_iterator it_k, double ns);

		/** Tries to make x feasible.
		    @param res The vector to store the result in.
				@return True, if project_LP was feasible, false else.
		*/
		bool project(dvector& res);

		/** Gives key for x.
		*/
		double key(const dvector& x);

	public:
		LagHeu2(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob=false, double closeval_tol_=0., Pointer<dvector> diam_=NULL, Pointer<Param> param_=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		int solve(Pointer<MinlpNode> node_);

		void set_reform(Pointer<Reformulation> reform_);

#if (!defined(__GNUC__)) // icc
		using LagHeu::solve;
#endif
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif
};

/** Modification of LagHeu2
    @class LagHeu2b
    @param LagHeu2 max locopt
		%options integer $\geq 0$
		%default 100
		%level 1
		Number of best candidates to keep for local optimization.
    @param LagHeu2 max combinations
		%options integer $\geq 0$
		%default 10000
		%level 1
		Upper bound on the number of combinations of extreme points to explore.
*/
class LagHeu2b : public LagHeu {
	private:
		dvector x;
		vector<dvector> z;

		/** { ((key, k), { (z*_kj, w_kj) }_j }_k
		*/
		map<pair<double, int>, multimap<double, list<ExtremePoint>::iterator> > W;

		dvector W_val;
		list<int> unfixed_blocks;

		Pointer<MinlpPenaltyFunc> penalty;

		multimap<double, dvector> candidates;
		int max_candidates;

		int max_ns;

		double ns_done;
		double ns_all;
		int lastprint;

		int block_value(double& val, int k);

		void search(double ns);

		/** Gives key for x.
		*/
		double key();

		/** Tries to make x feasible.
		*/
		bool project();

		bool switch_binaries(double val, bool eq, map<int,double>& coeff);
		void make_binaries_feasible();

		void switch_binaries2(double val, bool eq, map<int,double>& coeff);
		void switch_binaries2_rek(double val, bool eq, double& bestval, dvector& pt, map<int,double>::iterator coeff_it, map<int,double>& coeff);
		void move_continuous(double& val, bool eq, map<int,double>& cont_coeff);		
	public:
		LagHeu2b(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob=false, double closeval_tol_=0., Pointer<dvector> diam_=NULL, Pointer<Param> param_=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		int solve(Pointer<MinlpNode> node_);

		void set_reform(Pointer<Reformulation> reform_);

#if (!defined(__GNUC__)) // icc
		using LagHeu::solve;
#endif
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif
};

#endif
