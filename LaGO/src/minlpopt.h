// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef MINLPOPT_H
#define MINLPOPT_H

#include "standard.h"
#include "MINLPData.h"
#include "MINLP.h"
#include "problem.h"
#include "opt.h"

#include "decomp.h"
#include "relax.h"

class MinlpOpt;
#include "linrelax.h"

/** A solution candidate.
    This is a pair of a double and a dvector, with a special compare-operator.
*/
class SolCandidate : public pair<double, dvector> {
	friend ostream& operator<<(ostream& out, const SolCandidate& s) {
		out << s.first << ": " << s.second;
		return out;
	}

	private:
		double closeval_tol;
		Pointer<dvector> diam;

	public:

		SolCandidate(double val, const dvector& x, double closeval_tol_, Pointer<dvector> diam_)
		: pair<double, dvector>(val, x), closeval_tol(closeval_tol_), diam(diam_)
		{ }

		bool operator<(const SolCandidate& y) const {
			if (first<y.first-closeval_tol*(fabs(first)+1)) return true;
			if (first>y.first+closeval_tol*(fabs(first)+1)) return false;
			for (int i=0; i<second.dim(); i++)
				if (fabs(second(i)-y.second(i))>(*diam)(i))
					return (second(i)<y.second(i));
			return false;
		}

};

/** A container for storing several reformulations. */

class Reformulation {
	private:
		MinlpOpt& opt;

		/** Adds a variable to a block of a problem.
		*/
		void add_var(vector<vector<Pointer<BlockMatrix> > >& A, vector<vector<Pointer<SepQcFunc> > >& s, MinlpProblem& prob, int k, int index, double lower, double upper, Pointer<char> name=NULL);
		/** Adds one block of a constraint as a seperated constraint, assuming, that the last added variable is the appropiate t.
		*/
		void add_con(vector<vector<Pointer<BlockMatrix> > >& A, vector<vector<Pointer<SepQcFunc> > >& s, MinlpProblem& prob, int c, int k);

	public:
		Pointer<MinlpProblem> ext_prob;
		Pointer<MinlpProblem> ext_quad_prob;
		Pointer<MinlpProblem> ext_convex_prob;

		/** To remember, which t was added for which constraint.
		    related_t should have split_prob->con.size().
				related_t[c] stores the index of the t inside it's block (!), which belongs to constraint split_prob.con[c]
		*/
		vector<int> related_t;

		Reformulation(MinlpOpt& opt_)
		: opt(opt_)
		{ }

		void reformulate();

		dvector get_short_vector(const dvector &x);
		dvector get_long_vector(const dvector &x);

		bool set_t_block(dvector& x, int k);

};

/** Updates level cuts. */

class LevelCutHandler {
	private:
		/** The problems, where we want to have a level cut for and the location of the level cut there.
		*/
		list<pair<Pointer<MinlpProblem>, int> > problems_with_levelcut;

		list<Pointer<LinearRelax> > linrelax_with_levelcut;

		double val; // the value of the level cut

	public:

		LevelCutHandler(double init_val=INFINITY)
		: val(init_val)
		{	}
		
		~LevelCutHandler();

		void add_problem(Pointer<MinlpProblem> prob);
		void add_problem(Pointer<LinearRelax> linear_relax);

		void update_level_cut(double newval);
};

/** Main class to solve a MINLP.
    @class MinlpOpt
    @param MinlpOpt mode
    %options BCP, off
		%level 0
		%default BCP
    Determines, which algorithm to use. BCP is the Branch and Cut algorithm. off stops after preprocessing.
		@param Boxreduce effort
		%options 0, 1, 2
		%default 1
		%level 2
		How much effort to spend in boxreduction.
		0 computes unknown variable bounds only, except of the IntervalReduction, which is still run on all variables.
		1 applies IntervalReduction and boxreduction by linear relaxations for all variables.
		2 uses all available techniques for all variables.
		@param Boxreduce type
		%options off, NLP, NLP2, MINLP
		%default off
		%level 0
		Which method to use for the box reduction phase II. NLP is (C)-method, NLP2 is (Cext)-method, MINLP is (Pext)-method, off, if switches this part off.
    @param Boxreduction limit for reconvexification
    %options [0,1]
    %default 0.8
		%level 1
    If the boxdiameter was reduced better than the value of this parameter, using (C), we compute a new convex relaxation from (Q).
		@param Decomposition
		%options 0, 1
		%default 1
		%level 1
		Whether to apply decomposition or not.
		@param Decomposition sample set Monte Carlo
		%options integer >= 0
		%default 20
		%level 1
		@param Decomposition sample set mid point
		%options integer >= 0
		%default 1
		%level 1
		@param Decomposition sample set vertices
		%options integer >= 0
		%default 20
		%level 1
		@param Relax check convex sample set Monte Carlo
		%options integer >= 0
		%default 200
		%level 1
		The number of sample points to use for convexity test.
		@param Check decomposition
		%options 0, 1
		%default 0
		%level 0
		Performs a test on the decomposition.
		@param Check polynomial underestimator
		%options 0, 1
		%default 0
		%level 0
		Performs a test on the polynomial underestimators.
		@param Check convexification
		%options 0, 1
		%default 0
		%level 0
		Performs a test on the convexification.
    @param heu close value tolerance
    %options double $\geq 0$
		%default 10E-8
		%level 0
    It the relative distance of two (objective) values is less than this tolerance, they are considered as equal.
    @param heu close points tolerance
    %options double $\geq 0$
		%default .0001
		%level 0
    If the maximal (over all variables) relative (to box diameter) distance is less than this tolerance, two points are considered as equal.
		@param Level Cuts
		%options 0, 1
		%default 1
		%level 1
		Whether to add a level cut.
		@param Quadratic Underestimator adaptive
		%options 0, 1
		%default 1
		%level 2
		Whether to use the new adaptive (but slower) quadratic underestimator or the old less reliable ones.
*/
class MinlpOpt : public Solver {
	friend class LinearRelax;
	friend class Reformulation;
	friend int main(int argc, char** argv);
	private:
		/** Parameters.
		*/
		Pointer<Param> param;
		
		Pointer<Timer> timer;
		
		Pointer<MINLPData> minlpdata;

		/** Indicates, whether the original problem comes from gams.
		*/
		bool is_gams_prob;

		/** The original problem.
		*/
		Pointer<MinlpProblem> orig_prob; // (P)

		/** The decomposed original problem.
		*/
		Pointer<MinlpProblem> split_prob; // (S)

		/** The nonconvex quadratic relaxation of the decomposed problem.
		*/
		Pointer<MinlpProblem> quad_prob; // (Q)

		/** The convexified quadratic decomposed problem.
		*/
		Pointer<MinlpProblem> convex_prob; // (C)

		/** Reformulated problems (Pext), (Qext), (Cext).
		*/
		Pointer<Reformulation> reform;

		vector<vector<double> > min_eigval, max_eigval;

		/** Indicates, whether the problem is convex.
		*/
		bool prob_is_convex;

		/** To remember, which constants were added to the objective during computation of (Q).
		*/
		SparseVector<double> quad_obj_c_add;
		/** To remember, which constants were added to the constraints during computation of (Q).
		*/
		vector<SparseVector<double> > quad_con_c_add;
		/** To remember, which constants were added to the objective during computation of (C).
		*/
		SparseVector<double> conv_obj_c_add;
		/** To remember, which constants were added to the constraints during computation of (C).
		*/
		vector<SparseVector<double> > conv_con_c_add;

		ivector ineq_index;

		/** Linear relaxation methods: Generates (R) and cuts.
		*/
		Pointer<LinearRelax> linear_relax;

		Pointer<LevelCutHandler> levelcuts;

		IntervalReduction intervalreduction;

		/** Computes decomposed formulation (S) of (P).
		*/
		void decompose(); // generates (S)

		/** Determines curvature.
		*/
		bool check_convex(MinlpProblem& prob);

		/** Computes the quadratic relaxation (Q) of (S).
		*/
		void quad_relax();

		/** Computes the convex relaxation (C) of (Q) or (S).
		*/
		void convex_relax();

		bool init_called;
		/** Initialization of (basic) relaxations.
		*/
		void init();
		/** Initialization of further relaxations.
		*/
		void init2();

		double sol_cand_closeval_tol;
		Pointer<dvector> sol_cand_diam;

		enum { box_off, box_C, box_Cext, box_Pext } boxreduce_type;
		int boxreduce_effort;
		double boxreduce_time;

		/** The variables, where still no bound is known.
		*/
		list<int> unbounded_var;

		/** Gives a problem, which containts only the convex constraints from a given problem.
		*/
		Pointer<MinlpProblem> get_convex_prob(Pointer<MinlpProblem> prob);
		/** Gives a problem, which containts only the convex constraints from a given problem. Uses another problem as reference to get curvature information.
		*/
		Pointer<MinlpProblem> get_convex_prob(Pointer<MinlpProblem> prob, Pointer<MinlpProblem> prob_curv_ref);
		double print_box_reduce_quality(dvector& oldlow, dvector& oldup, Pointer<MinlpProblem> prob, char* prefix);

	public:
		/** Solution candidates.
		*/
		set<SolCandidate> sol_cand;

		/** A solution of (C_ext).
		*/
		Pointer<dvector> sol_Cext;
		bool sol_Cext_is_solution;

		double low_bound;

		MinlpOpt(Pointer<MinlpProblem> prob, Pointer<Param> param_, bool is_gams_prob_=false, Pointer<ostream> out_solver_p=out_out_p, Pointer<ostream> out_solver_log_p=out_log_p);

		void box_reduce0();

		/** Reduces the box using BoundsFinder.
				Using NLP method with (C), if available, or convex constraints from (Q) or (P).
		*/
		void box_reduce1();

		double box_reduce2();

		void box_reduce3();

		/** Starts Rounding and Partition Heuristic.
		*/
		int start_round_part_heu();

		/** Starts Branch-and-Bound Algorithm.
		*/
		int start_bb();

		/** Checks, if the initial point is inside it's bound.
		    Exits, if the initial point of the original problem isn't inside the box.
		*/
		void check_initial_point();

		/** Tries to solve the problem.
		*/
		int solve();

		using Solver::solve;
};

#endif
