// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef BCP_H
#define BCP_H

#include "standard.h"
#include "problem.h"
#include "func.h"
#include "opt.h"
#include "relaxopt.h"
#include "cuts.h"
#include "boxfind.h"

class ExtremePoint : public dvector {
	public:
  	SparseVector<double> rmpcolumn;
		double rmpobjcoeff;
		const MIPSolver::ColItem* rmpcolitem;

		ExtremePoint(const dvector& point)
		: dvector(point), rmpcolitem(NULL)
		{ }

		ExtremePoint(const dvector& point, Pointer<LinearRelax> linear_relax, int k)
		: dvector(point), rmpcolitem(NULL)
		{ set_rmpdata(linear_relax, k);
		}

		void set_rmpdata(Pointer<LinearRelax> linear_relax, int k);
};

#include "node.h"

class ColumnGenerator;
class LagHeu;

/** A branch-cut-and-price algorithm for solving a
    nonconvex MINLP problem.
    @class MinlpBCP
    @param MinlpBCP max iter
    %options integer $\geq 0$
		%default 10000
		%level 2
    The maximum number of Branch and Bound iterations.
		@param MinlpBCP max time
		%options double $\geq 0$
		%default 3600
		%level 2
		The maximum amount of seconds, which can be used by the preprocessing and MinlpBCP. 
		@param MinlpBCP gap tol
		%options double $\geq 0$
		%default 0.01
		%level 2
		Gap tolerance. Stops, if gap between upper and lower bound is smaller than gap tol.
		@param Lagrangian cuts
		%options 0, 1
		%default 1
		%level 0
		Indicates, whether to use Lagrangian cuts. If you use RMP bounds, you would be well advised to let them switched on.
    @param BCP bound type
    %options NLP,RMP,LP,LP-RMP,stop
		%default LP
		%level 1
    Determines, which bounding method to use. Options other then LP or NLP are likely to fail.
		\begin{itemize}
		\item NLP: Using the convex relaxation (Cext).
		\item LP: Using the linear relaxation (R).
		\item RMP: Does Lagrangian decomposition based on a $\mu$, computed by an inner approximation.
		\item LP-RMP: Using LP bounds but constructs and solves also the inner approximation.
		\item stop: After the first BCP phase (preprocessing, RMP bounds), subdivision and bound computation is stopped, only upper bounds for the remaining nodes are computed.
		\end{itemize}
		@param BCP preprocess max iter
		%options integer $\geq 0$
		%default 0
		%level 0
		The maximum number of BCP preprocessing iterations.
    @param BCP subdiv typ
    %options Binary, Cost, Bisection, Violation
		%default Binary
		%level 2
    The branching method. First, binary subdivision is tried. If all binaries are fixed, further actions depend on the value of this parameter. 
		\begin{itemize}
		\item Binary: If all binary variables are fixed, no further subdivision is performed.
		\item Cost: Tries to subdivide w.r.t. a variable, for which a maximum improvement of the Lagrangian can be achieved. (not tested)
		\item Bisection: Subdivides w.r.t. a variables, which boxdiameter is maximal.
		\item Violation: Tries to subdivide w.r.t. a variable, which can minimize the violation of the reference point.
    If it fails to find a variable, bisection is used.
		\end{itemize}
		@param IntervalGradient cuts
    %options 0, 1
		%default 0
		%level 1
    Enables (relaxed) IntervalGradientCuts.
		@param MIP cuts
		%options 0, 1
		%default 1
		%level 1
		Indicates, whether to derive cuts from the LP relaxation by considering the binary restrictions in the original problem.
		Currently, MixedIntegerRoundingCuts from the Cgl are used.
		@param maxcut
		%options 0, 1
		%default 0
		Indicates, whether the problem is a MaxCut problem.
		@param BCP upper bound effort
		%options $\{0, 1, 2, 3\}$
		%default 0 (2 for MaxCut)
		%level 0
		How much effort to spend in computation of upper bounds.
		\begin{itemize}
		\item[$\geq 0$] performs just local optimization, starting from the reference point.
		\item[$\geq 1$] applies preswitching as well.
		\end{itemize}
		@param LagHeu
		%options first, second, second b, Simulated Annealing
		%default none
		%level 0
		Which Lagrangian Heuristic to use to compute upper bounds. Only available, when RMP bounds are used.
		"first" means the first one, we implemented, "second" is the third one, also called LagHeu2, "second b" is a modification of LagHeu2.
		@param stopping rho
		%options $\geq 0$
		%default 0.1
		For the convergence rate control.
		If the relative improvement over the last iterations falls under $\rho$ times the first relative improvement, the convergence rate is considered as too small.
		@param minor iterations
		%options integer
		%default 5
		How many iterations to consider for the computation of one relative improvement.
		@param BCP IntervalReduction
		%level 1
		%options 0, 1
		%default 1
		If we should apply boxreduction based on interval arithmetic after branching.
		@param Memory limit
		%level 2
		%options $\geq 0$
		%default 0
		The amount of totally allocated memory (swaped and non-swaped) in Megabytes, LaGO is allowed to use. If set to 0, no limit is used.
		When the limit is exceeded in the Branch and Cut, nodes from the branching tree are pruned until LaGO consumes less than 95% of the limit on memory or only one node is left in the tree. 
*/
class MinlpBCP : public RelaxationSolver {
	friend class ColumnGenerator;
	friend class MinlpNode;
	private:
		Pointer<dvector> sol_C;
		bool sol_C_is_solution;
		
		/** Indicated whether original problem is convex.
		 * Default: false.
		*/
		bool prob_is_convex;

		/** The original problem for each block (P_k), without objective function.
		*/
		vector<Pointer<MinlpProblem> > block_prob;
		/** The convexification for each block (C_k), without objective function.
		*/
		vector<Pointer<MinlpProblem> > block_convex_prob;

		/** A pool of extreme points and columns for each block.
		*/
		vector<list<ExtremePoint> > ExtremePoints;

		vector<Pointer<MinlpProblem> > block_sub_convex_prob;

		Pointer<MinlpProblem> sub_convex_prob;

		/** The Lagrangian Subproblems for each block.
		*/
		vector<Pointer<MinlpProblem> > lag_problem;

		Pointer<ColumnGenerator> colgen;

		/** A second Branch-and-Bound tree */
		multimap<double, Pointer<MinlpNode> > bb_tree;

		/** The current node.
		    Needed by add_sol_candidate.
		*/
		Pointer<MinlpNode> current_node;

		/** gap tolerance */
		double gap_tol;

		/** bound improvement tolerance */
		double bound_impr_tol;

		Pointer<ostream> bound_print;

		/** starts the BB-algorithm
		@param bb_tree A Branch-and-Bound tree
		*/
		double start_bb();
		
		/** Indicates, whether we use IntervalGradientCuts.
		*/
		bool intgrad_cuts;
		IntervalGradientCutGenerator intgrad_cutgen;
		
		LinearizedConCutGenerator linconcutgen;
		
		/** Indicates, whether we want to derive cuts from the LP relaxation.
		*/
		bool mip_cuts;

		bool add_sol_candidate(const dvector& x);

		/** Finds solution candidates using a heuristic defined by heu_type.
		    @param node The MINLP node.
		    @param sol_point The solution point.
		*/
		int find_sol_candidates(Pointer<MinlpNode> node);

		int preswitching(Pointer<MinlpNode> node);

		Pointer<LagHeu> lagheu;

		/** Initializes the BCP-algorithm */
		void init();

		/** Initializes the block problems. */
		void init_block_problems();

		void clean_sub_problems();

		/** Adds a point to the RMP point pool
		    @param w The RMP point.
		    @param k The index of the block.
				@return The iterator to the RMP point in the pool and true, if it was a new point. False, else.
		*/
		pair<list<ExtremePoint>::iterator, bool> add_ExtremePoint(const dvector &w, int k, Pointer<MinlpNode> = NULL);

		/** Initializes the RMP of the root node.
		*/
		int init_ExtremePoints(Pointer<MinlpNode> node);

		int primal_init_ExtremePoints(const dvector& x, Pointer<MinlpNode> node);

		//$\mu$ is a dual solution point of (C\ext) or (R)
		int dual_init_ExtremePoints(Pointer<MinlpNode> node);

		/** Removes RMP points, which are not used by any node anymore.
		*/
		void prune_ExtremePoints(multimap<double, Pointer<MinlpNode> >& bb_tree);

		void project_ExtremePoints(dvector& x, Pointer<MinlpNode> node);

		// heu_type
		// bound_type
		typedef enum { RMP_bound, NLP_bound, LP_bound, LP_RMP_bound, stop_bound } t_bound_type;
		t_bound_type pre_bound_type;
		t_bound_type maj_bound_type;
		/** The current bound type.
		*/
		t_bound_type bound_type;
		// lagsolve_type[k]
		enum { BranchCut } lagsolve_type; //BranchCut should be BCP
		// subdiv_type
		enum { BinSubdiv, CostSubdivLag, CostSubdivNewton, BisectSubdiv, RMPSubdiv, ViolSubdiv } subdiv_type;

		int upper_bound_effort_level;

		/** Do a preprocessing Branch & Bound.
		*/
		int pre_bb_max_iter;

		/** Indicates, whether we should add Lagrangian cuts.
		*/
		bool lag_cuts;

		/** The number of computations of NLP-bounds or RMP-bounds, which failed, while (R[U]) was solved.
		*/
		int bound_failed;
		int bound_computed;
		double bound_time;
		double init_RMP_time;

		/** The number of lagrangian subproblems, we (tried to) solved.
		*/
		int lagprob_solves;

		double find_solcand_time;
		int nr_solcand_found;

		double subdiv_time;

		int update_ExtremePoints_count;

		Pointer<Timer> timer;
		double max_time;

		bool is_maxcut;

		Pointer<IntervalReduction> intervalreduction;

		/** Reducing the box and updating the relaxation after subdivision.
		*/
		bool boxreduce(Pointer<MinlpNode> node, int index, IntervalReduction::which_bound_type which_bound);
		
		/** Checks by interval arithmetic whether the box in node is feasible by evaluation the constraints over the box.
		 */ 
		bool feasibility_check(Pointer<MinlpNode> node);

		// --------------------- bounding

		/** Computes a lower bound of the objective value according to the bound_type: RMP_bound, NLP_bound, LP_bound.
		    If a lower bound can be computed , low_bound is updated and ref_point is set to the appropiate point. And dual_point is set.
		    @return 0, if a bound was computed; other, if a bound could not be computed.
		*/
		int set_low_bound(Pointer<MinlpNode> node);

		/** Improves a lower bound without subdivision.
		    @return If the subproblem is still feasible, and $(\lb v_{new}(U)-\lb v(U))/(1+|\lb v(U)|)$.
		*/
		pair<bool, double> improve_bound(Pointer<MinlpNode> node);

		int set_NLP_bound(Pointer<MinlpNode> node, bool improve=false);
		int set_RMP_bound(Pointer<MinlpNode> node);
		int set_LP_bound(Pointer<MinlpNode> node);
		int set_LP_RMP_bound(Pointer<MinlpNode> node);
		int improve_LP_bound(Pointer<MinlpNode> node);

		/** Improves a lower bound after subdivision
		    @param k The block number of the branching variable.
				@param i The block index of the branching variable.
		    @param node The new node.
				@param 1, if Node is infeasible, 0 else.
		*/
		int update_subdiv_bound(int k, int i, Pointer<MinlpNode> node);

		// --------------- branching

		/** Subdivides the feasible set into several sets according to the subdiv_type:
		    BinSubdiv, CostSubdiv
		    @param nodes A list, where we can add the new nodes to.
				@param node The node to subdivide.
				@return The index of the subdivision variable.
		*/
		int subdivide(list<Pointer<MinlpNode> >& nodes, Pointer<MinlpNode> node);

		/** subdivision according to violated binary constraints */
		int bin_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node);

		/** subdivision according to pseudo costs I*/
		void cost_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node);

		/** subdivision according to pseudo costs II*/
//		void cost_subdiv2(list<Pointer<MinlpNode> >& nodes, Pointer<MinlpNode> node);

		/** subdivision at the midpoint of the largest edge */
		void bisect_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node);

		/** subdivision according to constraint violation */
		void viol_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node);

		/** branching at variable i of block k w.r.t the cut-value*/
		void rect_subdiv(list<Pointer<MinlpNode> >& nodes, Pointer<MinlpNode> node, int k, int i, double cut);

		/** number of branchings on continuous variables.
		*/
		int nr_subdiv_contvar;
		/** number of subdivision by bisection.
		*/
		int nr_subdiv_bisect;

	  // -----------  lagrangian subproblems

		/** Builds the Lagrangian sub-problems with empty objective function.
		*/
		void init_lag_problems(Pointer<MinlpNode> node);

		/** Updates the objective functions of the lagrangian subproblems.
		    @param dual_point The new dual point to use.
		*/
		void update_lag_problems(const dvector& dual_point);
		void update_lag_problem(int k, const dvector& a);

		typedef struct LagSolveStatus_s {
			set<SolCandidate> solset;
			int ret;
			int iter;
			double lowbound;
			double value;
			int new_points;
		} LagSolveStatus;

		/** Solves the k-th Lagrangian sub-problem and
		    adds Lagrangian cuts
		    @param solset A set, where we can store the solutions, we found, in.
		    @param k The block number.
				@param node An optional node to store the new cuts in.
				@param temp_cut Temporary cut.
		    @return The solver status.
		*/
		void solve_lag_problem(LagSolveStatus& status, int k, Pointer<MinlpNode> node, Pointer<SepQcFunc> temp_cut=NULL);

		// --------------------------- Convergence check

    /** If conv_rate_cntrl is set and improvement in last minor_iter iterations is less than stopping_rho * rel_imp1, check() breakes the solving process.
    */
    double conv_rate_cntrl_stopping_rho;
    /** The number of minor iterations for the convergence rate control.
    */
    int conv_rate_cntrl_minor_iter;
    /** The value in the last major iteration.
    */
    double conv_rate_cntrl_last_major_val;
    /** The value in the last iteration.
    */
    double conv_rate_cntrl_last_val;
    /** The value in the first iteration.
    */
    double conv_rate_cntrl_first_major_val;
    /** First relative improvement.
    */
    double conv_rate_cntrl_max_rel_improvement;
    /** Counter for the number of iterations with improvements (serious steps).
    */
    int conv_rate_cntrl_improve_iter;

    /** Checks the actual iteration.
        @param val The last value of the dual function.
				@return 0, if the solver should continue
				@return 1, if convergence rate is too low
    */
    int conv_rate_check(double val);

		// ---------------------- Memory control

		int mem_limit;
		void mem_check();

	public:
		MinlpBCP(Pointer<MinlpProblem> orig_prob_, Pointer<MinlpProblem> split_prob_, Pointer<LinearRelax> linear_relax_,
			bool is_gams_prob=false, double closeval_tol_=0., Pointer<dvector> diam_=NULL, Pointer<Param> param_=NULL,
			Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		virtual ~MinlpBCP();

		void set_convex_relax(Pointer<MinlpProblem> convex_prob_, Pointer<dvector> sol_C_, bool sol_C_is_solution_);

		void set_reform(Pointer<Reformulation> reform_, Pointer<dvector> sol_C_, bool sol_C_is_solution_);

		void set_reform(Pointer<Reformulation> reform_);

		void set_linear_relax(Pointer<LinearRelax> linrelax_);
		
		void set_problem_is_convex(bool prob_is_convex_) { prob_is_convex=prob_is_convex_; }

		void set_MINLPData(Pointer<MINLPData> minlpdata_);
		
		void set_timer(Pointer<Timer> timer_);

		/** Calls the branch-cut-and-price algorithm.
		*/
		int solve();

		int solve(dvector& start);

};

#endif
