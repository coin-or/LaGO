// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef LINRELAX_H
#define LINRELAX_H

#include "standard.h"
#include "problem.h"
#include "cuts.h"
#include "boxfind.h"

class MinlpNode;
class LinearRelaxSolver;
class LinearRelaxSolverGeneral;
class LinearRelaxSolverMIP;

/** A linear relaxation.
    @class LinearRelax
	@param deep cuts
	%options 0 or 1
	%default 1
	%level 1
    Indicates, whether to add deep cuts (1) or not (0).
	@param Level Cuts
	%options 0 or 1
	%default 1
	%level 1
	Indicates, whether to add a level cut (1) or not (0).
    @param max cuts nr
    %options integer $\geq 0$
    %default 500000
	%level 1
    The maximum number of cuts handled (over all nodes).
    @param Cuts inactive time limit global
    %options integer $\geq 0$
    %default 10
    %level 0
    The number of iterations a global (i.e., valid for all nodes) cut should be inactive before it is deleted.
    @param Cuts inactive time limit local
    %options integer $\geq 0$
    %default 3
    %level 0
    The number of iterations a local (i.e., valid for only one node and its successors) cut should be inactive before it is deleted.
*/
class LinearRelax {
	friend class LinearRelaxSolverGeneral;
	friend class LinearRelaxSolverMIP;
	public:
	/** Class to represent a linear (block) constraint.
	*/
	class LinConstraint {
		public:
			vector<Pointer<UserVector<double> > > b;
			double c;
			/** Indicates, whether an equality constraint is ment.
			*/
			bool eq;
			Pointer<char> name;

			LinConstraint(vector<Pointer<UserVector<double> > >& b_, double c_, bool eq_, Pointer<char> name_=NULL)
			: b(b_), c(c_), eq(eq_), name(name_)
			{ }

			LinConstraint(Pointer<UserVector<double> > b_, double c_, bool eq_, Pointer<char> name_=NULL)
			: b(1, b_), c(c_), eq(eq_), name(name_)
			{ }

			LinConstraint(const LinConstraint& lincon)
			: b(lincon.b), c(lincon.c), eq(lincon.eq), name(lincon.name)
			{ }
	};

	private:
		/** Parameters.
		*/
		Pointer<Param> param;

		/** The maximum number of cuts, we will add.
		*/
		int max_cutsnr;

		/** Given to CutPool's.
		*/
		int inactivetime_limit_global;
		int inactivetime_limit_local;

		vector<Pointer<CutPool> > cutpool;
		
		/** Cutpool of coupling constraints.
		*/
		Pointer<CutPool> cutpoolcoupling;

		/** Bounds on variables.
		*/
		dvector lower, upper;
		
		/** Indices of integer variables.
		*/
		vector<int> i_discr;

		Pointer<LinearRelaxSolver> solver;

		/** The address of the node, the LinearRelaxSolver uses currently.
				Do not dereference it, it could be an invalid pointer.
		*/
		MinlpNode* solver_node;

		/** The LinConstraint in couple_con or block_con, which keeps the level cut.
		    If no level cut is used, levelcut_pos points to the end of couple_con.
		*/
		list<LinConstraint>::iterator levelcut_pos;

		/** Reduces the box of a variable by minimizing/maximizing the variable using the linear relaxation.
		    @param oldlow Old lower bound of the variable.
		    @param oldup Old upper bound of the variable.
				@param k Block number of the variable.
				@param i Block index of the variable.
				@param discrete Indicates, whether the variable is a discrete one.
				@param unknown_only Indicates, whether we should only try to compute unknown bounds (instead of reducing known ones).
				@return The new box of the variable.
		*/
		bool box_reduce(pair<double,double>& newbox, const pair<double,double>& oldbox, int k, int i, bool discrete=false, bool unknown_only=false);

	public:
		/** The objective function of the linear relaxation.
		*/
		Pointer<SepQcFunc> obj;
		/** The constraints of the core linear relaxation, which are defined for one block only.
		    Here, the member b of a LinConstraint is a vector of length 1.
		*/
		vector<list<LinConstraint> > block_con;
		/** The constraints of the core linear relaxation, which are defined for more than one block.
		*/
		list<LinConstraint> couple_con;
		/** The size of the core of the linear relaxation.
				couple_con.size()+\sum_k block_con[k].size().
		*/
		int core_size() const;

		set<pair<int, int> > boxreduce_fixed_binaries;

		/** Constructor.
		*/
		LinearRelax(Pointer<Param> param_);

		/** Constructs a copy of one block of a LinearRelax which includes already all cuts for a specific node.
		    Copies the blockproblem from linrelax, adds the global cuts and the one for the given node (if not NULL) to this problem.
				So, this the core of the linear relax is (R_k[U]).
				@param alt_obj An alternative objective function.
		*/
		LinearRelax(const LinearRelax& linrelax, int k, Pointer<MinlpNode> node=NULL, Pointer<SepQcFunc> alt_obj=NULL);

		/** Initialize a linear relaxation from a convex problem.
		*/
		void init(Pointer<MinlpProblem> convex_prob, const vector<int>& i_discr, Pointer<dvector> feas_point=NULL, bool is_feas=true);

		/** Adds a level cut.
		*/
		void add_level_cut(double level);
		/** Updates the level cut.
		*/
		void update_level_cut(double newlevel);

		/** Adds a cut to the cutpool.
		    If the cut is global, or the node matches the one, the LinearRelaxSolver uses, the solver is notified.
		*/
		void add_cut(const Pointer<SimpleCut>& cut, int k, const Pointer<MinlpNode>& node=NULL);
		
		/** Adds a cut to the cutpool of coupling constraints.
		*/
		void add_cut(const Pointer<SimpleCut>& cut, const Pointer<MinlpNode>& node=NULL) {
			add_cut(cut, -1, node);
		}

		void add_cut(const Pointer<LinearizationCut>& cut, int k, const Pointer<MinlpNode>& node=NULL);

		/** Adds a IntervalGradientCut to the cutpool.
		    If the cut is global, or the node matches the one, the LinearRelaxSolver uses, the solver is notified.
		*/
		void add_cut(Pointer<IntervalGradientCut> intgradcut, int k, Pointer<MinlpNode> node=NULL);

		/** Integrated the cuts from another linear relaxation into this one.
				Adds the global cuts from the first cutpool in linrelax, which are not tagged, to the k'th cutpool here.
		*/
		void integrate(LinearRelax& linrelax, int k, Pointer<MinlpNode> node=NULL);

		/** Collects the cuts for a node.
		*/
		void get_cuts(list<CutPool::CutInfo>& cutinfos, Pointer<MinlpNode> node=NULL);

		void update_cuts(Pointer<MinlpNode> node, int k, IntervalGradientCutGenerator& generator, LinearizedConCutGenerator& linconcutgen);

		/** Computes an upper bound of the linear relaxation by maximizing the objective function.
		*/
		double get_upper_bound();

		bool cutlimit_reached() const { return nr_all_cuts()>=max_cutsnr; }

		/** The number of local cuts for a node.
		*/
		int nr_local_cuts(Pointer<MinlpNode> node) const;
		/** The number of global cuts.
		*/
		int nr_global_cuts() const;
		/** The number of all cuts, global and local.
		*/
		int nr_all_cuts() const;

		/** Removes a node from the cuts in the cutpool.
		*/
		void remove_node(Pointer<MinlpNode> node);

		void duplicate_nodeinfo(Pointer<MinlpNode> oldnode, Pointer<MinlpNode> newnode);

		/** Clears the solver such that it does not represent the linear relaxation of some node.
		*/
		void clear_solver();
		/** Checks the linear relaxation (R) or (R[U]) for feasibility.
			  @return True, if the linear relaxation is feasible, false else.
		*/
		bool feasible(Pointer<MinlpNode> node, double tol=1E-4);

		/** Checks, if a specific point is feasible for the linear relaxation.
				@return True, if the point is feasible according to the given tolerance. False, else.
		*/
		bool point_feasible(Pointer<MinlpNode> node, const dvector& x, double tol=1E-4) const;

		/** Solves the linear relaxation (R) or (R[U]).
		    @param sol_point To store the solution point.
				@param value To store the objective value.
				@param node To distinguish between (R) and (R[U]).
				@param dual_point A pointer to a dvector to store the dual variables in, or NULL, if not interested.
				@return A return code similiar to SNOPT. 0 for solved, 1 for infeasible, 3 else.
		*/
		int solve(dvector& sol_point, double& value, Pointer<MinlpNode> node, dvector* dual_point=NULL, double tol=1E-4);

		/** Reduces the box of a block of variables by minimizing/maximizing the variables using the linear relaxation.
		    @param node The node to use, can be NULL.
				@param k The block number.
				@param discr Indicators for discrete variables. Needs to be the one for the whole problem.
				@param unknown_only Indicates, whether we should try to determine unknown bounds only.
				@param changed_var A set to store the indices of variables in, which bounds were reduced by at least the percentage given in min_impr. Can be used as input for IntervalReduction.
				@param min_impr The relative minimum improvement a bound must achieve to be added to changed_var.
		*/
		double box_reduce(Pointer<MinlpNode> node, int k, const vector<bool>& discr, bool unknown_only=false, set<pair<int, IntervalReduction::which_bound_type> >* changed_var=NULL, double min_impr=0.01);
		/** Reduces the box of all variables by minimizing/maximizing the variables using the linear relaxation.
		    @param node The node to use, can be NULL.
				@param discr Indicators for discrete variables. Needs to be the one for the whole problem.
				@param unknown_only Indicates, whether we should try to determine unknown bounds only.
				@param changed_var A set to store the indices of variables in, which bounds were reduced by at least the percentage given in min_impr. Can be used as input for IntervalReduction.
				@param min_impr The relative minimum improvement a bound must achieve to be added to changed_var.
		*/
		double box_reduce(Pointer<MinlpNode> node, const vector<bool>& discr, bool unknown_only=false, set<pair<int, IntervalReduction::which_bound_type> >* changed_var=NULL, double min_impr=0.01);
		double box_reduce(dvector& newlow, dvector& newup, const vector<bool>& discr, bool unknown_only=false, set<pair<int, IntervalReduction::which_bound_type> >* changed_var=NULL, double min_impr=0.01);

		void set_box(const dvector& lower_, const dvector& upper_) {
			lower=lower_;
			upper=upper_;
		}
		
		void generate_cuts(Pointer<MinlpNode> node);
};

/** An abstract solver for the linear relaxation.
    Defines methods for adding constraints, resetting and constructing the linear relaxation depending on a MinlpNode.
		In addition to a solve()-method, it contains also a method for feasiblity-check.
*/
class LinearRelaxSolver : public Solver {
	protected:
		LinearRelax& linrelax;

	public:
		list<CutPool::CutInfo> cutinfos;
		
		/** Contains the dual variables of the core linear relaxation after a successful solve() or feasible() run.
		*/
		dvector duals;

		LinearRelaxSolver(LinearRelax& linrelax_);

		virtual ~LinearRelaxSolver();

		/** Initializes the model for the core linear relaxation.
		*/
		virtual void init()=0;

		/** Resets the solver, so we use only the core linear relaxation.
		*/
		virtual void reset()=0;

		/** Builds the linear relaxation for a specific node or uses global cuts only.
		*/
		virtual void construct(Pointer<MinlpNode> node=NULL)=0;

		/** Updates the linear relaxation for a specific node.
		*/
//		virtual void update(Pointer<MinlpNode> node=NULL)=0;

		virtual void add_cut(const CutPool::CutInfo& cutinfo)=0;

		virtual void set_objective(Pointer<SepQcFunc> obj)=0;

		virtual void update_levelcut(double newlevel)=0;

		/** Checks the linear relaxation only for feasiblity instead of solving it.
		*/
		virtual bool feasible()=0;

		virtual int solve(Pointer<MinlpNode> node)=0;

		/** According to the information, set in the remove-member of the cutinfos, remove cuts.
		*/
		virtual void remove_cuts()=0;

		int solve();

#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif

		virtual int generate_cuts(list<Pointer<SimpleCut> >& cuts)=0;
};

/** A general solver for the LinearRelaxation, which constructs a MinlpProblem and uses LocOpt::get_lp_solver to solve it.
*/
class LinearRelaxSolverGeneral : public LinearRelaxSolver {
	private:
		/** (R) or (R[U]).
		    In it's raw state (as constructed or set after reset), this problem just consists of the linear constraints of the convex relaxation.
				The method construct adds global and local cuts.
				Initialized in init().
				@see reset
				@see init
				@see construct
		*/
		Pointer<MinlpProblem> prob;

		/** Start of cuts in prob.
		    Set in init().
		*/
		int cutstart;
		/** The position of the level cut in prob.
		*/
		int levelcut_pos;

	public:
		LinearRelaxSolverGeneral(LinearRelax& linrelax_)
		: LinearRelaxSolver(linrelax_), cutstart(0), levelcut_pos(-1)
		{ tol=1E-4;
			init();
		}

		void init();

		void reset();

		void construct(Pointer<MinlpNode> node=NULL);

		void add_cut(const CutPool::CutInfo& cutinfo);

		void set_objective(Pointer<SepQcFunc> obj);

		void update_levelcut(double newlevel);

		bool feasible();

		int solve(Pointer<MinlpNode> node=NULL);

#if (!defined(__GNUC__)) // icc
		using LinearRelaxSolver::solve;
#endif
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif
		
		void remove_cuts();

		int generate_cuts(list<Pointer<SimpleCut> >& cuts) { return 0; }
};

class LinearRelaxSolverMIP : public LinearRelaxSolver {
	private:
		/** The core linear relaxation as MIP.
		*/
		MipProblem mip;

		Pointer<MIPSolver> solver;

		int levelcut_pos;
		
		/** Starting index of interval gradient cut variables.
		*/
		int intervalgradcutvars_start;

		void add_intervalgradientcut(CutPool::CutInfo& cutinfo);

	public:
		LinearRelaxSolverMIP(LinearRelax& linrelax_);
		void init();

		void reset();

		void construct(Pointer<MinlpNode> node=NULL);

		void add_cut(const CutPool::CutInfo& cutinfo);

		void set_objective(Pointer<SepQcFunc> obj);

		void update_levelcut(double newlevel);

		bool feasible();

		int solve(Pointer<MinlpNode> node=NULL);

#if (!defined(__GNUC__)) // icc
		using LinearRelaxSolver::solve;
#endif
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Solver::solve;
#endif
		
		void remove_cuts();
		
		int generate_cuts(list<Pointer<SimpleCut> >& cuts);
};

#endif
