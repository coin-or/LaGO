// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske


#ifndef BOXFIND_H
#define BOXFIND_H

#include "standard.h"
#include "param.h"
#include "problem.h"
#include "opt.h"
#include "relax.h"
#include "graph.h"

/** Finds bounds for variables.
    @class BoundsFinder
    @param BoundsFinder method
    %options guess, expensive, expensive2
    %default expensive
		%level 1
    The method, how we compute the bounds.
		@param Box Reduction skip binaries
		%options 0, 1
		%default 0
		%level 1
		Indicates, whether we should skip to perform box reduction for binary variables.
*/
class BoundsFinder {
	private:
		/** Parameters.
		*/
		Pointer<Param> param;

		/** The method to use to find and update bounds.
		*/
		int method;

		/** Computes a new bound for one variables.
		    Calls SNOPT for the problem min { (low ? 1 : -1) * x_(block,index) s.t. con's from conv_prob }.
		    @param conv_prob A convex problem, we can give to SNOPT after changing it's objective.
		    @param ret An integer to store the return value from SNOPT in.
		    @param index The index (in a block) of the variable to compute the bound for.
		    @param block The block number of the variables to compute the bound for.
		    @param low Indicates, wheter we should compute the lower (true) or upper bound.
		    @return The new bound.
		*/
		double compute_bound(MinlpProblem& conv_prob, int& ret, int index, int block, bool low);

	public:
		bool low, up, known;

		/** (Standard-)Constructor.
		    @param param_ Parameters, default is NULL.
		*/
		BoundsFinder(Pointer<Param> param_=NULL)
		: param(param_), method(3), low(true), up(true), known(true)
		{	if (!LocOpt::nlp_solver_available()) method=0;
			else if ((!param) || (strcmp(param->get("BoundsFinder method", "expensive"), "expensive")==0))	method=1;
			else if (!strcmp(param->get("BoundsFinder method"), "guess"))	method=0;
			else if (!strcmp(param->get("BoundsFinder method"), "expensive2")) method=2;
		};
		
		/** Tries to find or improve variables bounds.
		    Depending on the parameter "BoundsFinder method" calls the simple or expensive method.
		    @param prob A problem, for which variable bounds should be computed.
				@param discr Indicates, which variables were discrete in the original problem.
		    @return The number of bounds, where we had problems to compute a new bound and the number of guessed bounds.
		*/
		pair<int,int> compute_bounds(MinlpProblem& prob, vector<bool>& discr);

		/** Tries to find missing bounds by guessing.
		    Computes the minimum/maximum of all existing lower/upper bounds.
		    A guess for a missing lower bound is than 10*MIN(min_lower, -1).
		    A guess for a missing upper bound is than 10*MAX(min_lower, -1).
		    @param prob A MinlpProblem to guess missing bounds for.
		    @return The number guessed bounds in the second part.
		*/
		pair<int,int> compute_bounds_guess(MinlpProblem& prob);

		/** Tries to find or improve variable bounds.
			  For each variable x, tries to solve the problems max/min { x s.t. convex-constraints from conv_prob. } with SNOPT.
			  If SNOPT returns 0, 3 or 9, it takes the solution from snopt as a new bound, if it doesn't exceed a limit.
			  If SNOPT returns with other code, it guesses a new bound.

			  In normal, a new bound of variable j is not used when we compute the bounds of variables after j.
			  But if the variable was unbounded or SNOPT returns with 0, we use the new bounds also to compute other bounds.

			  So the bounds in conv_prob are changed, when we are sure.
		    @param conv_prob A convex problem, which objective is changed during this method.
		    @param lower A dvector to store the computed lower bounds in, even the one, where we are not sure.
		    @param upper A dvector to store the computed upper bounds in, even the one, where we are not sure.
				@param discr Indicates, which variables were discrete in the original problem.
		    @return The number of bounds, where we had problems to compute a new bound (LocOpt return != 0) and the number of guessed bounds.
		*/
		pair<int,int> compute_bounds_expensive(MinlpProblem& conv_prob, dvector& lower, dvector& upper, vector<bool>& discr);

		/** Tries to find or improve variable bounds.
		    Runs compute_bounds_expensive as often, as the number of errors is decreasing.
		    @param conv_prob A convex problem, which variable bounds should be improved.
		    @param lower A dvector to store the computed lower bounds in, even the one, where we are not sure.
		    @param upper A dvector to store the computed upper bounds in, even the one, where we are not sure.
				@param discr Indicates, which variables were discrete in the original problem.
		    @return The (final) number of bounds, where we had problems to compute a new bound (LocOpt return != 0) and the number of guessed bounds.
		*/
		pair<int,int> compute_bounds_expensive2(MinlpProblem& conv_prob, dvector& lower, dvector& upper, vector<bool>& discr);

		/** Computes bounds using linear constraints.
		    Computes bounds of variables in constraints of the form bx+c <= 0.
				@param prob The problem, for which we compute the new bounds.
		*/
		void compute_linbounds(MinlpProblem& prob);
};

/** Computes dual bounds.
*/
class DualBounds {
	private:
		Pointer<Param> param;

		Pointer<MinlpProblem> S, C, R;
		vector<Pointer<MinlpProblem> > lag_prob;

		/** Indices of linear constraints for each block.
	      Constraints, where s_k==A_k==NULL, but b_k!=NULL.
	  */
		vector<list<int> > lin_con;
		/** Indicates for each constraint of (S) or (C), whether it is linear or not.
		*/
		vector<bool> lin_con2;

		int type;

		double dual_bound(UserVector<double>& a);

	public:
		DualBounds(Pointer<MinlpProblem> S_, Pointer<MinlpProblem> C_, Pointer<Param> param_, int type_=1);
	
		int update_box(int n);
		int update_box() { return update_box(S->dim()); }
	
		double obj_bound();
};


class IntervalReduction {
	public:
		typedef enum {LOWER, UPPER, WHATEVER} which_bound_type;
	private:
		class NodeData {
			friend ostream& operator<<(ostream& out, const NodeData& nd) {
				if (nd.block_nr>=0) out << nd.block_nr << ',' << nd.block_index;
				return out;
			}
			public:
				int block_nr, block_index;

				NodeData(int block_nr_=-1, int block_index_=-1)
				: block_nr(block_nr_), block_index(block_index_)
				{ }

				NodeData& operator+=(const NodeData& nd) { return *this; }
		};

		class EdgeData {
			friend ostream& operator<<(ostream& out, const EdgeData& ed);
			public:
				which_bound_type which_bound;
				int con_nr;
				double coeff; // the coefficient before the linear variable, this edge is directed to

				EdgeData(int con_nr_, which_bound_type which_bound_, double coeff_)
				: con_nr(con_nr_), which_bound(which_bound_), coeff(coeff_)
				{ }

				bool operator<(const EdgeData& ed) const { return con_nr<ed.con_nr; }

				const EdgeData& operator+=(const EdgeData& ed) const;
		};
		friend ostream& operator<<(ostream& out, const EdgeData& ed);

		Pointer<MinlpProblem> prob;

		typedef Graph<NodeData,EdgeData,true,true> DependencyGraph;
		DependencyGraph dependency_graph;

		void run(dvector& newlow, dvector& newup, const dvector& oldlow, const dvector& oldup, set<pair<const DependencyGraph::NodeType*, which_bound_type> >& nodeset);

	public:
		bool do_print;
		bool empty_boxes;
		double min_impr;
		set<pair<int, int> > fixed_binaries;

		/** Stores for each block, how much it was reduced.
		*/
		dvector reduction_by_block;
		double reduction;

		IntervalReduction(Pointer<MinlpProblem> prob_=NULL)
		: do_print(false), empty_boxes(false), min_impr(.01)
		{ if (prob_) set_problem(prob_);
		}

		void set_problem(Pointer<MinlpProblem> prob_);

		void compute(dvector& newlow, dvector& newup, const dvector& oldlow, const dvector& oldup);
		void compute(dvector& newlow, dvector& newup, const dvector& oldlow, const dvector& oldup, set<pair<int, which_bound_type> >& startset);

		void print_small_boxes(dvector& low, dvector& up);
};


/** Computes missing bounds or reduces existing bounds by minimizing/maximizing a variable over some linear constraints.
*/
class BoundsFinderLinear {
	private:
		Pointer<Param> param;

		Pointer<MinlpProblem> prob;

		Pointer<MIPSolver> solver;

		MIPSolver::SolutionStatus compute_bound(double& bound, int index, bool low);

	public:
		bool low, up, known;

		BoundsFinderLinear(Pointer<MinlpProblem> prob_, Pointer<Param> param_);

		int compute(dvector& newlow, dvector& newup, set<pair<int, IntervalReduction::which_bound_type> >* changed_var=NULL, double min_impr=0.01);

};

#endif
