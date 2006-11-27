// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef COLUMN_H
#define COLUMN_H

#include "standard.h"
#include "problem.h"
#include "func.h"
#include "opt.h"
//#include "node.h"
#include "minlpopt.h"
#include "bcp.h"
#include "rmp.h"

class MinlpNode;


/** Class for generating and solving a restricted master problem
    @class ColumnGenerator
		@param ColumnGenerator max major iter
		%options integer $\geq 0$
		%default 5
		How often $\mu$ is maximally updated during one generate\_RMP run.
		@param ColumnGenerator max minor iter
		%options integer $\geq 0$
		%default 0
		For how many blocks the lagrangian subproblem is solved. Only iterations, where the gap could be significantly reduced, are counted.
		If 0, the number is computed as $\max(5, 1.39\frac{n}{p}-.05(\frac{n}{p})^2+0.15p)$
		@param ColumnGenerator max init RMP iter
		%options integer $\geq 0$
		%default 0
		%level 0
		Iteration limit to make RMP feasible. If 0, the number is computed as 2*(maxblocksize+1).
*/
class ColumnGenerator {
	private:
		/** A reference to MinlpNode, where the sub-problem is defined.
		*/
		Pointer<MinlpNode> node;

		MinlpBCP& bcp;

		Pointer<Param> param;

		/** The restricted master problem */
		Pointer<RMPManager> RMP;

		int max_major_iter;
		int max_minor_iter;
		int max_initRMP_iter;

		double rc_tol;

		void add_cut(Pointer<SepQcFunc> cut, int k);

		int get_search_dir2(dvector& a, double& c, int k);

	public:
		ColumnGenerator(Pointer<MinlpNode> node_, MinlpBCP& bcp_);

		/** Solves the relaxed RMP and tries to find new RMP points, if y is not 0.
		    @return 0, if we could make the RMP feasible. 1, if (R[U]) is feasible, but the (RMP) not. 2, if the (R[U]) is feasible, but its solution exceeds the upper bound. 3, if (R[U]) is infeasible.
		*/
		int init_RMP();
		int solve_RMP();

		void clear_RMP() { RMP=NULL; }

		void solve_lag_problem(MinlpBCP::LagSolveStatus& status, int k, Pointer<SepQcFunc> temp_cut=NULL);

		/** Update the RMP points after subdivision, and adds the columns to the RMP, if the RMP was built before. */
		void update_ExtremePoints();

		/** Computes Wz=sum_l z_{k,l} * w_{k,l} for a fixed k.
		    @param x The vector to store the result in. Must have (R).block[k].size().
		    @param z The vector which is multiplicated by W. Must have (RMP).block[k].size()==node.i_RMP_points[k].size().
		*/
		void mult_W(UserVector<double>& x, const dvector &z, int k);

		/** Generates the RMP by column generation */
		int generate_RMP();

		/** A key for sorting Lagrangian subproblems
		    @param k The number of the subproblem
		*/
		double subprobl_key(double rc, int k);

		void print(ostream& out) const;
};

#endif
