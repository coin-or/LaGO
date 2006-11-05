// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef RMP_H
#define RMP_H

#include "standard.h"
#include "linrelax.h"
#include "node.h"
#include "opt.h"

/** Class to manage the RMP.
    @class RMPManager
		@param RMP delta factor
		%options $\geq 0$
		%default 1E+6
		The penalty parameter for the y in the RMP objective function.
*/
class RMPManager {
	private:
		Pointer<LinearRelax> linear_relax;
		Pointer<MinlpNode> node;
		Pointer<Param> param;

		double y_penalty_delta, y_penalty_scale;

		int coresize;
		int blocknr;
		int linrelaxdim;

		Pointer<MIPSolver> RMP;

		Pointer<SparseVector<double> > compute_col(const ExtremePoint& w, int k);
		void init();

	public:
		RMPManager(Pointer<LinearRelax> linear_relax_, Pointer<MinlpNode> node_, Pointer<Param> param_);

		int dim() const { return RMP->nr_col(); }

		int colindex(const ExtremePoint& w);

		void add_column(ExtremePoint& w, int k);

		void update_column(const ExtremePoint& w, int k);

		void update_s(const dvector& x, bool set_start=false);
		void restrict_y(const dvector& y);

		MIPSolver::SolutionStatus solve() { return RMP->solve(); };
		MIPSolver::SolutionStatus solve(const UserVector<double>& start) { return RMP->solve(start); }

		void get_yz(dvector& yz) { RMP->get_primal(yz); }
		void get_dual(UserVector<double>& mu) { RMP->get_dual(mu); };
		double get_optval() { return RMP->get_optval(); };

		void get_x(UserVector<double>& x, const dvector &yz);

		int prune();
};

#endif
