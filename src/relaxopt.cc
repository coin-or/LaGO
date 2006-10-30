// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "relaxopt.h"
#include "bcp.h" // to get IntervalGradientCutInfo

// ------------------------------ RelaxationSolver -----------------------------

RelaxationSolver::RelaxationSolver(Pointer<MinlpProblem> orig_prob_, bool is_gams_prob_,
	double closeval_tol_, Pointer<dvector> diam_, Pointer<Param> param_,
	Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: Solver(orig_prob_->dim(), out_solver_p_, out_solver_log_p_), orig_prob(orig_prob_), is_gams_prob(is_gams_prob_),
  sol_cand_closeval_tol(closeval_tol_), sol_cand_diam(diam_), param(param_),
	last_impr_iter(-1), point_conversion_warnings(0), lower_bound(-INFINITY)
{	tol=1E-4;
	if (!sol_cand_closeval_tol) sol_cand_closeval_tol=param->get_d("heu close value tolerance", .001);
	if (!sol_cand_diam) {
		sol_cand_diam=new dvector(orig_prob->upper); *sol_cand_diam-=orig_prob->lower;
		*sol_cand_diam*=param->get_d("heu close points tolerance", .001);
	}
}

bool RelaxationSolver::add_sol_candidate(const dvector& x) {
	double val;

	int oldsize=sol_cand.size();
	if (x.dim()==orig_prob->dim()) { // short vector
		val=orig_prob->obj->eval(x);
		sol_cand.insert(SolCandidate(val, x, sol_cand_closeval_tol, sol_cand_diam));
	} else { // long vector
		assert(x.dim()==reform->ext_prob->dim());
		val=reform->ext_prob->obj->eval(x);
		sol_cand.insert(SolCandidate(val, reform->get_short_vector(x), sol_cand_closeval_tol, sol_cand_diam));
	}
	if (oldsize==sol_cand.size()) return true; // no new point

	if (val<opt_val_-rtol) { // better upper bound
		out_solver << "Added " << sol_cand.size() << ". Solution Candidate, value=" << val << endl;
		opt_val_=val;
		last_impr_iter=iter();

		if (levelcuts) levelcuts->update_level_cut(val);
	}

	return false;
}

pair<int, bool> RelaxationSolver::locopt_NLP(dvector& start) {
	pair<int, bool> ret;
	if (orig_prob->i_discr.size()==orig_prob->dim()) { // all binaries
		if (start.dim()>orig_prob->dim()) {
			assert(reform);
			dvector short_start(reform->get_short_vector(start));
			if (orig_prob->feasible(short_start, tol, NULL)) ret.first=1;
			else {
				ret.first=0;
				ret.second=add_sol_candidate(reform->get_long_vector(short_start));
			}
		} else {
			if (orig_prob->feasible(start, tol, NULL)) ret.first=1;
			else {
				ret.first=0;
				ret.second=add_sol_candidate(start);
			}
		}
		return ret;
	}

	if (!LocOpt_NLP) {
		if (is_gams_prob) LocOpt_NLP=LocOpt::get_solver_origprob(orig_prob, param, NULL, NULL, NULL);
		else LocOpt_NLP=LocOpt::get_solver(orig_prob, param, "LocOpt", NULL, NULL);
	}

	if (start.dim()>orig_prob->dim()) {
		assert(reform);
		dvector short_start(reform->get_short_vector(start));
//out_log << "starting from " << endl;
//for (int i=0; i<short_start.dim(); i++) out_log << orig_prob->var_names[i] << "\t " << short_start(i) << endl;
		ret.first=LocOpt_NLP->solve(short_start);
	} else {
		ret.first=LocOpt_NLP->solve(start);
	}

	if ((!ret.first) || (!orig_prob->feasible(LocOpt_NLP->sol_point, tol, NULL)))
		ret.second=add_sol_candidate(start.dim()>orig_prob->dim() ? reform->get_long_vector(LocOpt_NLP->sol_point) : LocOpt_NLP->sol_point);

	return ret;
}

pair<int, bool> RelaxationSolver::locopt_NLP(dvector& start, dvector& solpoint) {
	pair<int, bool> ret;
	if (orig_prob->i_discr.size()==orig_prob->dim()) { // all binaries
		if (start.dim()>orig_prob->dim()) {
			assert(reform);
			dvector short_start(reform->get_short_vector(start));
			if (orig_prob->feasible(short_start, tol, NULL)) ret.first=1;
			else {
				ret.first=0;
				solpoint=reform->get_long_vector(short_start);
			}
		} else {
			if (orig_prob->feasible(start, tol, NULL)) ret.first=1;
			else {
				ret.first=0;
				solpoint=start;
			}
		}
		if (!ret.first) ret.second=add_sol_candidate(solpoint);
		return ret;
	}

	if (!LocOpt_NLP) {
		if (is_gams_prob) LocOpt_NLP=LocOpt::get_solver_origprob(orig_prob, param, NULL, NULL, NULL);
		else LocOpt_NLP=LocOpt::get_solver(orig_prob, param, "LocOpt", NULL, NULL);
	}
	if (start.dim()>orig_prob->dim()) {
		assert(reform);
		dvector short_start(reform->get_short_vector(start));
		ret.first=LocOpt_NLP->solve(short_start);
		solpoint=reform->get_long_vector(LocOpt_NLP->sol_point);
		if ((!ret.first) || (!orig_prob->feasible(LocOpt_NLP->sol_point, tol, NULL))) {
			if (reform->ext_prob->feasible(solpoint, tol, NULL)) {
				out_err << "Warning: Point feasible for original problem, but not for extended formulation." << endl;
				point_conversion_warnings++;
//				reform->ext_prob->feasible(solpoint, tol, out_log_p);
			}
			ret.second=add_sol_candidate(solpoint);
		}
		return ret;
	} else {
		ret.first=LocOpt_NLP->solve(start);
		solpoint=LocOpt_NLP->sol_point;
		if ((!ret.first) || (!orig_prob->feasible(solpoint, tol, NULL)))
			ret.second=add_sol_candidate(solpoint);
	}
	return ret;
}

void RelaxationSolver::make_sub_prob(Pointer<MinlpProblem>& out, Pointer<MinlpProblem> in, Pointer<MinlpNode> node) {
	out=new MinlpProblem(*in);
	out->lower=node->lower;
	out->upper=node->upper;
	for (int i=0; i<out->dim(); ++i)
		if (out->primal_point(i)<out->lower(i)-rtol) out->primal_point[i]=out->lower(i);
		else if (out->primal_point(i)>out->upper(i)+rtol) out->primal_point[i]=out->upper(i);
	for (int k=0; k<node->part_con.size(); k++)
		for (list<Pointer<SepQcFunc> >::iterator it(node->part_con[k].begin()); it!=node->part_con[k].end(); it++)
			out->add_con(*it, false, "partition con");
}

