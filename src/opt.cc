// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

// opt.cc

#include "opt.h"

#include "snopt.h"
#include "ipopt2.h"
#include "osi.h"

extern bool snoptlicenceok;

// ----------------------------- LocOpt ----------------------------------------

bool LocOpt::nlp_solver_available() {
#ifdef SNOPT_AVAILABLE
	if (snoptlicenceok) return true;
#endif
#ifdef IPOPT_AVAILABLE
	return true;
#else
	return false;	
#endif
}

Pointer<LocOpt> LocOpt::get_solver(const Pointer<MinlpProblem> prob, Pointer<Param> param, char* param_prefix, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_) {
#ifdef SNOPT_AVAILABLE
	if (snoptlicenceok)
		return new SnOpt(prob, param, param_prefix, out_solver_p_, out_solver_log_p_);
#endif
#ifdef IPOPT_AVAILABLE
	return new IpOpt(prob, param, out_solver_p_, out_solver_log_p_);
#endif
	out_err << "Sorry, no Solver for local optimization available. Aborting." << endl;
	exit(-1);
	return NULL;
}

Pointer<LocOpt> LocOpt::get_lp_solver(const Pointer<MinlpProblem> prob, Pointer<Param> param, char* param_prefix, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_) {
#if defined(CPLEX_AVAILABLE) && defined(CONCERT_AVAILABLE)
	return new Cplex(prob, param, out_solver_p_, out_solver_log_p_);
#endif
#if defined(COIN_AVAILABLE)
	return new LPSolver(prob);
#endif
#if defined(SNOPT_AVAILABLE)
	return new SnOpt(prob, param, param_prefix, out_solver_p_, out_solver_log_p_);
#endif
  out_err << "Sorry, no LP Solver available. Aborting." << endl;
	exit(-1);
	return NULL;
}

Pointer<LocOpt> LocOpt::get_solver_origprob(const Pointer<MinlpProblem> prob, Pointer<Param> param, char* param_prefix, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_) {
#ifdef GAMS_AVAILABLE
	if (gamsptr)
		return new gamsLocOpt(prob, param, out_solver_p_, out_solver_log_p_);
#endif
	return get_solver(prob, param, param_prefix, out_solver_p_, out_solver_log_p_);
}

// ----------------------------- MIPSolver ----------------------------------------

Pointer<MIPSolver> MIPSolver::get_solver(const MipProblem& mip, Pointer<Param> param) {
#ifdef COIN_AVAILABLE
	return new OSISolver(mip);
#endif
  out_err << "Sorry, no MIPSolver available. Aborting." << endl;
	exit(-1);
	return NULL;
}


// --------------------------- DualSolver -------------------------

void DualSolver::do_log() {
  if (log_frequency && !(improve_iter % log_frequency)) {
    if (dual_vals) dual_vals->push_back(obj.get_dual_val());
    if (dual_points) dual_points->push_back(obj.get_dual_point());
	  if (orig_points) {
	    vector<dvector> op;
	    for (int i=0; i<obj.nr_of_orig_points(); i++)
	      op.push_back(obj.get_orig_point(i));
	    orig_points->push_back(op);
    }

	 	out_solver_log << "adding to sample set in iter " << improve_iter << ": " << obj.get_dual_val() << endl;
	}
}

int DualSolver::check(double val) {
	out_solver_log << "Iter: " << iter() << " Value: " << val << endl;
  if (val>last_val) {  // serious step ?
	  improve_iter++;
		last_val=val;  // update last serious value
	}
	else return 0;  // check only non-null-steps

	do_log();

	if (threshold_cntrl && val>=threshold) return 100;

	if (conv_rate_cntrl) {
		if (improve_iter==0) { // first iteration
		  out_solver_log << "Iter (improve): " << iter() << " (first) Value " << val << endl;
			last_major_val=first_major_val=val;
			return 0;
		}
		if (improve_iter%minor_iter == 0) { // major iteration
			out_solver_log << "Iter (improve): " << iter() << " (" << improve_iter << ") Value: " << val << "  ";

			if (improve_iter==minor_iter) { // first major iteration
				max_rel_improvement=(val-first_major_val)/(fabs(first_major_val)+1);
				last_major_val=val;
				out_solver_log << "First improvement: " << max_rel_improvement << endl;
				return 0;
			}
		  last_rel_improvement=(val-last_major_val)/(fabs(first_major_val)+1);  // compute relative improvment
			last_major_val=val;
		  if (last_rel_improvement > max_rel_improvement) {  // update max improvement
		    max_rel_improvement=last_rel_improvement;
		    out_solver_log << "Updated maximum improvement: " << max_rel_improvement << endl;
		    return 0;
		  }
			out_solver_log << last_rel_improvement << " < " << max_rel_improvement*stopping_rho << " ?" << endl;
			if (last_rel_improvement < max_rel_improvement*stopping_rho) return 10;
		}
	}

  return 0;
}

// ---------------------------------------- LPSolver ----------------------------------

void LPSolver::initmip() {
	mipsolver=NULL;
#ifdef COIN_AVAILABLE
	mipsolver=new OSISolver(MipProblem(*prob));
#else
	out_err << "Sorry, no MIP Solver available." << endl;
#endif
//	mipsolver->set_tol(tol);
	mipsolver->set_maxiter(iter_max);
}

int LPSolver::solved(MIPSolver::SolutionStatus status) {
	status=mipsolver->solve(); // sometimes, a resolve might be necessary; but what is "sometimes" :-(
	if (status==MIPSolver::SOLVED || status==MIPSolver::FEASIBLE) {
		mipsolver->get_primal(sol_point);
		opt_val_=mipsolver->get_optval();
	}
	iter_=mipsolver->get_iter();
	switch (status) { // translate to SNOPT like return codes
		case MIPSolver::SOLVED: return 0;
		case MIPSolver::FEASIBLE: return 4;
		case MIPSolver::UNBOUNDED: return 2;
		case MIPSolver::INFEASIBLE: return 1;
		case MIPSolver::ITERATIONLIMITEXCEEDED: return 3;
		case MIPSolver::ABORTED: return 5;
		case MIPSolver::UNKNOWN: return 5;
	}
	out_err << "LPSolver: SolutionStatus " << status << " not known!!"; return 5;
}

LPSolver::LPSolver(const Pointer<MinlpProblem> prob_/*, Pointer<Param> param_, char* param_prefix_*/)
: LocOpt(prob_->dim()), prob(prob_)//, param(param_), param_prefix(param_prefix_ ? strdup(param_prefix_) : NULL)
{ initmip();
}

dvector LPSolver::get_lag_multipliers() {
	dvector mu(prob->con.size());
	mipsolver->get_dual(mu);
	return mu;
}
