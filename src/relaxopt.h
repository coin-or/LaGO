// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef RELAXOPT_H
#define RELAXOPT_H

#include "standard.h"
#include "opt.h"
#include "func.h"
#include "minlpopt.h"
#include "MINLPData.h"

//-------------------------------------------------------------------------

class PartitionHeu;
class MinlpNode;
class IntervalGradientCutInfo;

/** A base class for defining a relaxation-based solver.
*/
class RelaxationSolver : public Solver {
	friend class PartitionHeu;
	friend class MinlpNode;
	private:
		Pointer<LocOpt> LocOpt_NLP;

	protected:
		/** Parameters.
		*/
		Pointer<Param> param;

		/** The  original problem.
		*/
		Pointer<MinlpProblem> orig_prob; // (P)

		/** The block seperable formulation of the original problem.
		*/
		Pointer<MinlpProblem> split_prob; // (S)

		/** The convexified (quadratic) decomposed problem.
		*/
		Pointer<MinlpProblem> convex_prob; // (C)

		/** The linear relaxation.
		*/
		Pointer<LinearRelax> linear_relax; // (R)

		/** The reformulated problems.
		*/
		Pointer<Reformulation> reform;
		
		Pointer<MINLPData> minlpdata;

		Pointer<LevelCutHandler> levelcuts;

		bool is_gams_prob;

		int last_impr_iter;
		int point_conversion_warnings;

		void make_sub_prob(Pointer<MinlpProblem>& out, Pointer<MinlpProblem> in, Pointer<MinlpNode> node);

		pair<int, bool> locopt_NLP(dvector& start);
		pair<int, bool> locopt_NLP(dvector& start, dvector& sol_point);

		virtual bool add_sol_candidate(const dvector& x);

		void clean_sub_problems() { }

		double sol_cand_closeval_tol;
		Pointer<dvector> sol_cand_diam;

	public:
		set<SolCandidate> sol_cand;

		double lower_bound;

		RelaxationSolver(Pointer<MinlpProblem> orig_prob_, bool is_gams_prob=false,
			double closeval_tol_=0., Pointer<dvector> diam_=NULL, Pointer<Param> param_=NULL,
			Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		void set_split_prob(Pointer<MinlpProblem> split_prob_) {
			split_prob=split_prob_;
		}

		virtual void set_reform(Pointer<Reformulation> reform_) {
			reform=reform_;
		}

		void set_convex_prob(Pointer<MinlpProblem> convex_prob_) {
			convex_prob=convex_prob_;
		}

		void set_linear_relax(Pointer<LinearRelax> linear_relax_) {
			linear_relax=linear_relax_;
		}

		void set_upper_bound(double upper_bound) {
			opt_val_=upper_bound;
		}

		void set_levelcut_handler(Pointer<LevelCutHandler> levelcuts_) {
			levelcuts=levelcuts_;
		}
		
		void set_MINLPData(Pointer<MINLPData> minlpdata_) {
			minlpdata=minlpdata_;
		}

};

#endif
