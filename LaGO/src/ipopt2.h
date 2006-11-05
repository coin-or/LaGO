// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#define IPOPT_AVAILABLE
#ifdef IPOPT_AVAILABLE
#ifndef IPOPT_H
#define IPOPT_H

#include "standard.h"
#include "opt.h"
#include "problem.h"
#include "param.h"

#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"
#include "IpSolveStatistics.hpp"

using namespace Ipopt;

class IpOpt;

class IpOptProblem : public TNLP {
	friend class IpOpt;
	private:
		const Pointer<MinlpProblem> prob;

		IpOpt& ipopt;

		vector<Pointer<SparsityInfo> > sparsity;
		Index nnz_jac_g;
		Index nnz_h_lag;

		/** Sets or adds the hessian of one block of a function.
		    @param values Where to put the block of the hessian to.
		    @param factor Scalar factor for the hessian.
		    @param xx The vector to evaluate the hessian in.
		    @param y Scratch vector.
		    @param z Scratch vector.
		*/
		template <bool add>
		void set_hessian(const SepQcFunc& func, int blocknr, Number* values, double factor, dvector& xx, dvector& y, dvector& z);

		/** Sets or adds the hessian of one block of a quadratic function multiplied by a vector.
		    @param values Where to put the block of the hessian to.
		    @param factor Scalar factor for the hessian.
		*/
		template <bool add>
		void set_hessianquad(const SepQcFunc& func, int blocknr, Number* values, double factor);

	public:
		IpOptProblem(const Pointer<MinlpProblem> prob_, IpOpt& ipopt);

		/** Method to return some info about the nlp */
		virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
						Index& nnz_h_lag, IndexStyleEnum &index_style);

		/** Method to return the bounds for my problem */
		virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
							Index m, Number* g_l, Number* g_u);

		/** Method to return the starting point for the algorithm */
		virtual bool get_starting_point(Index n, bool init_x, Number* x,
																			bool init_z, Number* z_L, Number* z_U,
																			Index m, bool init_lambda, Number* lambda);

		/** Method to return the objective value */
		virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

		/** Method to return the gradient of the objective */
		virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

		/** Method to return the constraint residuals */
		virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

		/** Method to return:
		*   1) The structure of the jacobian (if "values" is NULL)
		*   2) The values of the jacobian (if "values" is not NULL)
		*/
		virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
					Index m, Index nele_jac, Index* iRow, Index *jCol,
					Number* values);

		/** Method to return:
		*   1) The structure of the hessian of the lagrangian (if "values" is NULL)
		*   2) The values of the hessian of the lagrangian (if "values" is not NULL)
		*/
		virtual bool eval_h(Index n, const Number* x, bool new_x,
						Number obj_factor, Index m, const Number* lambda,
						bool new_lambda, Index nele_hess, Index* iRow,
						Index* jCol, Number* values);

		/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
		virtual void finalize_solution(SolverReturn status,
					Index n, const Number* x, const Number* z_L, const Number* z_U,
					Index m, const Number* g, const Number* lambda,
					Number obj_value);

};

class IpOpt : public LocOpt {
	friend class IpOptProblem;
	private:
		Pointer<Param> param;

		SmartPtr<IpoptApplication> ipopt;
		SmartPtr<TNLP> ipoptproblem;

		dvector lambda;

	public:
		IpOpt(const Pointer<MinlpProblem> prob_, Pointer<Param> param_, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

		~IpOpt() {}

		void reinit();

		int solve();
		int solve(dvector& start);
		dvector get_lag_multipliers() { return lambda; }
	
};

#endif
#endif
