// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOBOXMINIMIZATIONPROBLEM_HPP_
#define LAGOBOXMINIMIZATIONPROBLEM_HPP_

#include "LaGObase.hpp"
#include "LaGOSparsity.hpp"
#include "IpTNLP.hpp"

using Ipopt::Index;
using Ipopt::Number;

namespace LaGO {

/** Defines a problem which is the minimization of a function over a box.
 */
class BoxMinimizationProblem : public Ipopt::TNLP {
private:
	const Function& func;
	const DenseVector& lower;
	const DenseVector& upper;
	
	SmartPtr<SparsityGraph> sparsitygraph;
	
	DenseVector solution;
	double funcvalue;
	
public:
	const DenseVector* startpoint;

	BoxMinimizationProblem(const Function& func_, const DenseVector& lower_, const DenseVector& upper_, const SmartPtr<SparsityGraph>& sparsitygraph_=NULL)
	: func(func_), lower(lower_), upper(upper_), sparsitygraph(sparsitygraph_), solution(lower.getNumElements()), funcvalue(getInfinity()), startpoint(NULL)
	{ }
	
	bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style);
	
	bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);
	
	bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);
	
	bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

	bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values);

	bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values);

	void finalize_solution(Ipopt::SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq);
	
	const DenseVector& getSolution() { return solution; }
	double getOptimalValue() { return funcvalue; } 
}; // class BoxMinimizationProblem
	
} // namespace LaGO

#endif /*LAGOBOXMINIMIZATIONPROBLEM_HPP_*/
