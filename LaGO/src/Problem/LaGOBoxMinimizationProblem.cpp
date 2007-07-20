// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOBOXMINIMIZATIONPROBLEM_CPP_
#define LAGOBOXMINIMIZATIONPROBLEM_CPP_

#include "LaGOBoxMinimizationProblem.hpp"
#include "LaGOSymSparseMatrix.hpp"

using namespace Ipopt;

namespace LaGO {

bool BoxMinimizationProblem::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {
	n=lower.getNumElements();
	m=0;
	nnz_jac_g=0;
	index_style=Ipopt::TNLP::C_STYLE;
	if (IsValid(sparsitygraph)) {
		nnz_h_lag=0;
		// unfortunately, in the sparsity graph, coupling between different variables are stored as two directed arcs
		for (SparsityGraph::arc_iterator it(sparsitygraph->arc_begin()); it!=sparsitygraph->arc_end(); ++it) {
			int var1=(**(*it).head()).varindex;
			int var2=(**(*it).tail()).varindex;
			if (var1<=var2) ++nnz_h_lag;
		}
	}
	else
		nnz_h_lag=n*(n+1)/2;
	return true;
}

bool BoxMinimizationProblem::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
	CoinCopyN(lower.getElements(), n, x_l);
	CoinCopyN(upper.getElements(), n, x_u);
	return true;
}

bool BoxMinimizationProblem::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
	if (init_x) {
		if (startpoint) {
			assert(startpoint->getNumElements()==n);
			CoinCopyN(startpoint->getElements(), n, x);
		} else {
			for (int i=0; i<n; ++i)
				x[i]=(lower[i]+upper[i])/2;
		}
	}
	init_z=false;
	init_lambda=false;
	return true;
}

bool BoxMinimizationProblem::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
	DenseVector x_(n, x);
	try {
		obj_value=func.eval(x_);
	}	catch (FunctionEvaluationError error) {
		return false;
	}
	return true;
}

bool BoxMinimizationProblem::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
	DenseVector x_(n, x);
	DenseVector grad_(n);
	try {
		func.gradient(grad_, x_);
		CoinCopyN(grad_.getElements(), n, grad_f);
	}	catch (FunctionEvaluationError error) {
		return false;
	}
	return true;
}

bool BoxMinimizationProblem::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
	return true;
}

bool BoxMinimizationProblem::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values) {
	return true;
}

bool BoxMinimizationProblem::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values) {
	if (values==NULL) { // give hessian sparsity info
		if (IsValid(sparsitygraph)) {
			for (SparsityGraph::arc_iterator it(sparsitygraph->arc_begin()); it!=sparsitygraph->arc_end(); ++it) {
				int var1=(**(*it).head()).varindex;
				int var2=(**(*it).tail()).varindex;
				if (var1<=var2) { *jCol++=var1; *iRow++=var2; }
				else continue;
			}
		} else
			for (int row=0; row<n; ++row)
				for (int col=0; col<=row; ++col, ++iRow, ++jCol) {
					*iRow=row;
					*jCol=col;
				}
	} else {
		if (obj_factor==0.) {
			CoinZeroN(values, nele_hess);
			return true;
		}
		SymSparseMatrixCreator hessian(n);
		try {
			DenseVector x_(n, x);
			func.fullHessian(hessian, x_);
		} catch (FunctionEvaluationError error) {
			return false;
		}
		hessian.cleanUpperDiagonal();

		if (IsValid(sparsitygraph)) {
			for (SparsityGraph::arc_iterator it(sparsitygraph->arc_begin()); it!=sparsitygraph->arc_end(); ++it) {
				int var1=(**(*it).head()).varindex;
				int var2=(**(*it).tail()).varindex;
				if (var1>var2) continue;
				pair<int,int> col_row(var2, var1); // var2<=var1
				map<pair<int,int>, double>::iterator it_h(hessian.find(col_row));
				if (it_h!=hessian.end()) *values=obj_factor*it_h->second;
				else *values=0.;
				++values;
			}
		} else {
			CoinZeroN(values, nele_hess);
			for (map<pair<int,int>, double>::iterator it_h(hessian.begin()); it_h!=hessian.end(); ++it_h) {
				int row=it_h->first.first;
				int col=it_h->first.second;
				int index=row*(row+1)/2+col;
				assert(index<nele_hess);
				values[index]=obj_factor*it_h->second;
			}
		}
	}
	return true;
}

void BoxMinimizationProblem::finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq) {
	if (x)
		CoinCopyN(x, n, solution.getElements());
	funcvalue=obj_value;
}

} // namespace LaGO

#endif /*LAGOBOXMINIMIZATIONPROBLEM_CPP_*/
