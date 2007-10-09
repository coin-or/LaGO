// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOQUADRATICORCONVEXAPPROXIMATION_HPP_
#define LAGOQUADRATICORCONVEXAPPROXIMATION_HPP_

#include "BonTMINLP.hpp" 

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"


namespace LaGO {
	
/** Representing the approximation of the MINLP where each function is either quadratic or nonquadratic but convex in a "scalar" format that should be easier to process then the "condensed" MinlpData object.
 */
class QuadraticOrConvexApproximation : public Bonmin::TMINLP {
	friend ostream& operator<<(ostream& out, const QuadraticOrConvexApproximation& quad) {
		quad.print(out);
		return out;
	} 
private:
	MINLPData& data;
	/** Whether we should build the convex quad. approx., or just the quad. approx.
	 */
	bool use_convex;
	
	/** Lower bounds on variables.
	 */
	DenseVector var_lb;
	/** Upper bounds on variables.
	 */
	DenseVector var_ub;
	/** Lower bounds on constraints.
	 */
	DenseVector con_lb;
	/** Upper bounds on constraints.
	 */
	DenseVector con_ub;
	/** Quadratic (and linear and constant) part of objective function.
	 */
	SmartPtr<QuadraticFunction> objQuad;
	/** Convex nonquadratic part of objective function.
	 */
	SmartPtr<Function> objNonQuad;
	/** Sparsity graph of nonquadratic part of constraints.
	 */
	SmartPtr<SparsityGraph> objNonQuadSparsityGraph;
	/** Quadratic (and linear and constant) part of constraints.
	 */
	vector<SmartPtr<QuadraticFunction> > conQuad;
	/** Convex nonquadratic part of constraints.
	 */
	vector<SmartPtr<Function> > conNonQuad;
	/** Sparsity graph of nonquadratic part of constraints.
	 */
	vector<SmartPtr<SparsityGraph> > conNonQuadSparsityGraph;
	// TODO: something for discrete var
	
	/** The number of original variables, i.e., variables that appear in MINLPData.
	 */
//	int nr_orig_var;
//	int nr_aux_var;
	
	/** Gives for each block function in each constraint the aux. var. that was introduced for it in the MIQQP.
	 * -1 if no auxvar was needed.
	 */ 
	vector<vector<int> > auxvar;
	
	/** Gives for each aux. var. that was introduced the constraint and block function index which it represents.
	 */
	vector<pair<int,int> > blockfunc;
	
	map<pair<int,int>, int> sparsity_hessian;
	int nnz_jac;



	/** Constructs the quadratic approximation from a MINLPData object.
	 */
	void construct();
	
	void addQuadEstConstraint(int con_nr, SymSparseMatrix& A, DenseVector& alpha, vector<int>& indices, int auxvar_index);
	void addQuadEstConstraint(int con_nr, QuadraticFunction& quad, DenseVector* alpha, vector<int>& indices, int auxvar_index);  
	void addConvexificationTerm(SymSparseMatrixCreator& A, SparseVectorCreator& b, double& constant, const DenseVector& alpha, const vector<int>& indices);
	/** Initialize sparsity_hessian.
	 */
	void initSparsityStructures();
	
	void addToHessian(double* values, double factor, SymSparseMatrix& A);
	void addToHessian(double* values, double factor, SymSparseMatrixCreator& A);
		
public:
	QuadraticOrConvexApproximation(MINLPData& data_, bool use_convex_);
	
	int numVariables() const { return data.numVariables()+blockfunc.size(); }
	int numConstraints() const { return conQuad.size(); }
	
	void print(ostream& out) const;	

	bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);

	bool get_variables_types(Index n, VariableType* var_types);

	bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);

	bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);

	bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);

	bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
	
	bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values);
        
	bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values);

	bool eval_gi(Index n, const Number* x, bool new_x, Index i, Number& gi);

	bool eval_grad_gi(Index n, const Number* x, bool new_x, Index i, Index& nele_grad_gi, Index* jCol, Number* values);

	void finalize_solution(Bonmin::TMINLP::SolverReturn status, Index n, const Number* x, Number obj_value);

	const Bonmin::TMINLP::BranchingInfo* branchingInfo() const;

	const Bonmin::TMINLP::SosInfo* sosConstraints() const;
	
}; // class QuadraticOrConvexApproximation 
	
} // namespace LaGO

#endif /*LAGOQUADRATICORCONVEXAPPROXIMATION_HPP_*/
