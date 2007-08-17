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
	/** Quadratic (and linear and constant) part of constraints.
	 */
	vector<SmartPtr<QuadraticFunction> > conQuad;
	/** Convex nonquadratic part of constraints.
	 */
	vector<SmartPtr<Function> > conNonQuad;
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
	
	/** Constructs the quadratic approximation from a MINLPData object.
	 */
	void construct();
	
	void addQuadEstConstraint(int con_nr, SymSparseMatrix& A, DenseVector& alpha, vector<int>& indices, int auxvar_index);
	void addQuadEstConstraint(int con_nr, QuadraticFunction& quad, DenseVector* alpha, vector<int>& indices, int auxvar_index);  
	void addConvexificationTerm(SymSparseMatrixCreator& A, SparseVectorCreator& b, double& constant, const DenseVector& alpha, const vector<int>& indices);
	
public:
	QuadraticOrConvexApproximation(MINLPData& data_, bool use_convex_);
	
	int numVariables() const { return data.numVariables()+blockfunc.size(); }
	int numConstraints() const { return conQuad.size(); }
	
	
	void print(ostream& out) const;	
	
}; // class QuadraticOrConvexApproximation 
	
} // namespace LaGO

#endif /*LAGOQUADRATICORCONVEXAPPROXIMATION_HPP_*/
