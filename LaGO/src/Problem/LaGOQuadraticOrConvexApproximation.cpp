// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOQuadraticOrConvexApproximation.hpp"
#include "LaGOSimpleBlockFunction.hpp"
#include "LaGOSumFunction.hpp"

namespace LaGO {

QuadraticOrConvexApproximation::QuadraticOrConvexApproximation(MINLPData& data_, bool use_convex_)
: data(data_), use_convex(use_convex_)
{ construct();	
}

void QuadraticOrConvexApproximation::construct() {
//	int nr_orig_var=data.numVariables();
	
	// setup auxiliary variables
	// one for each nonconvex nonquadratic block or convexified quadratic block if use_convex=true
	int nr_aux_var=0;
	// count auxiliary constraints
	// one for each under- and overestimator
	int nr_aux_con=0;
	auxvar.resize(data.numConstraints()+1);
	
	for (int c=0; c<=data.numConstraints(); ++c) {
		const MINLPData::ObjCon& objcon(c<data.numConstraints() ? (const MINLPData::ObjCon&)data.getConstraint(c) : (const MINLPData::ObjCon&)data.getObjective());
		auxvar[c].resize(objcon.decompfuncNL.size(), -1);
		bool need_lower = (c==data.numConstraints() || data.getConstraint(c).upper<+getInfinity());
		bool need_upper = (c <data.numConstraints() && data.getConstraint(c).lower>-getInfinity());	

		for (unsigned int k=0; k<objcon.decompfuncNL.size(); ++k) {
			BlockFunction& func(*objcon.decompfuncNL[k]);
			
			assert(IsNull(func.nonquad) || IsNull(func.quad)); // maybe we relax this later
			
			if ((!use_convex) && IsNull(func.nonquad)) continue; 
			
			if ( (need_lower && !(func.curvature&CONVEX)) || (need_upper && !(func.curvature&CONCAVE)) ) {
//				clog << "con " << objcon.name << " block " << k << ": " << func.curvature << endl;
//				clog << func.alpha_convexify.size() << '\t' << func.alpha_concavify.size() << '\t' << func.underestimators.size() << '\t' << func.overestimators.size() << endl;
				assert(func.alpha_convexify.size() || func.alpha_concavify.size() || !func.underestimators.empty() || !func.overestimators.empty()); 
				auxvar[c][k]=nr_aux_var;
				blockfunc.push_back(pair<int,int>(c,k));
				++nr_aux_var;
			}
			
			if (use_convex && (func.alpha_convexify.size() || func.alpha_concavify.size()))
				nr_aux_con+=(need_lower?1:0)+(need_upper?1:0);
			nr_aux_con+=func.underestimators.size();
			nr_aux_con+=func.overestimators.size();
			if (func.underestimators.empty() && !func.overestimators.empty() && need_lower)
				++nr_aux_con; // we will need to add func.nonquad as underestimator
			if (!func.underestimators.empty() && func.overestimators.empty() && need_upper)
				++nr_aux_con; // we will need to add func.nonquad as overestimator
		}
	}
	clog << "Number of aux. variables: " << nr_aux_var << "\t Number of aux. constraints: " << nr_aux_con << endl;
	
	var_lb.resize(data.numVariables()+nr_aux_var, -getInfinity());
	var_ub.resize(data.numVariables()+nr_aux_var, +getInfinity());
	con_lb.resize(data.numConstraints()+nr_aux_con, -getInfinity());
	con_ub.resize(data.numConstraints()+nr_aux_con, +getInfinity());
	conQuad.resize(data.numConstraints()+nr_aux_con);
	conNonQuad.resize(data.numConstraints()+nr_aux_con);
	
	// orig. variable bounds
	for (int i=0; i<data.numVariables(); ++i) {
		var_lb[i]=data.getVariable(i).getLower();
		var_ub[i]=data.getVariable(i).getUpper();
	}
	
	int aux_con_index=data.numConstraints(); // index of next aux. constraint
	// setup objective and constraints
	for (int c=0; c<=data.numConstraints(); ++c) {
		const MINLPData::ObjCon& objcon(c<data.numConstraints() ? (const MINLPData::ObjCon&)data.getConstraint(c) : (const MINLPData::ObjCon&)data.getObjective());
		bool need_lower = (c==data.numConstraints() || data.getConstraint(c).upper<+getInfinity());
		bool need_upper = (c <data.numConstraints() && data.getConstraint(c).lower>-getInfinity());	

		double constant=objcon.decompfuncConstant;
		
		SparseVectorCreator lin;
		if (IsValid(objcon.decompfuncLin)) lin.add(*objcon.decompfuncLin);
		
		SymSparseMatrixCreator quad(numVariables());
		
		list<SmartPtr<Function> > nonquad;
		int nonquad_length=0;
		
		for (unsigned int k=0; k<objcon.decompfuncNL.size(); ++k) {
			BlockFunction& func(*objcon.decompfuncNL[k]);
			
			if (auxvar[c][k]>=0) { // use aux. variable for this block func.
				lin.insert(data.numVariables()+auxvar[c][k], 1.);
				
				if (IsValid(func.quad)) {
					assert(use_convex); // otherwise why should we have added an aux.var.?
					if (need_lower) {
						addQuadEstConstraint(aux_con_index, *func.quad, func.alpha_convexify, func.indices, auxvar[c][k]);
						con_ub[aux_con_index]=0.;
						++aux_con_index;						
					}
					if (need_upper) {
						addQuadEstConstraint(aux_con_index, *func.quad, func.alpha_concavify, func.indices, auxvar[c][k]);
						con_lb[aux_con_index]=0.;
						++aux_con_index;
					}					
				} else {
					assert(IsValid(func.nonquad));
					for (list<SmartPtr<QuadraticEstimator> >::iterator it(func.underestimators.begin()); it!=func.underestimators.end(); ++it) {
						addQuadEstConstraint(aux_con_index, *(**it).func, use_convex ? &(**it).alpha : NULL, func.indices, auxvar[c][k]); 
						con_ub[aux_con_index]=0.;
						++aux_con_index;						
					}
					if (need_lower && func.underestimators.empty()) { 
						assert(func.curvature&CONVEX); // we would need underestimators if the function isn't convex
						conQuad[aux_con_index]=new QuadraticFunction(new SymSparseMatrix(numVariables()), new SparseVector(auxvar[c][k]+data.numVariables(), -1.), 0.);
						conNonQuad[aux_con_index]=new SimpleBlockFunction(GetRawPtr(func.nonquad), func.indices);
						con_ub[aux_con_index]=0.;
						++aux_con_index;
					}
					
					for (list<SmartPtr<QuadraticEstimator> >::iterator it(func.overestimators.begin()); it!=func.overestimators.end(); ++it) {
						addQuadEstConstraint(aux_con_index, *(**it).func, use_convex ? &(**it).alpha : NULL, func.indices, auxvar[c][k]); 
						con_lb[aux_con_index]=0.;
						++aux_con_index;						
					}
					if (need_upper && func.overestimators.empty()) { 
						assert(func.curvature&CONCAVE); // we would need overestimators if the function isn't concave
						conQuad[aux_con_index]=new QuadraticFunction(new SymSparseMatrix(numVariables()), new SparseVector(auxvar[c][k]+data.numVariables(), -1.), 0.);
						conNonQuad[aux_con_index]=new SimpleBlockFunction(GetRawPtr(func.nonquad), func.indices);
						con_lb[aux_con_index]=0.;
						++aux_con_index;
					}
				}				
			} else {
				if (IsValid(func.quad)) {
					quad.addBlockMatrix(1., *func.quad, func.indices);
					if (use_convex) {
						assert(!(func.alpha_convexify.size() && func.alpha_concavify.size()));
						if (func.alpha_convexify.size())
							addConvexificationTerm(quad, lin, constant, func.alpha_convexify, func.indices);
						else if (func.alpha_concavify.size())
							addConvexificationTerm(quad, lin, constant, func.alpha_concavify, func.indices);
					}					
				}
				if (IsValid(func.nonquad)) {
					nonquad.push_back(new SimpleBlockFunction(GetRawPtr(func.nonquad), func.indices));
					++nonquad_length;
				}
			}
		}
		
		QuadraticFunction* quadfunc=new QuadraticFunction(new SymSparseMatrix(quad), new SparseVector(lin), constant);
		if (c<data.numConstraints()) { // constraint
			conQuad[c]=quadfunc;
			if (nonquad_length>1)
				conNonQuad[c]=new SumFunction(nonquad);
			else if (nonquad_length==1)
				conNonQuad[c]=nonquad.front();
			con_lb[c]=data.getConstraint(c).lower;				
			con_ub[c]=data.getConstraint(c).upper;				
		} else { // objective
			objQuad=quadfunc;
			if (nonquad_length>1)
				objNonQuad=new SumFunction(nonquad);
			else if (nonquad_length==1)
				objNonQuad=nonquad.front();	
		}
	}	
}

// adds the constraint xAx + alpha*(x-lb x)(x-ub x) - t, lower and upper bound are set somewhere else 
void QuadraticOrConvexApproximation::addQuadEstConstraint(int con_nr, SymSparseMatrix& A, DenseVector& alpha, vector<int>& indices, int auxvar_index) {  
	assert(con_nr<(int)conQuad.size());  
	SymSparseMatrixCreator myA(numVariables());
	SparseVectorCreator b;
	double constant=0.;
	
	myA.addBlockMatrix(1., A, indices);
	b.insert(data.numVariables()+auxvar_index, -1);
	if (alpha.size())
		addConvexificationTerm(myA, b, constant, alpha, indices);
		
	conQuad[con_nr]=new QuadraticFunction(new SymSparseMatrix(myA), new SparseVector(b), constant);
}

void QuadraticOrConvexApproximation::addQuadEstConstraint(int con_nr, QuadraticFunction& quad, DenseVector* alpha, vector<int>& indices, int auxvar_index) {
	assert(con_nr<(int)conQuad.size());  
	SymSparseMatrixCreator A(numVariables());
	SparseVectorCreator b;
	double constant=quad.constant;
	
	if (IsValid(quad.A))
		A.addBlockMatrix(1., *quad.A, indices);
	b.insert(data.numVariables()+auxvar_index, -1);
	if (IsValid(quad.b))
		b.addBlockVector(*quad.b, indices);
	if (alpha)
		addConvexificationTerm(A, b, constant, *alpha, indices);
		
	conQuad[con_nr]=new QuadraticFunction(new SymSparseMatrix(A), new SparseVector(b), constant);
}


void QuadraticOrConvexApproximation::addConvexificationTerm(SymSparseMatrixCreator& A, SparseVectorCreator& b, double& constant, const DenseVector& alpha, const vector<int>& indices) {
	for (int i=0; i<alpha.size(); ++i) {
		if (alpha[i]==0.) continue;
		int orig_index=indices[i];
		
//		clog << "addConvexTerm: " << orig_index << "\t " << alpha[i] << ' ' << data.getVariable(orig_index).getLower() << ' ' << data.getVariable(orig_index).getUpper() << endl; 
		
		A[pair<int,int>(orig_index, orig_index)]-=alpha[i];
		b[orig_index]+=alpha[i]*(data.getVariable(orig_index).getUpper()-data.getVariable(orig_index).getLower());
		constant-=alpha[i]*data.getVariable(orig_index).getLower()*data.getVariable(orig_index).getUpper();
	}
}


void QuadraticOrConvexApproximation::print(ostream& out) const {
	out << "Variables: " << numVariables() << "\t Constraints: " << numConstraints();
	for (int i=0; i<data.numVariables(); ++i)
		out << data.getVariable(i); 
	for (int i=0; i<(int)blockfunc.size(); ++i) {
		out << data.numVariables()+i << ": Aux.var. " << i << " is for con. " << blockfunc[i].first << " (";
		if (blockfunc[i].first<data.numConstraints()) out << data.getConstraint(blockfunc[i].first).name;
		else out << data.getObjective().name;
		out << ") block " << blockfunc[i].second << endl;
	}
	
	out << "Objective: ";
	if (IsValid(objQuad)) out << "quad: " << *objQuad;
	if (IsValid(objNonQuad)) out << "nonquad: " << *objNonQuad;
	out << endl;
	
	for (int c=0; c<numConstraints(); ++c) {
		if (c<data.numConstraints())
			out << data.getConstraint(c).name;
		else
			out << "aux. constraint " << c-data.numConstraints();
		out << ": [" << con_lb[c] << ", " << con_ub[c] << "] ";
		if (IsValid(conQuad[c])) out << "quad: " << *conQuad[c];
		if (IsValid(conNonQuad[c])) out << "nonquad: " << *conNonQuad[c];
		out << endl; 		
	}	
}
	
} // namespace LaGO
