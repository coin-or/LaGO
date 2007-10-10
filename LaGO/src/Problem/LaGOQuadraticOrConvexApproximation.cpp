// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOQuadraticOrConvexApproximation.hpp"
#include "LaGOSimpleBlockFunction.hpp"
#include "LaGOSumFunction.hpp"

namespace LaGO {

QuadraticOrConvexApproximation::QuadraticOrConvexApproximation(MINLPData& data_, bool use_convex_)
: data(data_), use_convex(use_convex_), solver_return(Bonmin::TMINLP::MINLP_ERROR), solution_objective(-getInfinity())
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
	conNonQuadSparsityGraph.resize(data.numConstraints()+nr_aux_con);
	
	// orig. variable bounds
	for (int i=0; i<data.numVariables(); ++i) {
		var_lb[i]=data.getVariable(i).getLower();
		var_ub[i]=data.getVariable(i).getUpper();
//		clog << "bounds var. " << i << ": " << var_lb[i] << '\t' << var_ub[i] << endl; 
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
		
		SmartPtr<SparsityGraph> nonquad_sparsitygraph;
		
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
						assert(IsValid(func.sparsitygraph));
						conNonQuadSparsityGraph[aux_con_index]=new SparsityGraph(*func.sparsitygraph, func.indices);
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
						assert(IsValid(func.sparsitygraph));
						conNonQuadSparsityGraph[aux_con_index]=new SparsityGraph(*func.sparsitygraph, func.indices);
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
					assert(IsValid(func.sparsitygraph));
					if (nonquad_length==1)
						nonquad_sparsitygraph=new SparsityGraph(*func.sparsitygraph, func.indices);
					else
						nonquad_sparsitygraph->add(*func.sparsitygraph, func.indices);
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
			conNonQuadSparsityGraph[c]=nonquad_sparsitygraph;
			con_lb[c]=data.getConstraint(c).lower;				
			con_ub[c]=data.getConstraint(c).upper;				
		} else { // objective
			objQuad=quadfunc;
			if (nonquad_length>1)
				objNonQuad=new SumFunction(nonquad);
			else if (nonquad_length==1)
				objNonQuad=nonquad.front();
			objNonQuadSparsityGraph=nonquad_sparsitygraph;
		}
	}
	
	initSparsityStructures();
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
		b[orig_index]+=alpha[i]*(data.getVariable(orig_index).getUpper()+data.getVariable(orig_index).getLower());
		constant-=alpha[i]*data.getVariable(orig_index).getLower()*data.getVariable(orig_index).getUpper();
	}
}

void QuadraticOrConvexApproximation::getSolution(DenseVector& x) const {
	x=solution;
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
	if (IsValid(objQuad)) out << *objQuad;
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
	
	out << "Number of nonzeros in Jacobian: " << nnz_jac << endl;
	out << "Number of nonerzos in Hessian: " << sparsity_hessian.size() << endl;
	out << "Hessian sparsity elements: ";
	for (map<pair<int,int>, int>::const_iterator it(sparsity_hessian.begin()); it!=sparsity_hessian.end(); ++it)
		out << '(' << it->first.first << ',' << it->first.second << ") ";
	out << endl;	 
}

void QuadraticOrConvexApproximation::initSparsityStructures() {
	sparsity_hessian.clear();
	nnz_jac=0;
	
	for (int i=0; i<objQuad->A->getNumNonzeros(); ++i)
		sparsity_hessian.insert(pair<pair<int,int>, int>(pair<int,int>(objQuad->A->getRowIndices()[i], objQuad->A->getColIndices()[i]), 0)); 
	if (IsValid(objNonQuad)) {
		for (SparsityGraph::const_arc_iterator it_arc(objNonQuadSparsityGraph->arc_begin()); it_arc!=objNonQuadSparsityGraph->arc_end(); ++it_arc) {
			int tail=(**(*it_arc).tail()).varindex;
			int head=(**(*it_arc).head()).varindex;
			if (tail<head)
				sparsity_hessian.insert(pair<pair<int,int>, int>(pair<int,int>(head, tail), 0));
			else
				sparsity_hessian.insert(pair<pair<int,int>, int>(pair<int,int>(tail, head), 0));
		}
	}

	for (int c=0; c<numConstraints(); ++c) {
//		clog << "Con. " << c << ": ";
		assert(conQuad[c]->haveSparsity());
		nnz_jac+=conQuad[c]->getSparsity().size();
//		if (conQuad[c]->A->getNumNonzeros()) clog << "\t quad: ";
		for (int i=0; i<conQuad[c]->A->getNumNonzeros(); ++i) {
			sparsity_hessian.insert(pair<pair<int,int>, int>(pair<int,int>(conQuad[c]->A->getRowIndices()[i], conQuad[c]->A->getColIndices()[i]), 0));
//			clog << conQuad[c]->A->getRowIndices()[i] << ',' << conQuad[c]->A->getColIndices()[i] << ' ';
		}
		if (IsValid(conNonQuad[c])) {
//			clog << "\t nonquad: ";
			if (!conNonQuad[c]->haveSparsity()) {
				cerr << "Do not have sparsity for nonquad. func. in con. " << c << endl;
				cerr << *conNonQuad[c];
			}
			assert(conNonQuad[c]->haveSparsity());
			// let's hope that this is also ok if the sparsity of the quadratic and the nonquad. part overlap
			nnz_jac+=conNonQuad[c]->getSparsity().size();
			for (SparsityGraph::const_arc_iterator it_arc(conNonQuadSparsityGraph[c]->arc_begin()); it_arc!=conNonQuadSparsityGraph[c]->arc_end(); ++it_arc) {
				int tail=(**(*it_arc).tail()).varindex;
				int head=(**(*it_arc).head()).varindex;
				if (tail<head)
					sparsity_hessian.insert(pair<pair<int,int>, int>(pair<int,int>(head, tail), 0));
				else
					sparsity_hessian.insert(pair<pair<int,int>, int>(pair<int,int>(tail, head), 0));
//				clog << tail << ',' << head << ' ';
			}
		}
//		clog << endl;
	}
	
}

bool QuadraticOrConvexApproximation::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style) {
	n=numVariables();
	m=numConstraints();
	nnz_jac_g=nnz_jac;
	nnz_h_lag=sparsity_hessian.size();
	index_style=Ipopt::TNLP::C_STYLE;
	return true;
}

bool QuadraticOrConvexApproximation::get_variables_types(Index n, VariableType* var_types) {
	for (int i=0; i<n; ++i)
		if (i<data.numVariables() && data.getVariable(i).isDiscrete())
			var_types[i]=Bonmin::TMINLP::INTEGER;
		else
			var_types[i]=Bonmin::TMINLP::CONTINUOUS;
	return true;
}

bool QuadraticOrConvexApproximation::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types) {
	for (int i=0; i<m; ++i)
		const_types[i]=Ipopt::TNLP::NON_LINEAR;
	return true;
}

bool QuadraticOrConvexApproximation::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
	CoinCopyN(var_lb.getElements(), n, x_l);
	CoinCopyN(var_ub.getElements(), n, x_u);
	CoinCopyN(con_lb.getElements(), m, g_l);
	CoinCopyN(con_ub.getElements(), m, g_u);
	return true;
}

bool QuadraticOrConvexApproximation::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
	if (init_x) {
		if (!data.getStartPoints().empty())
			CoinCopyN(data.getStartPoints().front().getElements(), data.numVariables(), x);
		else {
			CoinCopyN(var_lb.getElements(), data.numVariables(), x);
			CoinZero(x+data.numVariables(), x+n);
		}
	}
	assert(init_z==false);
	assert(init_lambda==false);
	return true;
}

bool QuadraticOrConvexApproximation::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
	DenseVector x_(n, x);
	obj_value=objQuad->eval(x_);
	if (IsValid(objNonQuad))
		obj_value+=objNonQuad->eval(x_);
	return true;
}

bool QuadraticOrConvexApproximation::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
	DenseVector x_(n, x);
	DenseVector grad(n);
	objQuad->gradient(grad, x_);
	CoinCopyN(grad.getElements(), n, grad_f);
	if (IsValid(objNonQuad)) {
		objNonQuad->gradient(grad, x_);
		for (int i=0; i<n; ++i, ++grad_f)
			*grad_f+=grad[i];
	}
	return true;
}

bool QuadraticOrConvexApproximation::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
	DenseVector x_(n, x);
	for (int c=0; c<m; ++c) {
		g[c]=conQuad[c]->eval(x_);
		if (IsValid(conNonQuad[c]))
			g[c]+=conNonQuad[c]->eval(x_);
	}
	return true;
}

bool QuadraticOrConvexApproximation::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values) {
	if (values==NULL) {
		for (int c=0; c<numConstraints(); ++c) {
			const vector<int>& sparsity(conQuad[c]->getSparsity());
			for (int i=0; i<(int)sparsity.size(); ++i, ++iRow, ++jCol) {
				*iRow=c;
				*jCol=sparsity[i];
			}
			if (IsValid(conNonQuad[c])) {
				const vector<int>& sparsity(conNonQuad[c]->getSparsity());
				for (int i=0; i<(int)sparsity.size(); ++i, ++iRow, ++jCol) {
					*iRow=c;
					*jCol=sparsity[i];
				}
			}
		}		
	} else {
		DenseVector grad(n);
		DenseVector x_(n, x);
		for (int c=0; c<numConstraints(); ++c) {
			conQuad[c]->gradient(grad, x_);
			const vector<int>& sparsity(conQuad[c]->getSparsity());
			for (int i=0; i<(int)sparsity.size(); ++i, ++values)
				*values=grad[sparsity[i]];
			if (IsValid(conNonQuad[c])) {
				conNonQuad[c]->gradient(grad, x_);
				const vector<int>& sparsity(conNonQuad[c]->getSparsity());
				for (int i=0; i<(int)sparsity.size(); ++i, ++values)
					*values=grad[sparsity[i]];
			}
		}
	}
	return true;
}

bool QuadraticOrConvexApproximation::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values) {
	if (values==NULL) {
		int k=0;
		for (map<pair<int,int>, int>::iterator it(sparsity_hessian.begin()); it!=sparsity_hessian.end(); ++it, ++iRow, ++jCol, ++k) {
			*iRow=it->first.first;
			*jCol=it->first.second;
			it->second=k;
		}
		assert(k==nele_hess);
	} else {
		CoinZeroN(values, nele_hess);
		DenseVector x_(n, x);
		if (obj_factor) { 
			addToHessian(values, 2*obj_factor, *objQuad->A);
			if (IsValid(objNonQuad)) {
				SymSparseMatrixCreator mat;
				objNonQuad->fullHessian(mat, x_);
				addToHessian(values, obj_factor, mat);
			}
		}
		for (int c=0; c<numConstraints(); ++c) {
			if (!lambda[c]) continue;
			addToHessian(values, 2*lambda[c], *conQuad[c]->A);
			if (IsValid(conNonQuad[c])) {
				SymSparseMatrixCreator mat;
				conNonQuad[c]->fullHessian(mat, x_);
				addToHessian(values, lambda[c], mat);
			}
		}
	}
	
	return true;
}

void QuadraticOrConvexApproximation::addToHessian(double* values, double factor, SymSparseMatrix& A) {
	for (int i=0; i<A.getNumNonzeros(); ++i) {
		map<pair<int,int>, int>::iterator it(sparsity_hessian.find(pair<int,int>(A.getRowIndices()[i], A.getColIndices()[i])));
		assert(it!=sparsity_hessian.end());
		values[it->second]=factor*A.getValues()[i];
	}
}

void QuadraticOrConvexApproximation::addToHessian(double* values, double factor, SymSparseMatrixCreator& A) {
	for (SymSparseMatrixCreator::iterator it_A(A.begin()); it_A!=A.end(); ++it_A) {
		map<pair<int,int>, int>::iterator it(sparsity_hessian.find(it_A->first));
		assert(it!=sparsity_hessian.end());
		values[it->second]=factor*it_A->second;
	}
}

bool QuadraticOrConvexApproximation::eval_gi(Index n, const Number* x, bool new_x, Index i, Number& gi) {
	DenseVector x_(n, x);
	gi=conQuad[i]->eval(x_);
	if (IsValid(conNonQuad[i]))
		gi+=conNonQuad[i]->eval(x_);
	return true;
}

bool QuadraticOrConvexApproximation::eval_grad_gi(Index n, const Number* x, bool new_x, Index i, Index& nele_grad_gi, Index* jCol, Number* values) {
	if (values==NULL) {
		const vector<int>& sparsity(conQuad[i]->getSparsity());
		for (int j=0; j<(int)sparsity.size(); ++j, ++jCol)
			*jCol=sparsity[j];
		if (IsValid(conNonQuad[i])) {
			const vector<int>& sparsity(conNonQuad[i]->getSparsity());
			for (int j=0; j<(int)sparsity.size(); ++j, ++jCol)
				*jCol=sparsity[j];
		}		
	} else {
		DenseVector grad(n);
		DenseVector x_(n, x);
		conQuad[i]->gradient(grad, x_);
		const vector<int>& sparsity(conQuad[i]->getSparsity());
		for (int j=0; j<(int)sparsity.size(); ++j, ++values)
			*values=grad[sparsity[j]];
		if (IsValid(conNonQuad[i])) {
			conNonQuad[i]->gradient(grad, x_);
			const vector<int>& sparsity(conNonQuad[i]->getSparsity());
			for (int j=0; j<(int)sparsity.size(); ++j, ++values)
				*values=grad[sparsity[j]];
		}
	}
	return true;
}

void QuadraticOrConvexApproximation::finalize_solution(Bonmin::TMINLP::SolverReturn status, Index n, const Number* x, Number obj_value) {
	solver_return=status;
	if (x) solution.setVector(data.numVariables(), x);
	solution_objective=obj_value;
}

} // namespace LaGO
