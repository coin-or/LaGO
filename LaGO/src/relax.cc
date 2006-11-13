// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "relax.h"


// ---------------------------------------- Convexify -------------------------------------------------

void Convexify::new_sampleset(UserVector<double>& lower, UserVector<double>& upper, SepQcFunc& f) {
	sample_set.clear();
	sample_set.resize(f.block.size());

	for (int k=0; k<f.block.size(); k++)
		if (f.s[k] || f.A[k])
			sampling.get_points(sample_set[k], lower(f.block[k]), upper(f.block[k]));
}


double Convexify::convexify(SepQcFunc& f, bool eq, vector<double>& min_eigval, vector<double>& max_eigval, UserVector<double>& lower, UserVector<double>& upper) {
	assert(sample_set.size());

	convexify_c.resize(f.block.size());
	bool try_shift_for_s=true;
	for (int k=0; try_shift_for_s && k<f.block.size(); k++) {
		if (f.A[k] || f.s[k]) { // check, if min-eig * max-eig is high enough to try cgu-convexification
  		try_shift_for_s&=(min_eigval[k] * max_eigval[k] >= -very_large);
			assert(sample_set[k].size());
		}
  	convexify_c[k].first=0.;
  	convexify_c[k].second=0.;
 	}

	convexify_a.first=NULL;
	convexify_b.first=NULL;
	convexify_a.second=NULL;
	convexify_b.second=NULL;
 	if (!(f.get_curvature()&Func::CONVEX)) {
 		convexify_a.first=new dvector(f.dim());
	 	convexify_b.first=new dvector(f.dim());
	}
 	if (eq && !(f.get_curvature()&Func::CONCAVE)) {
	 	convexify_a.second=new dvector(f.dim());
 		convexify_b.second=new dvector(f.dim());
	}
	characteristica1.clear();
	characteristica2.clear();	

 	convexified_by_cgu.resize(f.block.size());

 	pair<double, double> alpha_norm_sq;
	double max_alpha_norm_sq=0.;
 	bool rerun_with_cgu=false;
 	for (int k=0; k<f.block.size(); k++) {
		if (f.s[k]) {
			if (try_shift_for_s && (!rerun_with_cgu)) {
				alpha_norm_sq=convexify_alpha2(f, k, lower(f.block[k]), upper(f.block[k]));
				convexified_by_cgu[k]=false;
				if (alpha_norm_sq.first > very_large*very_large) rerun_with_cgu=true;
				else if (alpha_norm_sq.second > very_large*very_large) rerun_with_cgu=true;
				else max_alpha_norm_sq=MAX(max_alpha_norm_sq, MAX(alpha_norm_sq.first, alpha_norm_sq.second));
			}
 			else {
 				convexify_cgu(f, k);
 				convexified_by_cgu[k]=true;
 			}
		}	else if (f.A[k]) {
 			if (min_eigval[k] * max_eigval[k] < -very_large) {
 				convexify_cgu(f, k);
 				convexified_by_cgu[k]=true;
 			}
 			else {
 				alpha_norm_sq=convexify_alpha2(f, k, lower(f.block[k]), upper(f.block[k]));
 				if ((alpha_norm_sq.first > very_large*very_large) || (alpha_norm_sq.second > very_large*very_large)) {
					convexify_cgu(f, k);
					convexified_by_cgu[k]=true;
	 			} else {
					convexified_by_cgu[k]=false;
					max_alpha_norm_sq=MAX(max_alpha_norm_sq, MAX(alpha_norm_sq.first, alpha_norm_sq.second));
				}
			}
		}
	}

  if (rerun_with_cgu)
  	for (int k=0; k<f.block.size(); k++)
  		if (f.s[k] && (!convexified_by_cgu[k])) {
  			convexify_cgu(f, k);
 				convexified_by_cgu[k]=true;
 			}

	return sqrt(max_alpha_norm_sq);
}



void Convexify::get_decomposed_functions(pair<Pointer<SepQcFunc>, Pointer<SepQcFunc> >& f) {
	assert(!f.second);

 	SepQcFunc& orig_f(*f.first);

 	if (convexify_a.first) f.first=new SepQcFunc(orig_f);
	if (convexify_a.second) f.second=new SepQcFunc(orig_f, true);

	for (int k=0; k<orig_f.block.size(); k++) {
		if (convexify_a.first) {
			Pointer<SparsityInfo> si(new SparsityInfo(2));

			Pointer<SparseMatrix2> A(new SparseMatrix2(orig_f.block[k].size()));
			for (int i=0; i<orig_f.block[k].size(); i++)
				A->AddElement(i, i, (*convexify_a.first)(orig_f.block[k][i]));
			A->finish();

 			Pointer<SparseVector<double> > B(new SparseVector<double>(*convexify_b.first, orig_f.block[k]));

			if (convexified_by_cgu[k]) {
				f.first->A[k]=NULL;
				f.first->b[k]=NULL;
				f.first->s[k]=NULL;
			}

			if (A->nonzeros())
				if (f.first->A[k]) f.first->A[k]=new SumMatrix(&*f.first->A[k], &*A);
				else f.first->A[k]=A;
			if (f.first->A[k]) si->add(*f.first->A[k]);
			
			if (!(*B==0))
				if (f.first->b[k]) {
					f.first->b[k]=f.first->b[k]->getcopy();
					*f.first->b[k]+=*B;
				}	else f.first->b[k]=B;
			if (f.first->b[k]) si->add(*f.first->b[k]);
			
			f.first->c+=convexify_c[k].first;

			if (f.first->s[k]) si->add(((const Func*)(Func*)f.first->s[k])->get_sparsity());

			f.first->set_curvature(k, (f.first->A[k] || f.first->s[k]) ? Func::CONVEX : Func::LINEAR);
			f.first->set_sparsity(k, si);
		}

		if (convexify_a.second) { // and the same fun for the second function
			Pointer<SparsityInfo> si(new SparsityInfo(2));

			Pointer<SparseMatrix2> A(new SparseMatrix2(orig_f.block[k].size()));
			for (int i=0; i<orig_f.block[k].size(); i++)
				A->AddElement(i, i, (*convexify_a.second)(orig_f.block[k][i]));
			A->finish();

 			Pointer<SparseVector<double> > B(new SparseVector<double>(*convexify_b.second, orig_f.block[k]));

			if (convexified_by_cgu[k]) {
				f.second->A[k]=NULL;
				f.second->b[k]=NULL;
				f.second->s[k]=NULL;
			} else
				f.second->A[k]=orig_f.A[k];

			if (A->nonzeros()) {
				if (f.second->A[k]) f.second->A[k]=new SumMatrix(&*f.second->A[k], &*A, -1.);
				else f.second->A[k]=A;
			}	else if (f.second->A[k]) f.second->A[k]=new MinusMatrix(f.second->A[k]);
			if (f.second->A[k]) si->add(*f.second->A[k]);

			if (!(*B==0))
				if (f.second->b[k]) *f.second->b[k]+=*B;
				else f.second->b[k]=B;
			if (f.second->b[k]) si->add(*f.second->b[k]);
			
			f.second->c+=convexify_c[k].second;

			if (f.second->s[k]) si->add(((const Func*)(Func*)f.second->s[k])->get_sparsity());

			f.second->set_curvature(k, (f.second->A[k] || f.second->s[k]) ? Func::CONVEX : Func::LINEAR);
			f.second->set_sparsity(k, si);
		}
	}
}

void Convexify::convexify_cgu(SepQcFunc& f, int k) {
	MinlpProblem lp(MinlpProblem::QQP);
	int i=0;
	for (; i<f.block[k].size(); i++)
		lp.add_var(i, 0, false, 0., INFINITY); // a's
	for (; i<2*f.block[k].size(); i++)
		lp.add_var(i, 0, false, -INFINITY, INFINITY); // b's
	lp.add_var(i, 0, false, -INFINITY, INFINITY, "c"); // c
	
	dvector start_point(2*f.block[k].size()+1);
		
	dvector xb(f.block[k].size()); // coefficients for a-variables
	dvector xa(f.block[k].size()); // coefficients for b-variables
	for (vector<dvector>::iterator it(sample_set[k].begin()); it!=sample_set[k].end(); it++) {
		xb+=*it;
		xa+=it->diagmult(*it);
	}
		
	Pointer<SepQcFunc> obj(new SepQcFunc(lp.block));
	obj->b[0]=new dvector(obj->dim());
	for (int j=0; j<f.block[k].size(); j++) {
		(*obj->b[0])[j]=-0.5*xa[j]+0.5;
		(*obj->b[0])[f.block[k].size()+j]=-xb[j];
	}
	(*obj->b[0])[2*f.block[k].size()]=-0.5*sample_set[k].size();
		
	lp.add_obj(obj);
		
	for (vector<dvector>::iterator it(sample_set[k].begin()); it!=sample_set[k].end(); it++) {
		Pointer<SepQcFunc> con(new SepQcFunc(lp.block));
		con->b[0]=new dvector(con->dim());
		for (int j=0; j<f.block[k].size(); j++) {
			(*con->b[0])[j]=0.5* ((*it)(j) * (*it)(j));
			(*con->b[0])[f.block[k].size()+j]=(*it)(j);
		}
		(*con->b[0])[2*f.block[k].size()]=0.5;
			
		con->c=-f.eval(*it, k)-cgu_eps;

		lp.add_con(con, false);
	}
	
	int ret;
	Pointer<SparseMatrix2> H1, H2;
	if (f.s[k])
		H1=new SparseMatrix2(HessMatrix(*f.s[k], sample_set[k][0]));
	if (f.A[k])
		H2=new SparseMatrix2(*f.A[k]);
	dvector grad(f.block[k].size());
	grad=f.grad(sample_set[k][0], k);
				
	for (int i=0; i<f.block[k].size(); i++) {
		if (H1) start_point[i]=.5*(*H1)(i,i);
		if (H2) start_point[i]+=.5*(*H2)(i,i);
		if (start_point[i]<0) start_point[i]=0;
			
		start_point[f.block[k].size()+i]=0.5*(grad[i]-start_point[i]*sample_set[k][0][i]);

		start_point[2*f.block[k].size()]-=sample_set[k][0][i]*sample_set[k][0][i]*start_point[i];
		start_point[2*f.block[k].size()]-=2*sample_set[k][0][i]*start_point[f.block[k].size()+i];
	}
	start_point[2*f.block[k].size()]+=f.eval(sample_set[k][0], k);

	if (convexify_a.first) {
		Pointer<LocOpt> lpsolver(LocOpt::get_lp_solver(Pointer<MinlpProblem>(&lp, false), param, "CGU", NULL, NULL));
		ret=lpsolver->solve(start_point);
//		out_log << "CGU block " << k << ": SnOpt: " << ret << ";\t ";
//		out_log << snopt.sol_point;

		convexify_a.first->set_block(lpsolver->sol_point(0, f.block[k].size()-1), f.block[k]);
		convexify_b.first->set_block(lpsolver->sol_point(f.block[k].size(), 2*f.block[k].size()), f.block[k]);
		convexify_c[k].first=lpsolver->sol_point(2*f.block[k].size());
	}

	if (convexify_a.second) {

		for (int i=0; i<lp.con.size(); i++) {
			lp.con[i]->c*=-1;
			lp.con[i]->c-=2*cgu_eps;
		}

		if (H1) *H1*=-1;
		if (H2) *H2*=-1;
		grad*=-1;
		start_point=0.;
		for (int i=0; i<f.block[k].size(); i++) {
			if (H1) start_point[i]=.5*(*H1)(i,i);
			if (H2) start_point[i]+=.5*(*H2)(i,i);
			if (start_point[i]<0) start_point[i]=0;

			start_point[f.block[k].size()+i]=0.5*(grad[i]-start_point[i]*sample_set[k][0][i]);

			start_point[2*f.block[k].size()]-=sample_set[k][0][i]*sample_set[k][0][i]*start_point[i];
			start_point[2*f.block[k].size()]-=2*sample_set[k][0][i]*start_point[f.block[k].size()+i];
		}
		start_point[2*f.block[k].size()]-=f.eval(sample_set[k][0], k);

		Pointer<LocOpt> lpsolver(LocOpt::get_lp_solver(Pointer<MinlpProblem>(&lp, false), param, "CGU", NULL, NULL));
		ret=lpsolver->solve(start_point);
//		out_log << "CGU -block " << k << ": SnOpt: " << ret << ";\t " << endl;
//		out_log << snopt.sol_point;

		convexify_a.second->set_block(lpsolver->sol_point(0, f.block[k].size()-1), f.block[k]);
		convexify_b.second->set_block(lpsolver->sol_point(f.block[k].size(), 2*f.block[k].size()), f.block[k]);
		convexify_c[k].second=lpsolver->sol_point(2*f.block[k].size());
	}
//	out_log << endl;
}

pair<double, double> Convexify::convexify_alpha(SepQcFunc& f, int k, const dvector& low, const dvector& up, Pointer<set<int> > i_quad, Pointer<set<int> > i_nonquadlin) {
	dvector diam(up-low);  // diameter

	double lowest_lambda1=0.; // lambda_min of inequality-constraint
	double lowest_lambda2=0.; // lambda_min of -inequality-constraint
	dvector eig_vec(f.block[k].size());

	dvector shift(f.block[k].size());
	for (int i=0; i<shift.dim(); i++)
		if ((i_quad && i_quad->count(i)) || (i_nonquadlin && i_nonquadlin->count(i)))
			if (diam(i)>rtol)
				shift[i]=1./(diam(i)*diam(i));

//	out_log << "shift block " << k << ": " << shift;

	SparseMatrix2 A(f.block[k].size());
	if (f.A[k]) {
		A+=*f.A[k];
		A*=2; // hessian of quadratic part
	}
	if (f.s[k]) { // nonquadratic function
		double lambda;
		for (vector<dvector>::iterator it(sample_set[k].begin()); it!=sample_set[k].end(); it++) {
			SparseMatrix2 H(A);
			H+=HessMatrix(*f.s[k], *it);
			H.finish();

			const int* row_ind=H.GetRowInd();
			const int* col_ptr=H.GetColPtr();
			double* val=H.GetVal();
			int j=0;  // ==col_ptr[0]
			for (int i=0; i<H.dim() && j<H.dim(); i++) // rows
			  for (; j<col_ptr[i+1]; j++, val++) // columns
			    *val*=diam(i)*diam(row_ind[j]);

			if (convexify_a.first) {
				int ret=H.eig(eig_vec, lambda, param);
				if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for nonquad. func block " << k << endl;
				if ((!ret) && lambda<lowest_lambda1) lowest_lambda1=lambda;
				eig_vec=0.;
			}
			if (convexify_a.second) {
				val=H.GetVal();
				for (int i=0; i<H.nonzeros(); i++, val++) *val*=-1;
				int ret=H.eig(eig_vec, lambda, param);
				if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for nonquad. func block " << k << endl;
				if ((!ret) && lambda<lowest_lambda2) lowest_lambda2=lambda;
				eig_vec=0.;
			}
		}
	} else if (f.A[k]) { // quadratic function
		A.finish();
		const int* row_ind=A.GetRowInd();
		const int* col_ptr=A.GetColPtr();
		double* val=A.GetVal();
		int j=0;
		for (int i=0; i<A.dim(); i++) // columns
			for (; j<col_ptr[i+1]; j++, val++) // rows
				*val *= diam(i) * diam(row_ind[j]);

		if (convexify_a.first) {
			int ret=A.eig(eig_vec, lowest_lambda1, param);
			if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for quad. func block " << k << endl;
		}

		if (convexify_a.second) {
			val=A.GetVal();
			for (int i=0; i<A.nonzeros(); i++, val++) *val*=-1;
			eig_vec=0.;
			int ret=A.eig(eig_vec, lowest_lambda2, param);
			if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for quad. func block " << k << endl;
		}
	}

	pair<double, double> alpha_norm_sq(0., 0.);

// to give some space
//	lowest_lambda1*=1.05;
//	lowest_lambda2*=1.05;

	// shift is now (-.5*lowest_lambda)*shift
	if (lowest_lambda1<-rtol) {
		alpha_norm_sq.first=.25*lowest_lambda1*lowest_lambda1*shift.sq_norm2();

//		out_log << "shift block " << k << " by " << -.5*lowest_lambda1 << ";\t ";
		// add shift (x-low) diag(shift) (x-up) = x diag(shift) x - diag(shift)(low+up) * x + low*diag(shift)*up
		convexify_a.first->set_block((-.5*lowest_lambda1)*shift, f.block[k]);
		convexify_b.first->set_block((.25*lowest_lambda1)*shift.diagmult(low+up), f.block[k]);
		convexify_c[k].first=-.5*lowest_lambda1*low*shift.diagmult(up);
	}
	if (lowest_lambda2<-rtol) {
		alpha_norm_sq.second=.25*lowest_lambda2*lowest_lambda2*shift.sq_norm2();

//		out_log << "shift -block " << k << " by " << -.5*lowest_lambda2 << ";\t ";
		convexify_a.second->set_block((-.5*lowest_lambda2)*shift, f.block[k]);
		convexify_b.second->set_block((.25*lowest_lambda2)*shift.diagmult(low+up), f.block[k]);
		convexify_c[k].second=-.5*lowest_lambda2*low*shift.diagmult(up);
	}
//	out_log << endl;
	return alpha_norm_sq;
}

pair<double, double> Convexify::convexify_alpha(SepQcFunc& f, int k, const dvector& low, const dvector& up) {
	dvector diam(up-low);  // diameter

	double lowest_lambda1=0.; // lambda_min of inequality-constraint
	double lowest_lambda2=0.; // lambda_min of -inequality-constraint
	dvector eig_vec(f.block[k].size());

	dvector shift(f.block[k].size());
	for (VariableIterator it(f.get_sparsity(k), false, true); it; ++it) // only nonlinear variables
		if (diam(it())>rtol)
			shift[it()]=1./(diam(it())*diam(it()));

//	out_log << "shift block " << k << ": " << shift;

	SparseMatrix2 A(f.block[k].size());
	if (f.A[k]) {
		A+=*f.A[k];
		A*=2; // hessian of quadratic part
	}
	if (f.s[k]) { // nonquadratic function
		double lambda;
		for (vector<dvector>::iterator it(sample_set[k].begin()); it!=sample_set[k].end(); it++) {
			SparseMatrix2 H(A);
			H+=HessMatrix(*f.s[k], *it);
			H.finish();

			const int* row_ind=H.GetRowInd();
			const int* col_ptr=H.GetColPtr();
			double* val=H.GetVal();
			int j=0;  // ==col_ptr[0]
			for (int i=0; i<H.dim() && j<H.dim(); i++) // rows
			  for (; j<col_ptr[i+1]; j++, val++) // columns
			    *val*=diam(i)*diam(row_ind[j]);
		
			if (convexify_a.first) {
				int ret=H.eig(eig_vec, lambda, param);
				if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for nonquad. func block " << k << endl;
				if ((!ret) && lambda<lowest_lambda1) lowest_lambda1=lambda;
				eig_vec=0.;
			}
			if (convexify_a.second) {
				val=H.GetVal();
				for (int i=0; i<H.nonzeros(); i++, val++) *val*=-1;
				int ret=H.eig(eig_vec, lambda, param);
				if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for nonquad. func block " << k << endl;
				if ((!ret) && lambda<lowest_lambda2) lowest_lambda2=lambda;
				eig_vec=0.;
			}
		}
	} else if (f.A[k]) { // quadratic function
		A.finish();
		const int* row_ind=A.GetRowInd();
		const int* col_ptr=A.GetColPtr();
		double* val=A.GetVal();
		int j=0;
		for (int i=0; i<A.dim(); i++) // columns
			for (; j<col_ptr[i+1]; j++, val++) // rows
				*val *= diam(i) * diam(row_ind[j]);

		if (convexify_a.first) {
			int ret=A.eig(eig_vec, lowest_lambda1, param);
			if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for quad. func block " << k << endl;
		}

		if (convexify_a.second) {
			val=A.GetVal();
			for (int i=0; i<A.nonzeros(); i++, val++) *val*=-1;
			eig_vec=0.;
			int ret=A.eig(eig_vec, lowest_lambda2, param);
			if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for quad. func block " << k << endl;
		}
	}

	pair<double, double> alpha_norm_sq(0., 0.);

// to give some space
//	lowest_lambda1*=1.05;
//	lowest_lambda2*=1.05;

	// shift is now (-.5*lowest_lambda)*shift
	if (lowest_lambda1<-rtol) {
		alpha_norm_sq.first=.25*lowest_lambda1*lowest_lambda1*shift.sq_norm2();

//		out_log << "shift block " << k << " by " << -.5*lowest_lambda1 << ";\t ";
		// add shift (x-low) diag(shift) (x-up) = x diag(shift) x - diag(shift)(low+up) * x + low*diag(shift)*up
		convexify_a.first->set_block((-.5*lowest_lambda1)*shift, f.block[k]);
		convexify_b.first->set_block((.25*lowest_lambda1)*shift.diagmult(low+up), f.block[k]);
		convexify_c[k].first=-.5*lowest_lambda1*low*shift.diagmult(up);
		for (unsigned int i=0; i<f.block[k].size(); ++i)
			if (shift(i)) characteristica1[f.block[k][i]]=-.5*lowest_lambda1*shift(i);
	}
	if (lowest_lambda2<-rtol) {
		alpha_norm_sq.second=.25*lowest_lambda2*lowest_lambda2*shift.sq_norm2();

//		out_log << "shift -block " << k << " by " << -.5*lowest_lambda2 << ";\t ";
		convexify_a.second->set_block((-.5*lowest_lambda2)*shift, f.block[k]);
		convexify_b.second->set_block((.25*lowest_lambda2)*shift.diagmult(low+up), f.block[k]);
		convexify_c[k].second=-.5*lowest_lambda2*low*shift.diagmult(up);
		for (unsigned int i=0; i<f.block[k].size(); ++i)
			if (shift(i)) characteristica2[f.block[k][i]]=-.5*lowest_lambda2*shift(i);
	}
//	out_log << endl;
	return alpha_norm_sq;
}

pair<double, double> Convexify::convexify_alpha2(SepQcFunc& f, int k, const dvector& low, const dvector& up) {
	if (f.s[k]) return convexify_alpha(f, k, low, up);

	pair<double, double> alpha_norm_sq(0., 0.);

	Pointer<SepQcFunc> fk(new SepQcFunc(f.A[k], NULL, NULL, 0., new SparsityInfo(f.get_sparsity(k))));
	Pointer<SepQcFunc> fkb(decomp.decompose((Pointer<Func>)fk, sample_set[k][0]));

	for (int l=0; l<fkb->block.size(); ++l) {
		if (!fkb->A[l]) continue;

		dvector shift(fkb->block[l].size());
		dvector diam(shift.dim()); // diameter
		for (int i=0; i<shift.dim(); ++i) {
			diam[i]=up(fkb->block[l][i])-low(fkb->block[l][i]);
			if (diam(i)>rtol) shift[i]=1./(diam(i)*diam(i));
		}

		SparseMatrix2 A(*fkb->A[l]);
		A.finish();
		const int* row_ind=A.GetRowInd();
		const int* col_ptr=A.GetColPtr();
		double* val=A.GetVal();
		int j=0;
		for (int i=0; i<A.dim(); ++i) // columns
			for (; j<col_ptr[i+1]; ++j, ++val) // rows
				*val *= 2.*diam(i) * diam(row_ind[j]);

		dvector eigvec(shift.dim());
		double eigval;
		if (convexify_a.first) {
			int ret=A.eig(eigvec, eigval, param);
			if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for quad. func block " << k << ',' << l << endl;

			if (eigval<-rtol) {
				alpha_norm_sq.first+=eigval*eigval*shift.sq_norm2();
				
//		out_log << "shift block " << k << ',' << l << " by " << -.5*eigval << ";\t ";
		// add shift (x-low) diag(shift) (x-up) = x diag(shift) x - diag(shift)(low+up) * x + low*diag(shift)*up
				for (int i=0; i<shift.dim(); ++i) {
					int ib=fkb->block[l][i];
					int i0=f.block[k][ib];
					(*convexify_a.first)[i0]=-.5*eigval*shift(i);
					(*convexify_b.first)[i0]=.25*eigval*shift(i)*(low(ib)+up(ib));
					convexify_c[k].first+=-.5*eigval*low(ib)*shift(i)*up(ib);
					characteristica1[i0]=-.5*eigval*shift(i);
				}
			}
		}

		if (convexify_a.second) {
			val=A.GetVal();
			for (int i=0; i<A.nonzeros(); i++, val++) *val*=-1;
			eigvec=0.;
			int ret=A.eig(eigvec, eigval, param);
			if (ret) out_out << "Convexify::convexify_alpha(): Error " << ret << " in eigval computation for quad. func block " << k << ',' << l << endl;

			if (eigval<-rtol) {
				alpha_norm_sq.second+=eigval*eigval*shift.sq_norm2();
				
//		out_log << "shift block " << k << ',' << l << " by " << -.5*eigval << ";\t ";
		// add shift (x-low) diag(shift) (x-up) = x diag(shift) x - diag(shift)(low+up) * x + low*diag(shift)*up
				for (int i=0; i<shift.dim(); ++i) {
					int ib=fkb->block[l][i];
					int i0=f.block[k][ib];
					(*convexify_a.second)[i0]=-.5*eigval*shift(i);
					(*convexify_b.second)[i0]=.25*eigval*shift(i)*(low(ib)+up(ib));
					convexify_c[k].second+=-.5*eigval*low(ib)*shift(i)*up(ib);
					characteristica2[i0]=-.5*eigval*shift(i);
				}
			}
		}
	}

	return alpha_norm_sq;
}


void Convexify::check_convex2(vector<double>& min_eigval, vector<double>& max_eigval, SepQcFunc& f, const vector<vector<dvector> >& sample_set) {
	min_eigval.resize(f.block.size());
	max_eigval.resize(f.block.size());
	int ret;

	for (int k=f.block.size()-1; k>=0; --k) {
		if (!(f.get_curvature(k)&Func::UNKNOWN)) continue; // curvature already known
		if ((!f.s[k]) && (!f.A[k])) { f.set_curvature(k, Func::LINEAR); continue; }

		dvector eigvec(f.block[k].size());
		if (!f.s[k]) { // quadratic function
			ret=f.A[k]->eig(eigvec, min_eigval[k], param);
			if (ret) {
				out_err << "Convexify::check_convex(): Error " << ret << " computing minimum eigval for quad. func. block " << k << "! Aborting." << endl;
				exit(-1);
			}
			min_eigval[k]*=2.;

			eigvec=0; // using eigvec for mineig as approximation for maxeig gives trouble for some examples
			MinusMatrix MinusA(f.A[k]);
			ret=MinusA.eig(eigvec, max_eigval[k], param);
			if (ret) {
				out_err << "Convexify::check_convex(): Error " << ret << " computing maximum eigval for quad. func. block " << k << "! Aborting." << endl;
				exit(-1);
			}
			max_eigval[k]*=-2.;
		} else {
			Pointer<SepQcFunc> fk(new SepQcFunc(f.A[k], NULL, f.s[k], 0., new SparsityInfo(f.get_sparsity(k))));
			Pointer<SepQcFunc> fkb(decomp.decompose((Pointer<Func>)fk, sample_set[k][0]));

			min_eigval[k]=INFINITY; max_eigval[k]=-INFINITY;
			bool ss_checked=false;
			for (int l=fkb->block.size()-1; l>=0; --l) {
				if ((!fkb->s[l]) && (!fkb->A[l])) continue;
				
				out_log << fkb->block[l].size();
				dvector eigvec(fkb->block[l].size()); double eigval;
				if (!fkb->s[l]) {
					ret=fkb->A[l]->eig(eigvec, eigval, param);
					if (ret) {
						out_err << "Convexify::check_convex(): Error " << ret << " computing eigval for quad. func. block " << k << ',' << l << "! Aborting." << endl;
						exit(-1);
					}
					eigval*=2.;
					if (eigval<min_eigval[k]) min_eigval[k]=eigval;

					eigvec=0; // using eigvec for mineig as approximation for maxeig gives trouble for some examples
					ret=MinusMatrix(fkb->A[l]).eig(eigvec, eigval, param);
					if (ret) {
						out_err << "Convexify::check_convex(): Error " << ret << " computing eigval for quad. func. block " << k << ',' << l << "! Aborting." << endl;
						exit(-1);
					}
					eigval*=-2.;
					if (eigval>max_eigval[k]) max_eigval[k]=eigval;
					continue;
				}

				bool onlyerror1=true; bool onlyerror2=true;
				for (vector<dvector>::const_iterator it(sample_set[k].begin()); it!=sample_set[k].end(); ++it) {
					if (!ss_checked) if (!finite(fkb->eval(*it))) continue; // skipping sample points with nan

					dvector x(*it, fkb->block[l]);

					Pointer<UserMatrix> H(new HessMatrix(*fkb->s[l], x));
					if (fkb->A[l]) ret=SumMatrix(H, fkb->A[l]).eig(eigvec, eigval, param);
					else ret=H->eig(eigvec, eigval, param);
					if (ret) out_out << "Convexify::check_convex(): Error " << ret << " computing eigval for nonquad. func. block " << k << ',' << l << endl;
					if ((!ret) && eigval<min_eigval[k]) min_eigval[k]=eigval;
					if (!ret) onlyerror1=false;

					eigvec=0; // using eigvec for mineig as approximation for maxeig gives trouble for some examples
					if (fkb->A[l]) ret=SumMatrix(H, fkb->A[l], -1., -1.).eig(eigvec, eigval, param);
					else ret=MinusMatrix(H).eig(eigvec, eigval, param);
					if (ret) out_out << "Convexify::check_convex(): Error " << ret << " computing eigval for nonquad. func. block " << k << ',' << l << endl;
					if ((!ret) && -eigval>max_eigval[k]) max_eigval[k]=-eigval;
					if (!ret) onlyerror2=false;
				}
				ss_checked=true;
				if (onlyerror1 || onlyerror2) {
					out_err << "Couldn't compute any eigenvalue! Aborting." << endl;
					exit(-1);
				}
			}
		}

		if (min_eigval[k]<-rtol && max_eigval[k]>rtol) f.set_curvature(k, Func::INDEFINITE);
		else if (min_eigval[k]<-rtol) f.set_curvature(k, Func::CONCAVE); // so max_eigval<rtol
		else if (max_eigval[k]>rtol) f.set_curvature(k, Func::CONVEX); // so min_eigval>-rtol
		else f.set_curvature(k, Func::LINEAR); // so min_eigval>-rtol and max_eigval<rtol
		if (!f.A[k]) f.s[k]->set_curvature(f.get_curvature(k));
	}
}

void Convexify::check_convexification(MinlpProblem& conv_prob, MinlpProblem& orig_prob, const ivector& ineq_index, int test_pts) {
	dvector x(orig_prob.dim());
	double eig_val;
	for (int c=0; c<=conv_prob.con.size(); c++) {
		int nr=0, nr2=0;
		double diff=0., diff2=0.;
		for (int i=0; i<test_pts; i++) {
			x.set_random(conv_prob.lower, conv_prob.upper);
			SepQcFunc& f(c ? *conv_prob.con[c-1] : *conv_prob.obj);
			double val=f.eval(x);
			if (val!=val) continue;
			for (int k=0; k<f.block.size(); k++) { // checking eigenvalue for each block
				if ((!f.A[k]) && (!f.s[k])) continue;
				dvector eig_vec(f.block[k].size());
				int ret=HessMatrix(SepQcFunc(f.A[k], NULL, f.s[k]), x(f.block[k])).eig(eig_vec, eig_val, param);
				if ((!ret) && (eig_val < -10*rtol)) {
					out_log << c << ": " << (c ? &*conv_prob.con_names[c-1] : "objective") << " block " << k << " eigval = " << eig_val << endl;
					f.print(*out_log_p, conv_prob.var_names);
				}
			}

			if (c<orig_prob.con.size()+1) { // checking, if feas(C)>feas(P)
				double val1=(c ? orig_prob.con[c-1] : orig_prob.obj)->eval(x);
				double val2=f.eval(x);
				if (val1<1E-4 && val2>1E-4) {
					diff+=fabs(val2-val1);
					nr++;
					f.print(*out_log_p, conv_prob.var_names);
				}
//				if (val<1E-4 && diff < -rtol) out_log << (c ? orig_prob.con_names[c-1] : "objective") << " diff = " << diff << endl;
				if (c && orig_prob.con_eq[c-1] && (!conv_prob.con_eq[c-1])) {
					val2=conv_prob.con[ineq_index[c-1]]->eval(x);
					if (fabs(val1)<1E-4 && fabs(val2)>1E-4) {
						diff2+=fabs(-val1-val2);
						nr2++;
//					if (diff < -rtol) out_log << conv_prob.con_names[ineq_index[c-1]] << ", " << orig_prob.con_names[c-1] << " diff = " << diff << endl;
					}
				}
			}
		}
		if (nr) out_log << (c ? &*conv_prob.con_names[c-1] : "objective") << ": " << 100*(double)nr/(double)test_pts << "% points cut by (C), |(C)-(P)| in this points: " << diff << endl;
		if (nr2) out_log << conv_prob.con_names[ineq_index[c-1]] << ": " << 100*(double)nr2/(double)test_pts << "% points cut by (C), |(C)-(P)| in this points: " << diff2 << endl;
	}
}

