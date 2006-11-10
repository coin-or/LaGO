// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "quaduest.h"
#include "opt.h"


// ----------------------------------------- QuadraticUnderestimator -------------------------------------------

QuadraticUnderestimator::QuadraticUnderestimator(Pointer<Param> param_)
: param(param_), decomp(param),
	eps(param_->get_d("Quadratic Underestimator epsilon", 1E-3)),
	iter_max(param_->get_i("Quadratic Underestimator iteration limit", 100)),
	sampling(param_, "Quadratic Underestimator"),
	sampling_vertices(param_, "Quadratic Underestimator"), sampling_minimizer(param_, "Quadratic Underestimator"),
	sampling_initial(param_->get_i("Quadratic Underestimator sample set initial", 1))
{ }

void QuadraticUnderestimator::new_multiindices(const SparsityInfo& si, int n) {
	multiindices.clear();

	multiindices.push_back(MultiIndex());

	for (VariableIterator it(si); it; ++it)
		multiindices.push_back(it());

	for (VariableIterator it(si, false, true, true); it; ++it) // only nonlinear variables
		multiindices.push_back(MultiIndex(it(), it()));
	for (map<pair<int,int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it)
		multiindices.push_back(MultiIndex(it->first.first, it->first.second));

	multiindices.sort();

	monoms.clear();
	for (list<MultiIndex>::iterator it(multiindices.begin()); it!=multiindices.end(); it++)
		monoms.push_back(Monom(n, *it));
}

void QuadraticUnderestimator::new_sampleset(const dvector& lower, const dvector& upper) {
	sample_set.clear();
	
	vector<dvector> sample_set_;

	sampling.get_points(sample_set_, lower, upper);
	sampling_vertices.get_points(sample_set_, lower, upper);
	
	for (int i=0; i<sample_set_.size(); ++i)
		add_point_to_sampleset(sample_set_[i]);
}

void QuadraticUnderestimator::check_for_nan(const Func& f) {
	multimap<double, dvector>::iterator it(sample_set.begin());
	while (it!=sample_set.end()) {
  	double val=f.eval(it->second);
		if (finite(val)) { ++it; continue; }
		out_log << "Value nan found for sample point " << it->second;
		multimap<double, dvector>::iterator next(it); ++next;
		sample_set.erase(it);
		it=next;
	}
}

multimap<double, dvector>::iterator QuadraticUnderestimator::add_point_to_sampleset(const dvector& point) {
	double norm=point.sq_norm2();
	
	// check whether point already in sample set
	multimap<double, dvector>::iterator it(sample_set.lower_bound(norm));
	while (it!=sample_set.end() && it->first>norm-1E-4 && it->first<norm+1E-4) {
		bool equal=true;
		for (int i=0; i<point.dim(); ++i)
			if (fabs(it->second(i)-point(i))>rtol) equal=false;
		if (equal) return it;
		++it;
	}

	return sample_set.insert(it, pair<double, dvector>(norm, point));
}

multimap<double, dvector>::iterator QuadraticUnderestimator::add_minimizer_to_sample(Pointer<Func> f, const dvector& lower, const dvector& upper) {
	vector<dvector> sample_set_;
	bool ok=sampling_minimizer.add_minimizer(sample_set_, f, lower, upper);
	if (ok) return add_point_to_sampleset(sample_set_.front());
	return sample_set.end();
}

void QuadraticUnderestimator::quadratic_underestimator(MinlpProblem& prob, MINLPData& minlpdata, ivector& ineq_index, SparseVector<double>& obj_c_add, vector<SparseVector<double> >& con_c_add) {
	ineq_index.resize(prob.con.size());

	c_add1.resize(prob.block.size());
	c_add2.resize(prob.block.size());

	bool nonquadblocks=false;
	for (int k=prob.block.size()-1; (!nonquadblocks) && k>=0; --k) nonquadblocks=prob.obj->s[k];
	if (nonquadblocks) {
		prob.obj=new SepQcFunc(*prob.obj);
		out_log << "objective: ";
		quadratic_underestimator(prob.obj, false, prob.lower, prob.upper, prob.primal_point);
		if (obj_c_add.size()) obj_c_add=c_add1;
		minlpdata.obj.polynomial_underestimator=prob.obj; //TODO: this should only be set if the nonquadratic parts of the objective are not convex
		for (int k=0; k<c_add1.size(); ++k)
			if (c_add1(k))
				minlpdata.obj.polynomial_approx_constants_lower.insert(pair<int, double>(k, c_add1(k)));
	}

	int consize=prob.con.size();
	for (int c=0; c<consize; c++) {
		nonquadblocks=false;
		for (int k=prob.block.size()-1; (!nonquadblocks) && k>=0; --k) nonquadblocks=prob.con[c]->s[k];
		if (!nonquadblocks) continue; // linear or quadratic function

		prob.con[c]=new SepQcFunc(*prob.con[c]);
		out_log << "con " << c << ": "; if (prob.con_names[c]) out_log << prob.con_names[c] << ": ";

		Pointer<SepQcFunc> fm=quadratic_underestimator(prob.con[c], prob.con_eq[c], prob.lower, prob.upper, prob.primal_point);
		if (con_c_add.size()) con_c_add[c]=c_add1;
		minlpdata.con[c].polynomial_underestimator=prob.con[c]; //TODO: this should only be set if the nonquadratic parts of the constraints are not convex
		for (int k=0; k<c_add1.size(); ++k)
			if (c_add1(k))
				minlpdata.con[c].polynomial_approx_constants_lower.insert(pair<int, double>(k, c_add1(k)));

		if (prob.con_eq[c]) {
			ineq_index[c]=prob.con.size();
			char* name=NULL;
			if (prob.con_names[c]) sprintf(name=new char[strlen(prob.con_names[c])+2], "-%s", (char*)prob.con_names[c]);
			prob.add_con(fm, false, name);
			if (name) delete name;
			prob.con_eq[c]=false;
			if (con_c_add.size()) con_c_add.push_back(c_add2);
			minlpdata.con[c].polynomial_overestimator=fm; //TODO: this should only be set if the nonquadratic parts of the objective are not concave
			for (int k=0; k<c_add2.size(); ++k)
				if (c_add2(k))
					minlpdata.con[c].polynomial_approx_constants_upper.insert(pair<int, double>(k, c_add2(k)));
		}
	}
}

Pointer<SepQcFunc> QuadraticUnderestimator::quadratic_underestimator(Pointer<SepQcFunc> f, bool eq, dvector& lower, dvector& upper, dvector& primal) {
	c_add1=0.;
	Pointer<SepQcFunc> fm;
	if (eq) {
		fm=new SepQcFunc(*f, true);
		c_add2=0.;
	}

	for (int k=0; k<f->block.size(); ++k) {
		if (!f->s[k]) continue;
		if ((f->get_curvature(k)&Func::CONVEX) && ((!eq) || f->get_curvature(k)&Func::CONCAVE)) continue;
// 		out_log << 'C' << f->get_curvature(k) << ' ';

		dvector prim(primal, f->block[k]);

		Pointer<SepQcFunc> fk=decomp.decompose(f->s[k], prim);
		out_log << 'D' << fk->block.size() << ' ';

		Pointer<SepQcFunc> fkm;
		if (eq) fkm=new SepQcFunc(*fk, true); // minus f_k

		Pointer<SparseMatrix2> A(new SparseMatrix2(fk->dim()));
		Pointer<SparseVector<double> > b(new SparseVector<double>(fk->b, fk->block));
		double c=fk->c; fk->c=0.;
		Pointer<SparsityInfo> si(new SparsityInfo(2));
		for (int l=0; l<fk->block.size(); ++l) {
			fk->b[l]=NULL;
			if (!fk->A[l]) continue;
			SparseMatrix2* Al=dynamic_cast<SparseMatrix2*>((UserMatrix*)fk->A[l]); assert(Al);
			A->set_block(*Al, fk->block[l]);
			fk->A[l]=NULL;
		}

		Pointer<SparseMatrix2> Am;
		Pointer<SparseVector<double> > bm;
		double cm;
		Pointer<SparsityInfo> sim;
		if (eq) {
			cm=fkm->c; fkm->c=0.;
			Am=new SparseMatrix2(fk->dim());
			bm=new SparseVector<double>(fkm->b, fkm->block);
			for (int l=0; l<fkm->block.size(); ++l) {
				fkm->b[l]=NULL;
				if (!fkm->A[l]) continue;
				SparseMatrix2* Al=dynamic_cast<SparseMatrix2*>((UserMatrix*)fkm->A[l]); assert(Al);
				Am->set_block(*Al, fkm->block[l]);
				fkm->A[l]=NULL;
			}
			sim=new SparsityInfo(2);
		}

		for (int l=0; l<fk->block.size(); ++l) {
			if (!fk->s[l]) continue;
			if ((fk->s[l]->get_curvature()&Func::CONVEX) && (!eq || (fkm->s[l]->get_curvature()&Func::CONVEX))) continue;

			Pointer<dvector> low=new dvector(fk->block[l].size());
			Pointer<dvector> up=new dvector(low->dim());
			int i0;
			for (int i=0; i<low->dim(); ++i) {
				i0=f->block[k][fk->block[l][i]];
				(*low)[i]=lower(i0);
				(*up)[i]=upper(i0);
			}

			out_log << 'B' << low->dim() << ' ';

			new_sampleset(*low, *up);
			check_for_nan(*fk->s[l]);
			multimap<double, dvector>::iterator initial=sample_set.end();
			if (sampling_initial) initial=add_point_to_sampleset(prim(fk->block[l]));
			multimap<double, dvector>::iterator minimizer=add_minimizer_to_sample(fk->s[l], *low, *up);
			
			if (minimizer!=sample_set.end()) enforce_tightness=minimizer;
			else if (initial!=sample_set.end()) enforce_tightness=initial;
			else if (!sample_set.empty()) enforce_tightness=--sample_set.end(); // set to last element in sample set
			else { out_err << "Empty sampleset. Aborting." << endl; exit(-1); }

			out_log << "SS" << sample_set.size();

			new_multiindices(((const Func*)(Func*)fk->s[l])->get_sparsity(), fk->s[l]->dim());
			out_log << 'C' << fk->s[l]->get_curvature() << ' ';

			if (!(fk->s[l]->get_curvature()&Func::CONVEX)) {
				quadratic_underestimator(*A, *b, c, fk->s[l], fk->block[l], low, up);
				fk->s[l]=NULL;
			}
			if (eq && !(fkm->s[l]->get_curvature()&Func::CONVEX)) {
			  if (minimizer!=sample_set.end()) {
					sample_set.erase(minimizer);
					minimizer=sample_set.end();
				}
				minimizer=add_minimizer_to_sample(fkm->s[l], *low, *up);
			
				if (minimizer!=sample_set.end()) enforce_tightness=minimizer;
				else if (initial!=sample_set.end()) enforce_tightness=initial;
				else if (!sample_set.empty()) enforce_tightness=--sample_set.end();
				else { out_err << "Empty sampleset. Aborting." << endl; exit(-1); }
				
 				quadratic_underestimator(*Am, *bm, cm, fkm->s[l], fkm->block[l], low, up);
				fkm->s[l]=NULL;
			}
		}

		if (c_add1.size()) c_add1[k]=c;
		if (eq && c_add2.size()) c_add2[k]=cm;

		// merge fk into f
		if (f->A[k]) si->add(*f->A[k]);
		if (A->nonzeros()) {
			A->finish();
			si->add(*A);
			if (f->A[k]) f->A[k]=new SumMatrix(f->A[k], (Pointer<UserMatrix>)A);
			else f->A[k]=A;
		}
		if (f->b[k]) si->add(*f->b[k]);
		if (!(*b==0)) {
			si->add(*b);
			if (f->b[k]) {  f->b[k]=f->b[k]->getcopy(); *f->b[k]+=*b; }
			else f->b[k]=b;
		}
		f->c+=c;
		bool s_left=false;
		for (int l=0; l<fk->block.size(); ++l)
			if (fk->s[l]) s_left=true;
		if (s_left) {
			fk->set_sparsity();
			f->s[k]=fk;
			si->add(fk->get_sparsity());
		} else f->s[k]=NULL;
		f->set_sparsity(k, si);

		// merge fkm into fm (=ret.second.first)
		if (eq) {
			if (fm->A[k]) sim->add(*fm->A[k]);
			if (Am->nonzeros()) {
				Am->finish();
				sim->add(*Am);
				if (fm->A[k]) fm->A[k]=new SumMatrix(fm->A[k], (Pointer<UserMatrix>)Am);
				else fm->A[k]=Am;
			}
			if (fm->b[k]) sim->add(*fm->b[k]);
			if (!(*bm==0)) {
				sim->add(*bm);
				if (fm->b[k]) *fm->b[k]+=*bm;
				else fm->b[k]=bm;
			}
			fm->c+=cm;
			bool s_left=false;
			for (int l=0; l<fkm->block.size(); ++l)
				if (fkm->s[l]) s_left=true;
			if (s_left) {
				fkm->set_sparsity();
				fm->s[k]=fkm;
				sim->add(fkm->get_sparsity());
			} else fm->s[k]=NULL;
			fm->set_sparsity(k, sim);
		}
		out_log << "; \t";
	}
	out_log << endl;

	return fm;
}

void QuadraticUnderestimator::quadratic_underestimator(SparseMatrix2& A, SparseVector<double>& b, double& c, const Pointer<Func>& f, ivector& indices, const Pointer<dvector>& lower, const Pointer<dvector>& upper) {
// create auxiliary LP with initial sample set
	int multiindices_size=multiindices.size(); // sparsity of p(x)
	int aux_vars=3;
	MipProblem lp(multiindices_size+aux_vars, sample_set.size());

	Pointer<UserVector<double> > obj=new dvector(multiindices_size+aux_vars);
	dvector b1(multiindices_size+aux_vars);
	list<double> rhs;
	double obj_const=0;

	int j=0;
	for (multimap<double, dvector>::iterator it_sample_point(sample_set.begin()); it_sample_point!=sample_set.end(); ++it_sample_point, ++j) {
		double f_val=f->eval(it_sample_point->second);
		assert(f_val==f_val);

		double scale=1;
		int i=0; // number of multiindex alpha
		for (list<Monom>::iterator it_monom(monoms.begin()); it_monom!=monoms.end(); it_monom++, i++) {
			b1[i]=it_monom->eval(it_sample_point->second);
			if (fabs(b1[i])>scale) scale=fabs(b1[i]);
		}
		b1/=scale;
		rhs.push_back(f_val/scale);

/*		if (j==0 || j==sample_set.size()-1 || it_sample_point==enforce_tightness) {
			*obj-=b1;
			obj_const+=rhs.back();
		}*/
		double lhs=-INFINITY;
		if (j==0) {
			b1[multiindices_size]=1.;
			(*obj)[multiindices_size]=1.;
			lhs=rhs.back();
			lp.setColBounds(multiindices_size, 0, INFINITY);
		} else if (j==sample_set.size()-1) {
			b1[multiindices_size+1]=1.;
			(*obj)[multiindices_size+1]=1.;
			lhs=rhs.back();
			lp.setColBounds(multiindices_size+1, 0, INFINITY);
		} else if (it_sample_point==enforce_tightness) {
			b1[multiindices_size+2]=1.;
			(*obj)[multiindices_size+2]=1.;
			lhs=rhs.back();
			lp.setColBounds(multiindices_size+2, 0, eps*fabs(rhs.back()));
		}
		lp.setRow(j, b1, lhs, rhs.back());
		
		b1=0.;
		
/*		if (it_sample_point==enforce_tightness) lp.setRow(j, b1, rhs.back()-eps*fabs(rhs.back()), rhs.back());
		else lp.setRow(j, b1, -INFINITY, rhs.back());*/
	}

	lp.setObj(obj, obj_const);
	lp.finish();
	Pointer<MIPSolver> solver=MIPSolver::get_solver(lp, param);
// 	solver->set_tol(1E-5);
//	solver->set_maxiter(10*varnr);
	
	dvector best_coeff(multiindices_size);
	double best_violation=-INFINITY;
	double best_U3_val;
	
	bool finished;
	int iter=0;
	do {
		MIPSolver::SolutionStatus ret=solver->solve();
	
		if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) {
			out_err << "U3 returned " << ret << '.';
			if (best_violation>-INFINITY) {
				out_log << " Keeping prior found solution with violation " << best_violation << ". ";
				break;
			} else {
				out_err	<< " Aborting." << endl;
				exit(-1);
			}
		}
		double U3_val=solver->get_optval();
// 		if (ret==MIPSolver::FEASIBLE) out_log << ret << ' ';
	
		dvector coeff(multiindices_size);
		solver->get_primal(coeff);
		
		dvector rowact(solver->nr_row());
		solver->get_rowactivity(rowact);
		
		// setup f(x)-p(x) as SepQcFunc
		SepQcFunc fpdiff(f->dim());
		fpdiff.s[0]=f;
		fpdiff.b[0]=new dvector(f->dim());
		Pointer<SparseMatrix2> A1=new SparseMatrix2(f->dim());
		int i=0;
		for (list<MultiIndex>::iterator it_mind(multiindices.begin()); it_mind!=multiindices.end(); ++it_mind, ++i) {
			switch(it_mind->size()) {
				case 0: fpdiff.c-=coeff[i]; break;
				case 1:	(*fpdiff.b[0])[*it_mind->begin()]-=.5*coeff[i]; break;
				case 2: int second=*(++it_mind->begin());
					if (*it_mind->begin()!=second) {
						A1->AddToElement(*it_mind->begin(), second, -.5*coeff[i]);
						A1->AddToElement(second, *it_mind->begin(), -.5*coeff[i]);
					} else {
						A1->AddToElement(second, second, -coeff[i]);
					}
			}
		}
		A1->finish();
		fpdiff.A[0]=A1;
//  		out_log << fpdiff;
		BoxLocOpt locminsolver(fpdiff, lower, upper);
	
		finished=true;
		double maxviol=0;
		double maxviol_unscaled=0;
		// check active constraints to determine ''active'' sample points
		// and start local minimization of f-p from active sample points
		multimap<double, dvector>::iterator it_sample_point(sample_set.begin()); 
		list<double>::iterator it_rhs(rhs.begin());
		for (int j=0; j<rowact.size(); ++j, ++it_sample_point, ++it_rhs) {
			if (rowact(j)-*it_rhs<-1E-4) continue;
// 			out_log << "U3 constraint active " << rowact(j)-*it_rhs << " for sample point " << it_sample_point->second;
			if (!finished) continue;
			
			int ret=locminsolver.solve(it_sample_point->second);
			if (ret) out_log << 'm' << ret;
			double optval=locminsolver.opt_val();
			if (!(optval==optval) || optval>=INFINITY) continue; // locmin give no useful result
			
			double f_val=f->eval(locminsolver.sol_point);
			double viol=optval/max(1.,fabs(f_val));
			if (maxviol>viol) { maxviol=viol; maxviol_unscaled=optval; }
			
			double viol2=0.;
			double scale2=1.;
			if (viol<-1E-4) { // compute violation according to second (LP) scaling
				int i=0; // number of multiindex alpha
				for (list<Monom>::iterator it_monom(monoms.begin()); it_monom!=monoms.end(); it_monom++, i++) {
					b1[i]=it_monom->eval(locminsolver.sol_point);
					if (fabs(b1[i])>scale2) scale2=fabs(b1[i]);
				}
				viol2=optval/scale2;
			}

// 			out_log << "\tLocMin from sample point returns " << ret << " and value " << optval << " scaled " << viol << ' ' << viol2 << " in point " << locminsolver.sol_point;

			if (viol<-1E-4 && viol2<-1E-4) { // too large, add point to sample set and LP and restart
				multimap<double, dvector>::iterator newpoint(add_point_to_sampleset(locminsolver.sol_point));
				if (newpoint!=sample_set.end()) {
					b1/=scale2;
					double buffer=eps*MIN(1, sqrt(locminsolver.sol_point.dist(it_sample_point->second)/MAX(1, it_sample_point->first)))*fabs(f_val/scale2);
					rhs.push_back(f_val/scale2-buffer);
					solver->add_row(b1, -INFINITY, rhs.back());
					finished=false;
// 	 				out_log << "\tadded row with rhs " << rhs.back() << " and coeff " << b1;
				} else {
					out_log << "Local minimizer already in sample set, not adding again." << endl;
				}
			}
		}
		out_log << 'v' << -maxviol << ' ';

		if (maxviol_unscaled>best_violation || finished) {
			best_violation=maxviol_unscaled;
			best_coeff=coeff;
			best_U3_val=U3_val;
		}

		if (iter>=iter_max) {
			out_log << "Iteration limit (" << iter_max << ") of quadratic underest. improvement exceeded. Keeping prior found solution with violation " << best_violation << ' ';
			finished=true;
		}
		
		++iter;
	} while (!finished);

	c+=best_violation; // lower underestimator by maximum known violation
	out_log << "U3:" << best_U3_val << ' ';

	// set A, b, and c.
	int i=0;
	list<MultiIndex>::iterator it_mind(multiindices.begin());
	if (it_mind!=multiindices.end() && it_mind->size()==0) {
		c+=best_coeff[i++];
		it_mind++;
	}

	while (it_mind!=multiindices.end() && it_mind->size()==1) {
		b[indices(*it_mind->begin())]+=.5*best_coeff[i++];
		it_mind++;
	}

	while (it_mind!=multiindices.end()) {
		int second=*(++it_mind->begin());
		if (*it_mind->begin()!=second) {
		  A.AddToElement(indices(*it_mind->begin()), indices(second), .5*best_coeff[i]);
		  A.AddToElement(indices(second), indices(*it_mind->begin()), .5*best_coeff[i]);
		} else {
		  A.AddToElement(indices(second), indices(second), best_coeff[i]);
		}
		it_mind++;
		i++;
	}

}
