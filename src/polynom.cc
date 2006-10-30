// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "polynom.h"
#include "opt.h"

void Monom::grad(UserVector<double>& grad, const UserVector<double>& x) const {
	double val=eval(x);
	int count;
	for (int i=0; i<dim(); i++) {
		count=indices.count(i);
		grad[i] = (count && val) ? count*val/x(i) : 0.;
	}
}

void Monom::HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
	y=0;
	MultiIndex ind;
	double z_val;
	for (int j=0; j<dim(); j++) {
		z_val=z(j);
		if (fabs(z_val)>rtol)
			for (int i=0; i<dim(); i++) {
				ind.clear();
				ind.insert(i); ind.insert(j);
				y[i]+=z_val*part_derivate(x, ind);
			}
	}

}

double Monom::part_derivate_rek(const UserVector<double>& x, MultiIndex& alpha, MultiIndex& beta) const {
	if (!beta.size()) {
		double val=1.;
		for (MultiIndex::iterator it(alpha.begin()); it!=alpha.end(); it++) val*=x(*it);
		return val;
	}

	int count=alpha.count(*beta.begin());
	if (!count) return 0.;
	alpha.erase(alpha.find(*beta.begin()));
	beta.erase(beta.begin());
	return count*part_derivate_rek(x, alpha, beta);
}

// ----------------------------------------- PolynomialUnderestimator2 -------------------------------------------

PolynomialUnderestimator2::PolynomialUnderestimator2(Pointer<Param> param_)
: param(param_), decomp(param),
	sampling0(param_, "Polynomial Underestimator K0"), sampling1(param_, "Polynomial Underestimator K1"), sampling2(param_, "Polynomial Underestimator K2"),
	sampling_vertices(param_, "Polynomial Underestimator K0"), sampling_minimizer(param_, "Polynomial Underestimator K0"),
	sampling_initial(-1), ss_size(3)
{ max_degree=param ? param->get_i("Polynomial Underestimator max polynom degree", 2) : 2;
	if (param) {
		if (param->get_i("Polynomial Underestimator K2 sample set initial", 1)) sampling_initial=2;
		else if (param->get_i("Polynomial Underestimator K1 sample set initial", 0)) sampling_initial=1;
		else if (param->get_i("Polynomial Underestimator K0 sample set initial", 0)) sampling_initial=0;
	}
}

void PolynomialUnderestimator2::polynomial_underestimator(MinlpProblem& prob, MINLPData& minlpdata, ivector& ineq_index, SparseVector<double>& obj_c_add, vector<SparseVector<double> >& con_c_add) {
	ineq_index.resize(prob.con.size());

	c_add1.resize(prob.block.size());
	c_add2.resize(prob.block.size());

	bool nonquadblocks=false;
	for (int k=prob.block.size()-1; (!nonquadblocks) && k>=0; --k) nonquadblocks=prob.obj->s[k];
	if (nonquadblocks) {
		prob.obj=new SepQcFunc(*prob.obj);
		out_log << "objective: ";
		polynomial_underestimator(prob.obj, false, prob.lower, prob.upper, prob.primal_point);
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

		Pointer<SepQcFunc> fm=polynomial_underestimator(prob.con[c], prob.con_eq[c], prob.lower, prob.upper, prob.primal_point);
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

Pointer<SepQcFunc> PolynomialUnderestimator2::polynomial_underestimator(Pointer<SepQcFunc> f, bool eq, dvector& lower, dvector& upper, dvector& primal) {
	c_add1=0.;
	Pointer<SepQcFunc> fm;
	if (eq) {
		fm=new SepQcFunc(*f, true);
		c_add2=0.;
	}

	for (int k=0; k<f->block.size(); ++k) {
		if (!f->s[k]) continue;
		if ((f->get_curvature(k)&Func::CONVEX) && ((!eq) || f->get_curvature(k)&Func::CONCAVE)) continue;
		out_log << 'C' << f->get_curvature(k) << ' ';

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

			dvector low(fk->block[l].size());
			dvector up(low.dim());
			int i0;
			for (int i=0; i<low.dim(); ++i) {
				i0=f->block[k][fk->block[l][i]];
				low[i]=lower(i0);
				up[i]=upper(i0);
			}

			out_log << 'B' << low.dim() << ' ';

			new_sampleset(low, up);
			check_for_nan(*fk->s[l]);
			add_point_to_sampleset(prim(fk->block[l]));
			bool min_added=add_minimizer_to_sample(fk->s[l], low, up);

			out_log << "SS" << ss_size[0] << ',' << ss_size[1] << ',' << ss_size[2] << ' ';

			new_multiindices(((const Func*)(Func*)fk->s[l])->get_sparsity(), fk->s[l]->dim());
			out_log << 'C' << fk->s[l]->get_curvature() << ' ';

			if (!(fk->s[l]->get_curvature()&Func::CONVEX)) {
				polynomial_underestimator(*A, *b, c, *fk->s[l], fk->block[l]);
				fk->s[l]=NULL;
			}
			if (eq && !(fkm->s[l]->get_curvature()&Func::CONVEX)) {
			  if (min_added) remove_last_point_from_sample();
			  add_minimizer_to_sample(fkm->s[l], low, up);
				polynomial_underestimator(*Am, *bm, cm, *fkm->s[l], fkm->block[l]);
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
			} else fm->s[k]=NULL;
			fm->set_sparsity(k, sim);
		}
		out_log << "; \t";
	}
	out_log << endl;

	return fm;
}

void PolynomialUnderestimator2::polynomial_underestimator(SparseMatrix2& A, SparseVector<double>& b, double& c, Func& f, ivector& indices) {
	int multiindices_size=multiindices.size();
	int connr=ss_size[0]+(maxdegree1_size-1)*ss_size[1]+(maxdegree2_size-maxdegree1_size)*ss_size[2]; // number of constraints

	int varnr=multiindices_size+ // a's
		ss_size[0]+ // u's for ss0
		2*MIN(ss_size[0],ss_size[1])*(maxdegree1_size-1)+ // u's and v's for ss1
		2*MIN(ss_size[0],ss_size[1])*(maxdegree2_size-maxdegree1_size); // u's and v's for ss2

	MipProblem lp(varnr, connr);

	// u's and v's are positive
	for (int i=multiindices_size; i<varnr; ++i)
		lp.setColBounds(i, 0, INFINITY);

	Pointer<UserVector<double> > obj(new dvector(varnr));

	dvector grad(f.dim());
	SparseVector<double> e(f.dim()); // for unit vector
	dvector prod(f.dim()); // for product of hessian in *it_ss with unit vector

	double rhs;
	dvector b1(multiindices_size+1);
	b1[multiindices_size]=1.; // u
	dvector b2(multiindices_size+2);
	b2[multiindices_size]=1.; b2[multiindices_size+1]=-1.; // u and v
	ivector indices1(b1.dim());
	ivector indices2(b2.dim());
	for (int i=0; i<multiindices_size; ++i) {
		indices1[i]=i;
		indices2[i]=i;
	}

	int uv_ind=multiindices_size; // index of u and v variables
	int con_ind=0; // constraint index
	for (int j=0; j<ss_size[0]; j++) { // for all sample points
		list<MultiIndex>::iterator it_mind_beta(multiindices.begin());
		if (it_mind_beta->size()==0) { // constraint for closeneth of function value
			if (j<ss_size[0]) { // still degree0-constraints
				rhs=f.eval(sample_set[j]);
				assert(rhs==rhs);

				int i=0; // number of multiindex alpha
				for (list<Monom>::iterator it_monom(monoms.begin()); it_monom!=monoms.end(); it_monom++, i++)
					b1[i]=it_monom->eval(sample_set[j]);
				assert(i==multiindices_size);

				indices1[multiindices_size]=uv_ind;
				(*obj)[uv_ind]=200.; // delta_0 for u
				++uv_ind;

				lp.setRow(con_ind++, b1, indices1, rhs, rhs);
			}
			it_mind_beta++;
		}

		if (j<ss_size[1]) { // still degree1-constraints
			f.grad(grad, sample_set[j]);
			while (it_mind_beta->size()==1) { // constraints for closeneth of gradient values
				rhs=grad[*it_mind_beta->begin()];
				assert(rhs==rhs);

				int i=0; // number of multiindex alpha
				for (list<Monom>::iterator it_monom(monoms.begin()); it_monom!=monoms.end(); it_monom++, i++)
					b2[i]=it_monom->part_derivate(sample_set[j], *it_mind_beta);

				indices2[multiindices_size]=uv_ind;
				(*obj)[uv_ind]=2.;
				++uv_ind;

				indices2[multiindices_size+1]=uv_ind;
				(*obj)[uv_ind]=2.;
				++uv_ind;

				lp.setRow(con_ind++, b2, indices2, rhs, rhs);

				it_mind_beta++; // next multiindex finished
			}
		}

		if (j<ss_size[2]) { // still degree2-constraints
			int first; // first index in paritial derivate (of size 2)
			int second; // second index in paritial derivate (of size 2)
			while (it_mind_beta!=multiindices.end() && it_mind_beta->size()==2) { // constraints for closeneth of hessian values
				first=*it_mind_beta->begin();
				e.SetElement(first, 1., false);
				f.HessMult(prod, sample_set[j], e);
				do {
					second=*(++it_mind_beta->begin());

					double rhs=prod[second];
					assert(rhs==rhs);

					int i=0; // number of multiindex alpha
					for (list<Monom>::iterator it_monom(monoms.begin()); it_monom!=monoms.end(); it_monom++, i++)
						b2[i]=it_monom->part_derivate(sample_set[j], *it_mind_beta);

					indices2[multiindices_size]=uv_ind;
					(*obj)[uv_ind]=00.2;
					++uv_ind;

					indices2[multiindices_size+1]=uv_ind;
					(*obj)[uv_ind]=00.2;
					++uv_ind;

					lp.setRow(con_ind++, b2, indices2, rhs, rhs);

					it_mind_beta++; // next multiindex finished
				} while (it_mind_beta!=multiindices.end() && *it_mind_beta->begin()==first && it_mind_beta->size()==2);

				e.DelElement(first);
			}
		}
	}

	lp.setObj(obj);
	lp.finish();
	Pointer<MIPSolver> solver=MIPSolver::get_solver(lp, param);
//	solver->set_maxiter(10*varnr);
	MIPSolver::SolutionStatus ret=solver->solve();

	if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) {
		out_err << "U2 returned " << ret << ". Aborting" << endl;
		exit(-1);
	}
	out_log << "U2:" << solver->get_optval() << ' ';
	if (ret==MIPSolver::FEASIBLE) out_log << ret << ' ';

	dvector coeff(multiindices_size);
	solver->get_primal(coeff);

	int i=0;
	list<MultiIndex>::iterator it_mind(multiindices.begin());
	list<Monom>::iterator it_monom(monoms.begin());
	if (it_mind->size()==0) {
		c+=coeff[i++];
		it_mind++;
		it_monom++;
	}

	while (it_mind->size()==1) {
		b[indices(*it_mind->begin())]+=.5*coeff[i++];
		it_mind++;
		it_monom++;
	}

	while (it_mind->size()==2) {
		int second=*(++it_mind->begin());
		if (*it_mind->begin()!=second) {
		  A.AddToElement(indices(*it_mind->begin()), indices(second), .5*coeff[i]);
		  A.AddToElement(indices(second), indices(*it_mind->begin()), .5*coeff[i]);
		} else {
		  A.AddToElement(indices(second), indices(second), coeff[i]);
		}
		it_mind++;
		it_monom++;
		i++;
	}

	assert(max_degree<=2);
}


void PolynomialUnderestimator2::new_multiindices(const SparsityInfo& si, int n) {
	multiindices.clear();

	multiindices.push_back(MultiIndex());

	if (max_degree>=1)
		for (VariableIterator it(si); it; ++it)
			multiindices.push_back(it());

	maxdegree1_size=multiindices.size(); // the number of multiindices with degree <= 1

	if (max_degree>=2) {
		for (VariableIterator it(si, false, true, true); it; ++it) // only nonlinear variables
			multiindices.push_back(MultiIndex(it(), it()));
		for (map<pair<int,int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it)
			multiindices.push_back(MultiIndex(it->first.first, it->first.second));
	}

	maxdegree2_size=multiindices.size(); // the number of multiindices with degree <= 2

	multiindices.sort();

	monoms.clear();
	for (list<MultiIndex>::iterator it(multiindices.begin()); it!=multiindices.end(); it++)
		monoms.push_back(Monom(n, *it));
}

void PolynomialUnderestimator2::new_sampleset(const dvector& lower, const dvector& upper) {
	sample_set.clear();
	ss_size=0;

	sampling2.get_points(sample_set, lower, upper); // K2
	ss_size[2]=sample_set.size();

	sampling1.get_points(sample_set, lower, upper); // K1
	ss_size[1]=sample_set.size();

	sampling0.get_points(sample_set, lower, upper); // K0
	sampling_vertices.get_points(sample_set, lower, upper);
	ss_size[0]=sample_set.size();
}

void PolynomialUnderestimator2::check_for_nan(const Func& f) {
	int i=0;
	vector<dvector>::iterator it(sample_set.begin());
	while (it!=sample_set.end()) {
  	double val=f.eval(*it);
		if (finite(val)) { ++it; ++i; continue; }
		out_log << "Removing point which value is nan from sample set" << endl;
		it=sample_set.erase(it);
		if (i<ss_size[0]) --ss_size[0];
		if (i<ss_size[1]) --ss_size[1];
		if (i<ss_size[2]) --ss_size[2];
	}
}

bool PolynomialUnderestimator2::add_point_to_sampleset(const dvector& point) {
	if (sampling_initial<0) return false;
	vector<dvector>::iterator pos;
	switch(sampling_initial) {
		case 2: pos=sample_set.begin(); // start of K2
			++ss_size[0]; ++ss_size[1]; ++ss_size[2];
			break;
		case 1: pos=sample_set.end(); // end of K1
			for (int i=sample_set.size()-1; i>=ss_size[1]; --i) --pos;
			++ss_size[0]; ++ss_size[1];
			break;
		case 0: pos=sample_set.begin(); // end of K0
			for (int i=0; i<ss_size[1]; ++i) ++pos;
			++ss_size[0];
			break;
		default: return false;
	}
	sample_set.insert(pos, point);
	return true;
}

bool PolynomialUnderestimator2::add_minimizer_to_sample(Pointer<Func> f, const dvector& lower, const dvector& upper) {
	if (sampling_minimizer.add_minimizer(sample_set, f, lower, upper)) {
		++ss_size[0];
		return true;
	}
	return false;
}

void PolynomialUnderestimator2::remove_last_point_from_sample() {
	if (!sample_set.size()) return;
	sample_set.pop_back();
	--ss_size[0];
}

void PolynomialUnderestimator2::check(MinlpProblem& prob, MinlpProblem& quad, ivector& ineq_index) {
	dvector x(prob.dim());
	int errors=0;
	for (int c=0; c<prob.con.size(); ++c) {
		int err=0;
		for (int i=0; i<20; ++i) {
			x.set_random(prob.lower, prob.upper);
			double val1=prob.con[c]->eval(x);
			double val2=quad.con[c]->eval(x);
			double reldiff=(val1-val2)/(fabs(val1)+1);
			if (reldiff<-.01) {
				++err;
				out_log << "PolynomialUnderestimator check: " << prob.con_names[c] << " not underestimating. rel. diff = " << reldiff << "\t " << val1 << "\t " << val2 << endl;
			}

			if (ineq_index[c]) {
				val1*=-1;
				val2=quad.con[ineq_index[c]]->eval(x);
				reldiff=(val1-val2)/(fabs(val1)+1);
				if (reldiff<-.01) {
					++err;
					out_log << "PolynomialUnderestimator check: -" << prob.con_names[c] << " not underestimating. rel. diff = " << reldiff << "\t " << val1 << "\t " << val2 << endl;
				}
			}
		}
		if (err>5) {
			prob.con[c]->print(*out_log_p, prob.var_names);
			quad.con[c]->print(*out_log_p, prob.var_names);
			if (ineq_index[c]) quad.con[ineq_index[c]]->print(*out_log_p, prob.var_names);
			++errors;
		}
	}
	out_log << "PolynomialUnderestimator check: " << errors << " errors" << endl;
}
