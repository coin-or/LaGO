// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "minlpopt.h"
#include "sampling.h"
#include "bcp.h"
#include "relaxopt.h"
#include "linrelax.h"
#include "cuts.h"
#include "polynom.h"
#include "quaduest.h"

// --------------------------------- Reformulation -----------------------------

void Reformulation::add_var(vector<vector<Pointer<BlockMatrix> > >& A, vector<vector<Pointer<SepQcFunc> > >& s, MinlpProblem& prob, int k, int index, double lower, double upper, Pointer<char> name) {
	for (int c=0; c<=prob.con.size(); c++) {
		Pointer<SepQcFunc> f(c ? prob.con[c-1] : prob.obj);
		if (f->A[k] && (!A[c][k])) { // hide matrix in BlockMatrix
			A[c][k]=new BlockMatrix(f->A[k]);
			A[c][k]->block.push_back(ivector(0));
			A[c][k]->A.push_back(NULL);
		}
		f->A[k]=NULL;

		if (f->s[k] && (!s[c][k])) // hide function in SepQcFunc
			s[c][k]=new SepQcFunc(NULL, NULL, f->s[k]);
		f->s[k]=NULL;

		if (A[c][k]) {
			A[c][k]->block[1].resize(A[c][k]->block[1].size()+1); // increase size of block of overlapping variables
			A[c][k]->block[1][A[c][k]->block[1].size()-1]=prob.block[k].size();
			A[c][k]->set_dim();
		}

		if (s.size() && s[c][k])
			s[c][k]->add_var(s[c][k]->dim(), 1); // adding variable to empty block
	}

	prob.add_var(index, k, false, lower, upper, name);

	for (int c=0; c<=prob.con.size(); c++) {
		if (A[c][k]) (c ? prob.con[c-1] : prob.obj)->A[k]=A[c][k];
		if (s.size() && s[c][k]) (c ? prob.con[c-1] : prob.obj)->s[k]=s[c][k];
	}
}

void Reformulation::add_con(vector<vector<Pointer<BlockMatrix> > >& A, vector<vector<Pointer<SepQcFunc> > >& s, MinlpProblem& prob, int c, int k) {
	Pointer<SepQcFunc> con(new SepQcFunc(prob.block)); // new constraint for this block
	for (int l=0; l<prob.block.size(); ++l) {
		con->set_curvature(l, Func::LINEAR);
		if (k!=l) con->set_sparsity(l, new SparsityInfo(2));
	}

	SepQcFunc& f(c ? *prob.con[c-1] : *prob.obj);

	con->A[k]=f.A[k];
	if (f.b[k]) con->b[k]=new SparseVector<double>(*f.b[k]);
	con->s[k]=f.s[k];
	con->set_curvature(k, f.get_curvature(k));
	Pointer<SparsityInfo> si(new SparsityInfo(f.get_sparsity(k)));

	A.push_back(vector<Pointer<BlockMatrix> >(f.block.size()));
	A.back()[k]=A[c][k];
	A[c][k]=NULL;
	f.A[k]=NULL;

	if (s.size()) {
		s.push_back(vector<Pointer<SepQcFunc> >(f.block.size()));
		s.back()[k]=s[c][k];
		s[c][k]=NULL;
		f.s[k]=NULL;
	}

	f.b[k]=new SparseVector<double>(f.block[k].size(), f.block[k].size()-1, .5);
	if (!con->b[k]) con->b[k]=new SparseVector<double>(f.block[k].size(), f.block[k].size()-1, -.5);
	else con->b[k]->SetElement(f.block[k].size()-1, -.5);
	si->linear->insert(pair<int, SparsityInfo::LinearVariable>(f.block[k].size()-1, SparsityInfo::LinearVariable(-1.)));

	Pointer<SparsityInfo> si_f(new SparsityInfo(2));
	si_f->linear->insert(pair<int, SparsityInfo::LinearVariable>(f.block[k].size()-1, SparsityInfo::LinearVariable(1.)));

	f.set_curvature(k, Func::LINEAR);
	f.set_sparsity(k, si_f);
	con->set_sparsity(k, si);

	Pointer<char> name=new char[5+(c ? strlen((char*)prob.con_names[c-1]) : 10)];
	if (c) sprintf((char*)name, "%s_%d", (char*)prob.con_names[c-1], k);
	else sprintf((char*)name, "objective_%d", k);
	prob.add_con(con, false, name);

	int c2=c ? opt.ineq_index[c-1] : 0;
	if (c2) {
		prob.con_eq[c-1]=true;
		con=new SepQcFunc(prob.block); // new constraint for this block
		for (int l=0; l<prob.block.size(); ++l) {
			con->set_curvature(l, Func::LINEAR);
			if (k!=l) con->set_sparsity(l, new SparsityInfo(2));
		}

		con->A[k]=prob.con[c2]->A[k];
		if (prob.con[c2]->b[k]) con->b[k]=new SparseVector<double>(*prob.con[c2]->b[k]);
		con->s[k]=prob.con[c2]->s[k];
		con->set_curvature(k, prob.con[c2]->get_curvature(k));
		si=new SparsityInfo(((const SepQcFunc*)(SepQcFunc*)prob.con[c2])->get_sparsity(k));

		A.push_back(vector<Pointer<BlockMatrix> >(prob.block.size()));
		A.back()[k]=A[c2+1][k];
		A[c2+1][k]=NULL;

		if (s.size()) {
			s.push_back(vector<Pointer<SepQcFunc> >(prob.block.size()));
			s.back()[k]=s[c2+1][k];
			s[c2+1][k]=NULL;
		}

		if (!con->b[k]) con->b[k]=new SparseVector<double>(prob.block[k].size(), prob.block[k].size()-1, .5);
		else con->b[k]->SetElement(prob.block[k].size()-1, .5);
		si->linear->insert(pair<int, SparsityInfo::LinearVariable>(prob.block[k].size()-1, SparsityInfo::LinearVariable(1.)));

//		prob.con[c2]->set_curvature(k, Func::LINEAR);
		con->set_sparsity(k, si);

		Pointer<char> name=new char[5+strlen((char*)prob.con_names[c2])];
		sprintf((char*)name, "%s_%d", (char*)prob.con_names[c2], k);
		prob.add_con(con, false, name);
	}
}

void Reformulation::reformulate() {
	out_log << "Reformulating... " << endl;
	Timer t;
//	bool primal_feasible=opt.split_prob->feasible(opt.split_prob->primal_point, opt.tol, NULL)==0;
bool primal_feasible=false;

#ifdef FILIB_AVAILABLE
	vector<IntervalVector> box(opt.split_prob->block.size());
	for (int k=0; k<box.size(); k++)
		box[k]=IntervalVector(opt.split_prob->lower(opt.split_prob->block[k]), opt.split_prob->upper(opt.split_prob->block[k]));
#else
	vector<Pointer<dvector> > low(opt.split_prob->block.size());
	vector<Pointer<dvector> > up(opt.split_prob->block.size());
	for (int k=0; k<low.size(); k++) {
		low[k]=new dvector(opt.split_prob->lower, opt.split_prob->block[k]);
		up[k]=new dvector(opt.split_prob->upper, opt.split_prob->block[k]);
	}
#endif

	// upper bounds of objective
	dvector obj_upper(opt.split_prob->block.size());
	if (primal_feasible) {
		for (int k=0; k<obj_upper.dim(); k++) {
			if (opt.split_prob->obj->A[k] || opt.split_prob->obj->s[k])
				obj_upper[k]=opt.split_prob->obj->eval(opt.split_prob->primal_point(opt.split_prob->block[k]), k);
		}
	} else {
#ifdef FILIB_AVAILABLE
		for (int k=0; k<obj_upper.dim(); k++)
			if (opt.split_prob->obj->A[k] || opt.split_prob->obj->s[k])
				obj_upper[k]=opt.split_prob->obj->eval(box[k], k).sup();
#else
		out_out << "Computing upper bound of objective" << endl;
		Pointer<SepQcFunc> mobj(new SepQcFunc(*opt.split_prob->obj, true)); // minus objective
		dvector conv_mobj_c_add(mobj->block.size());
		if (!(opt.split_prob->obj->get_curvature()&Func::CONCAVE)) { // convexifying it
			Convexify convexify(opt.param);
			convexify.new_sampleset(opt.split_prob->lower, opt.split_prob->upper, *mobj);
			convexify.convexify(*mobj, false, opt.max_eigval.front(), opt.min_eigval.front(), opt.split_prob->lower, opt.split_prob->upper);
			pair<Pointer<SepQcFunc>, Pointer<SepQcFunc> > conv(mobj, NULL);
			convexify.get_decomposed_functions(conv);
			mobj=conv.first;
			for (int k=0; k<mobj->block.size(); k++)
				conv_mobj_c_add[k]=convexify.convexify_c[k].first;
		}
		for (int k=0; k<obj_upper.dim(); k++)
			if (opt.split_prob->obj->A[k] || opt.split_prob->obj->s[k]) {
				SepQcFunc f(mobj->A[k], mobj->b[k], mobj->s[k], conv_mobj_c_add(k));
				BoxLocOpt locopt(f, new dvector(opt.split_prob->lower, opt.split_prob->block[k]), new dvector(opt.split_prob->upper, opt.split_prob->block[k]));
				dvector start(opt.split_prob->primal_point, opt.split_prob->block[k]);
				int ret=locopt.solve(start);
				if (ret) out_log << "block " << k << " LocOpt return " << ret << " \t new bound: " << -locopt.opt_val() << endl;
				obj_upper[k]=-locopt.opt_val();
			}
#endif
	}

	// lower and upper bounds for t's
	vector<SparseVector<double> > lower(opt.orig_prob->con.size()+1);
	vector<SparseVector<double> > upper(opt.orig_prob->con.size()+1);
	for (int c=0; c<lower.size(); c++) {
		SepQcFunc& f(c ? *opt.split_prob->con[c-1] : *opt.split_prob->obj);
		lower[c].resize(f.block.size());
		upper[c].resize(f.block.size());
		for (int k=0; k<f.block.size(); k++) {
			// skip linear blocks
			if (!f.A[k] && !f.s[k]) continue;

#ifdef FILIB_AVAILABLE
			interval<double> tbox(f.eval(box[k], k));
			lower[c].SetElement(k, tbox.inf());
			upper[c].SetElement(k, tbox.sup());
#else
			SepQcFunc& fc(c ? *opt.convex_prob->con[c-1] : *opt.convex_prob->obj);
			if ((!fc.A[k]) && (!fc.b[k]) && (!fc.s[k])) {
				out_log << "Warning: Convex underestimator of block " << k << " of " << (c ? opt.split_prob->con_names[c-1] : " objective ") << " is a constant." << endl;
				lower[c].SetElement(k, 0);
			} else {
				SepQcFunc fb(fc.A[k], fc.b[k], fc.s[k]);
				BoxLocOpt locopt(fb, low[k], up[k]);
				dvector start(opt.convex_prob->primal_point, fc.block[k]);
				int ret=locopt.solve(start);
				if (ret) {
					out_log << "Computed lower bound of block " << k;
					if (c) { out_log << " con " << opt.convex_prob->con_names[c-1]; }
					else out_log << " objective";
					out_log << "\t LocOpt return: " << ret << "\t value: " << locopt.opt_val() << endl;
//					if (c && out_log_p) opt.convex_prob->con[c-1]->print(*out_log_p, opt.convex_prob->var_names);
				}
				lower[c].SetElement(k, locopt.opt_val());
				upper[c].SetElement(k, INFINITY);
				if (!c) assert(locopt.opt_val()<=obj_upper[k]);
			}
#endif
		}
	}

	related_t.resize(opt.split_prob->con.size(), -1);

	// initialize extended problems
	ext_prob=new MinlpProblem(*opt.split_prob);
	ext_prob->obj=new SepQcFunc(*ext_prob->obj);
	for (int c=0; c<ext_prob->con.size(); c++) ext_prob->con[c]=new SepQcFunc(*ext_prob->con[c]);
	if (opt.quad_prob) {
		ext_quad_prob=new MinlpProblem(*opt.quad_prob);
		for (int c=0; c<ext_quad_prob->con.size(); c++) ext_quad_prob->con[c]=new SepQcFunc(*ext_quad_prob->con[c]);
		ext_quad_prob->obj=new SepQcFunc(*ext_quad_prob->obj);
	}
	ext_convex_prob=new MinlpProblem(*opt.convex_prob);
	for (int c=0; c<ext_convex_prob->con.size(); c++) ext_convex_prob->con[c]=new SepQcFunc(*ext_convex_prob->con[c]);
	ext_convex_prob->obj=new SepQcFunc(*ext_convex_prob->obj);

	// introducing t's
	vector<vector<Pointer<BlockMatrix> > > A_split(opt.split_prob->con.size()+1);
	vector<vector<Pointer<SepQcFunc> > > s_split(opt.split_prob->con.size()+1);
	vector<vector<Pointer<BlockMatrix> > > A_quad(opt.split_prob->con.size()+1);
	vector<vector<Pointer<SepQcFunc> > > s_quad(opt.split_prob->con.size()+1);
	vector<vector<Pointer<BlockMatrix> > > A_convex(opt.split_prob->con.size()+1);
	vector<vector<Pointer<SepQcFunc> > > s_convex(opt.split_prob->con.size()+1);
	for (int c=0; c<A_split.size(); c++) {
		A_split[c].resize(opt.split_prob->block.size());
		s_split[c].resize(opt.split_prob->block.size());
		if (opt.quad_prob) {
			A_quad[c].resize(opt.split_prob->block.size());
			s_quad[c].resize(opt.split_prob->block.size());
		}
		A_convex[c].resize(opt.split_prob->block.size());
		s_convex[c].resize(opt.split_prob->block.size());
	}
	for (int c=0; c<lower.size(); c++) {
		SepQcFunc& f(c ? *opt.split_prob->con[c-1] : *opt.split_prob->obj);
		bool lin_block=false; int nonlinblocks=0;
		for (int k=0; k<f.block.size(); k++) {
			if (f.s[k] || f.A[k]) nonlinblocks++;
			else if (f.b[k]) lin_block=true;
		}
		if ((!lin_block) && (nonlinblocks<=1) && (c || (!nonlinblocks))) continue;

		for (int k=0; k<f.block.size(); k++) {
			if ((!f.s[k]) && (!f.A[k])) continue;

			double initial=f.eval(opt.split_prob->primal_point(f.block[k]), k);
/*			out_log << "add var. for block " << k << " of ";
			if (c) { out_log << opt.split_prob->con_names[c-1]; }	else out_log << "objective";
			out_log << "\t lower bound: " << lower[c](k) << "\t initial: " << initial << endl;
*/

			Pointer<char> t_name=new char[10+(c ? strlen((char*)opt.split_prob->con_names[c-1]) : 10)];
			if (c) sprintf((char*)t_name, "t%s_%d", (char*)opt.split_prob->con_names[c-1], k);
			else sprintf((char*)t_name, "tobjective_%d", k);
			add_var(A_split, s_split, *ext_prob, k, ext_prob->dim(), lower[c](k), c ? upper[c](k) : obj_upper[k], t_name);
			if (ext_quad_prob) add_var(A_quad, s_quad, *ext_quad_prob, k, ext_quad_prob->dim(), lower[c](k), c ? upper[c](k) : obj_upper[k], t_name);
			add_var(A_convex, s_convex, *ext_convex_prob, k, ext_convex_prob->dim(), lower[c](k), c ? upper[c](k) : obj_upper[k], t_name);

			ext_convex_prob->primal_point[ext_prob->dim()-1]=ext_prob->primal_point[ext_prob->dim()-1]=initial;
			if (ext_quad_prob) ext_quad_prob->primal_point[ext_prob->dim()-1]=initial;

			add_con(A_split, s_split, *ext_prob, c, k);
			if (ext_quad_prob) add_con(A_quad, s_quad, *ext_quad_prob, c, k);
			add_con(A_convex, s_convex, *ext_convex_prob, c, k);
			
			if (!c) {
				opt.minlpdata->obj.reformulation_constraints_lower[k]=ext_prob->con.size()-1;
			} else {
				if(opt.ineq_index[c-1]) {
					opt.minlpdata->con[c-1].reformulation_constraints_lower[k]=ext_prob->con.size()-2;
					opt.minlpdata->con[c-1].reformulation_constraints_upper[k]=ext_prob->con.size()-1;
				} else
					opt.minlpdata->con[c-1].reformulation_constraints_lower[k]=ext_prob->con.size()-1;
			}			

			related_t.push_back(ext_prob->block[k].size()-1);
			if (c && opt.ineq_index[c-1]) related_t.push_back(ext_prob->block[k].size()-1);

			if (!c) {
				ext_convex_prob->con.back()->c = opt.conv_obj_c_add(k);
				ext_convex_prob->obj->c-=opt.conv_obj_c_add(k);
#ifndef FILIB_AVAILABLE
				ext_convex_prob->lower[ext_convex_prob->lower.dim()-1] += opt.conv_obj_c_add(k);
#endif
				if (ext_quad_prob) {
					ext_quad_prob->con.back()->c = opt.quad_obj_c_add(k);
					ext_quad_prob->obj->c-=opt.quad_obj_c_add(k);
#ifndef FILIB_AVAILABLE
					ext_quad_prob->lower[ext_quad_prob->lower.dim()-1] += opt.quad_obj_c_add(k);
#endif
				}
			} else { // adjust constant parts
				if (opt.ineq_index[c-1]) {
					ext_convex_prob->con[ext_convex_prob->con.size()-2]->c = opt.conv_con_c_add[c-1](k);
					ext_convex_prob->con[c-1]->c -= opt.conv_con_c_add[c-1](k);
					ext_convex_prob->con.back()->c = opt.conv_con_c_add[opt.ineq_index[c-1]](k);
					ext_convex_prob->con[opt.ineq_index[c-1]]->c -= opt.conv_con_c_add[opt.ineq_index[c-1]](k);
#ifndef FILIB_AVAILABLE
					ext_convex_prob->lower[ext_convex_prob->lower.dim()-1] += opt.conv_con_c_add[c-1](k);
#endif
					if (ext_quad_prob) {
						ext_quad_prob->con[ext_quad_prob->con.size()-2]->c = opt.quad_con_c_add[c-1](k);
						ext_quad_prob->con[c-1]->c -= opt.quad_con_c_add[c-1](k);
						ext_quad_prob->con.back()->c = opt.quad_con_c_add[opt.ineq_index[c-1]](k);
						ext_quad_prob->con[opt.ineq_index[c-1]]->c -= opt.quad_con_c_add[opt.ineq_index[c-1]](k);
#ifndef FILIB_AVAILABLE
						ext_quad_prob->lower[ext_quad_prob->lower.dim()-1] += opt.quad_con_c_add[c-1](k);
#endif
					}
				} else {
					ext_convex_prob->con.back()->c = opt.conv_con_c_add[c-1](k);
					ext_convex_prob->con[c-1]->c -= opt.conv_con_c_add[c-1](k);
#ifndef FILIB_AVAILABLE
					ext_convex_prob->lower[ext_convex_prob->lower.dim()-1] += opt.conv_con_c_add[c-1](k);
#endif
					if (ext_quad_prob) {
						ext_quad_prob->con.back()->c = opt.quad_con_c_add[c-1](k);
						ext_quad_prob->con[c-1]->c -= opt.quad_con_c_add[c-1](k);
#ifndef FILIB_AVAILABLE
						ext_quad_prob->lower[ext_quad_prob->lower.dim()-1] += opt.quad_con_c_add[c-1](k);
#endif
					}
				}
			}
		}
	}

	int deleted_constraints=0;
	for (int c=0; c<opt.ineq_index.size(); c++) {
		if (opt.ineq_index[c]<=0) continue;
		ext_prob->del_con(opt.ineq_index[c]);
		ext_prob->con_eq[c]=true;
		if (ext_quad_prob) {
			ext_quad_prob->con_eq[c]=true;
			ext_quad_prob->del_con(opt.ineq_index[c]);
		}
		ext_convex_prob->del_con(opt.ineq_index[c]);
 		// make equality con only, if we added t's here; if it had the form g(x)+y<=0, where g is defined for one block, it needs to remain as inequality
		ext_convex_prob->con_eq[c]=related_t[c]<0 ? false : true;

		for (int d=c+1; d<opt.ineq_index.size(); d++)
			if (opt.ineq_index[d]>opt.ineq_index[c]) opt.ineq_index[d]--;

		for (int d=opt.ineq_index[c]; d<related_t.size()-1; d++) // move t's forward
			related_t[d]=related_t[d+1];

		opt.ineq_index[c]=-1;
//		opt.minlpdata->con[c].ineq_index=0;
		++deleted_constraints;
	}
	related_t.resize(ext_prob->con.size(), -1);
	
	for (map<int,int>::iterator it(opt.minlpdata->obj.reformulation_constraints_lower.begin()); it!=opt.minlpdata->obj.reformulation_constraints_lower.end(); ++it)
		it->second-=deleted_constraints;
	for (int d=0; d<opt.minlpdata->con.size(); ++d) {
		for (map<int,int>::iterator it(opt.minlpdata->con[d].reformulation_constraints_lower.begin()); it!=opt.minlpdata->con[d].reformulation_constraints_lower.end(); ++it)
			it->second-=deleted_constraints;
		for (map<int,int>::iterator it(opt.minlpdata->con[d].reformulation_constraints_upper.begin()); it!=opt.minlpdata->con[d].reformulation_constraints_upper.end(); ++it)
			it->second-=deleted_constraints;
	}

	out_log << "Reformulation time: " << t.stop() << endl;
	out_log << "Initial point violates " << ext_prob->feasible(ext_prob->primal_point, 1E-4, NULL) << " constraints in splitting formulation; Objective value: " << ext_prob->obj->eval(ext_convex_prob->primal_point) << endl;
	out_log << "Initial point violates " << ext_convex_prob->feasible(ext_convex_prob->primal_point, 1E-4, NULL) << " constraints in convexified problem; Objective value: " << ext_convex_prob->obj->eval(ext_convex_prob->primal_point) << endl;
	out_log << "Initial point violates (1=true, 0=false) lower bound: " << (ext_convex_prob->lower-ext_convex_prob->primal_point>=rtol) << " upper bound: " << (ext_convex_prob->primal_point-ext_convex_prob->upper>=rtol) << endl;
}

dvector Reformulation::get_short_vector(const dvector &x) {
	return x(0, opt.orig_prob->dim()-1);
}

dvector Reformulation::get_long_vector(const dvector &x) {
	dvector y(ext_prob->dim());
	for (int i=0; i<x.dim(); i++) y[i]=x(i);
	for (int c=0; c<related_t.size(); c++) { // because of level cut
		if (related_t[c]<0) continue;
		for (int k=0; k<ext_prob->block.size(); k++) {
			if (ext_prob->con[c]->A[k] || ext_prob->con[c]->b[k] || ext_prob->con[c]->s[k]) {
				if (y[ext_prob->block[k][related_t[c]]]) break; // already set before
				assert(ext_prob->con[c]->b[k] && ((*ext_prob->con[c]->b[k])(related_t[c])));
				y[ext_prob->block[k][related_t[c]]]=-ext_prob->con[c]->eval(y)/(2*(*ext_prob->con[c]->b[k])(related_t[c]));
				break;
			}
		}
	}
	return y;
}

bool Reformulation::set_t_block(dvector& x, int k) {
//	for (int i=opt.split_prob->block[k].size(); i<ext_prob->block[k].size(); i++) x[i]=0.;

	for (int c=0; c<related_t.size(); c++) { // because of level cut
		int tind=related_t[c];
		if (tind<0) continue; // no t
		if (!(ext_prob->con[c]->A[k] || ext_prob->con[c]->b[k] || ext_prob->con[c]->s[k])) continue; // linear constraint
//		if (x[tind]) continue; // already set before
		double val=ext_prob->con[c]->eval(x, k)+ext_prob->con[c]->c;
		if (!finite(val)) return false;
		if (val<1E-4 && (val>-1E-4 || !ext_prob->con_eq[c])) continue; // constraint is satisfied
		assert(ext_prob->con[c]->b[k] && ((*ext_prob->con[c]->b[k])(tind)));
		x[tind]-=val/(2*(*ext_prob->con[c]->b[k])(tind));
		if (x[tind]<ext_prob->lower[ext_prob->block[k][tind]]) x[tind]=ext_prob->lower[ext_prob->block[k][tind]];
		else if (x[tind]>ext_prob->upper[ext_prob->block[k][tind]]) x[tind]=ext_prob->upper[ext_prob->block[k][tind]];
	}

	return true;
}

// ----------------------------- LevelCutHandling ------------------------------

void LevelCutHandler::add_problem(Pointer<MinlpProblem> prob) {
	int pos=prob->con.size();
	prob->add_con(new SepQcFunc(*prob->obj), false, "levelcut");
	prob->con.back()->c-=val;
	problems_with_levelcut.push_back(pair<Pointer<MinlpProblem>, int>(prob, pos));
}

void LevelCutHandler::add_problem(Pointer<LinearRelax> linear_relax) {
	linear_relax->add_level_cut(val);
	linrelax_with_levelcut.push_back(linear_relax);
}

void LevelCutHandler::update_level_cut(double newval) {
	if (newval>val) return;
	val=newval;
	for (list<pair<Pointer<MinlpProblem>, int> >::iterator it(problems_with_levelcut.begin()); it!=problems_with_levelcut.end(); it++)
		it->first->con[it->second]->c=it->first->obj->c-val;
	for (list<Pointer<LinearRelax> >::iterator it(linrelax_with_levelcut.begin()); it!=linrelax_with_levelcut.end(); it++)
		(*it)->update_level_cut(val);
}
		
LevelCutHandler::~LevelCutHandler() { };

// ------------------------------------ MinlpOpt -------------------------------

MinlpOpt::MinlpOpt(Pointer<MinlpProblem> prob, Pointer<Param> param_, bool is_gams_prob_, Pointer<ostream> out_solver_p, Pointer<ostream> out_solver_log_p)
: Solver(prob->dim(), out_solver_p, out_solver_log_p),
	orig_prob(prob), param(param_), is_gams_prob(is_gams_prob_), ineq_index(prob->con.size()), prob_is_convex(false), boxreduce_time(0.), low_bound(-INFINITY), timer(new Timer())
{	tol=1E-4;
	boxreduce_effort=param->get_i("Boxreduce effort", 1);

	char* val=param->get("Boxreduce type", "off");
	if ((!val) || (!strcmp(val, "NLP"))) boxreduce_type=box_C;
	else if (!strcmp(val, "NLP2")) boxreduce_type=box_Cext;
	else if (!strcmp(val, "MINLP")) boxreduce_type=box_Pext;
	else if (!strcmp(val, "off")) boxreduce_type=box_off;
	else {
		out_err << "Boxreduce type " << val << " not known. Try other. Aborting." << endl;
		exit(-1);
	}
	if (val) free(val);
}

void MinlpOpt::decompose() {
	if (param->get_i("Decomposition", 1)) {
		vector<vector<dvector> > sample_set;
		sample_set.resize(orig_prob->block.size());
		if (!param->get("Decomposition sample set Monte Carlo"))
			param->add("Decomposition sample set Monte Carlo", "20");
		if (!param->get("Decomposition sample set mid point"))
			param->add("Decomposition sample set mid point", "1");
		if (!param->get("Decomposition sample set vertices"))
			param->add("Decomposition sample set vertices", "20");
		Sampling(param, "Decomposition").get_points(sample_set, orig_prob->lower, orig_prob->upper, orig_prob->block);

		split_prob=Decomposition(param).decompose(*orig_prob, sample_set);
		if (out_log_p && param->get_i("Check decomposition", 0)) { // with Decomposition check
			out_log << "Checking decomposition..." << endl;
/*	 		dvector x1(orig_prob->dim());
 			dvector x2(split_prob->dim());
	 		for (int i=0; i<20; i++) {
	 			x1.set_random(orig_prob->lower, orig_prob->upper);
 				ss->old2new(x2, x1);
				double rel=orig_prob->obj->eval(x1);
 				double diff=(rel-split_prob->obj->eval(x2))/(1+fabs(rel));
   			if (diff>rtol) out_log << "  |old obj - new obj|: " << diff << endl;
  		}
		 	for (int c=0; c<orig_prob->con.size(); c++) {
				double diff=0.;
				int nr=0;
				for (int i=0; i<20; i++) {
					x1.set_random(orig_prob->lower, orig_prob->upper);
					ss->old2new(x2, x1);
					double rel=orig_prob->con[c]->eval(x1);
					double difff=(rel-split_prob->con[c]->eval(x2))/(1+fabs(rel));
					if (fabs(difff)>rtol) nr++;
			  	diff+=difff;
				}
				if (nr) {
					out_log << c << ": " << orig_prob->con_names[c] << ": |old - new|/(1+|old|): " << diff << "\t " << nr << endl;
					out_log << "old: "; orig_prob->con[c]->print(*out_log_p, orig_prob->var_names);
					out_log << "new: "; split_prob->con[c]->print(*out_log_p, split_prob->var_names);
				}
			}
*/
			dvector x(orig_prob->dim());
	 		for (int i=0; i<20; i++) {
	 			x.set_random(orig_prob->lower, orig_prob->upper);
				double rel=orig_prob->obj->eval(x);
 				double diff=(rel-split_prob->obj->eval(x))/(1+fabs(rel));
   			if (diff>rtol) out_log << "  |old obj - new obj|: " << diff << endl;
  		}
		 	for (int c=0; c<orig_prob->con.size(); c++) {
				double diff=0.;
				int nr=0;
				for (int i=0; i<20; i++) {
					x.set_random(orig_prob->lower, orig_prob->upper);
					double rel=orig_prob->con[c]->eval(x);
					double difff=(rel-split_prob->con[c]->eval(x))/(1+fabs(rel));
					if (fabs(difff)>rtol) nr++;
			  	diff+=difff;
				}
				if (nr) {
					out_log << c << ": " << orig_prob->con_names[c] << ": |old - new|/(1+|old|): " << diff << "\t " << nr << endl;
					out_log << "old: "; orig_prob->con[c]->print(*out_log_p, orig_prob->var_names);
					out_log << "new: "; split_prob->con[c]->print(*out_log_p, split_prob->var_names);
				}
			}
		}

	} else split_prob=orig_prob;
/*
	int maxk=0; int maxsize=0;
	for (int k=0; k<split_prob->block.size(); k++) {
		if (maxsize<split_prob->block[k].size()) { maxsize=split_prob->block[k].size(); maxk=k; }

		if (split_prob->block[k].size()<=11) continue;
		out_log << endl << "block " << k << " has " << split_prob->block[k].size() << " variables:" << endl;
	for (int i=0; i<split_prob->block[k].size(); i++)
		out_log << split_prob->var_names[split_prob->block[k][i]] << endl;
	out_log << endl << "used in constraints " << endl;
	for (int c=0; c<split_prob->con.size(); c++)
		if (split_prob->con[c]->A[k] || split_prob->con[c]->s[k]) out_log << split_prob->con_names[c] << endl;

	}
//	out_log << endl << "biggest block " << maxk << " has " << maxsize << " variables: " << endl;
*/
	conv_obj_c_add.resize(split_prob->block.size());
	conv_obj_c_add=0;

	intervalreduction.set_problem(orig_prob);
}

bool MinlpOpt::check_convex(MinlpProblem& prob) {
	min_eigval.resize(prob.con.size()+1);
	max_eigval.resize(prob.con.size()+1);

	vector<vector<dvector> > sample_set(prob.block.size());
	if (!param->get("Relax check convex sample set Monte Carlo"))
		param->add("Relax check convex sample set Monte Carlo", "200");
	Sampling(param, "Relax check convex").get_points(sample_set, prob.lower, prob.upper, prob.block);

	Convexify convexify(param);

	out_solver << "Checking problem for convexity ";
	Timer t;
	bool allconvex=true;
	for (int c=0; c<=prob.con.size(); c++) {
		out_solver_log << '.';
//		if (c) out_solver_log << prob.con_names[c-1] << endl; 
		convexify.check_convex2(min_eigval[c], max_eigval[c], c ? *prob.con[c-1] : *prob.obj, sample_set);
		Func::CurvatureType ct=c ? prob.con[c-1]->get_curvature() : prob.obj->get_curvature();
		if ((!(ct&Func::CONVEX)) || (c && prob.con_eq[c-1] && !(ct&Func::CONCAVE))) allconvex=false;
	}
	out_solver << " done. time: " << t.stop() << endl;

	return allconvex;
}

void MinlpOpt::quad_relax() {
	quad_obj_c_add=0;
	quad_obj_c_add.resize(split_prob->block.size());
	quad_con_c_add.clear();
	quad_con_c_add.resize(split_prob->con.size(), SparseVector<double>(split_prob->block.size()));

	quad_prob=new MinlpProblem(*split_prob);
//	quad_prob->obj=new SepQcFunc(*quad_prob->obj);
//	for (int c=0; c<quad_prob->con.size(); c++) quad_prob->con[c]=new SepQcFunc(*quad_prob->con[c]);

	Timer t;
	if (param->get_i("Quadratic Underestimator adaptive", 0)) {
		if (!param->get("Quadratic Underestimator sample set Monte Carlo"))
			param->add("Quadratic Underestimator sample set Monte Carlo", "20");
		if (!param->get("Quadratic Underestimator sample set mid point"))
			param->add("Quadratic Underestimator sample set mid point", "0");
		if (!param->get("Quadratic Underestimator sample set box ends"))
			param->add("Quadratic Underestimator sample set box ends", "0");
		if (!param->get("Quadratic Underestimator sample set vertices2"))
			param->add("Quadratic Underestimator sample set vertices2", "200");
		if (!param->get("Quadratic Underestimator sample set minimizer"))
			param->add("Quadratic Underestimator sample set minimizer", LocOpt::nlp_solver_available() ? "1" : "0");	
	
		QuadraticUnderestimator quaduest(param);
		quaduest.quadratic_underestimator(*quad_prob, *minlpdata, ineq_index, quad_obj_c_add, quad_con_c_add);
		if (t.stop()>rtol) {
			out_log << "Time for U3: " << ((int)(1000*quaduest.U3_time/t))/10. << "\\%\t";
			out_log << "Time for locopt: " << ((int)(1000*quaduest.locopt_time/t))/10. << "\\%\t";
		}
		out_log << "Maximal coefficient: " << quaduest.max_abscoeff << '\t';
		out_log << "Maximal locmin: " << quaduest.max_locmin << endl;
	} else {
		if (!param->get("Polynomial Underestimator K0 sample set Monte Carlo"))
			param->add("Polynomial Underestimator K0 sample set Monte Carlo", "200");
		if (!param->get("Polynomial Underestimator K0 sample set mid point"))
			param->add("Polynomial Underestimator K0 sample set mid point", "0");
		if (!param->get("Polynomial Underestimator K0 sample set box ends"))
			param->add("Polynomial Underestimator K0 sample set box ends", "0");
		if (!param->get("Polynomial Underestimator K0 sample set vertices2"))
			param->add("Polynomial Underestimator K0 sample set vertices2", "200");
		if (!param->get("Polynomial Underestimator K0 sample set minimizer"))
			param->add("Polynomial Underestimator K0 sample set minimizer", LocOpt::nlp_solver_available() ? "1" : "0");
		if (!param->get("Polynomial Underestimator K1 sample set Monte Carlo"))
			param->add("Polynomial Underestimator K1 sample set Monte Carlo", "10");
		if (!param->get("Polynomial Underestimator K1 sample set mid point"))
			param->add("Polynomial Underestimator K1 sample set mid point", "1");
		if (!param->get("Polynomial Underestimator K1 sample set box ends"))
			param->add("Polynomial Underestimator K1 sample set box ends", "0");
		if (!param->get("Polynomial Underestimator K1 sample set vertices"))
			param->add("Polynomial Underestimator K1 sample set vertices", "0");
		if (!param->get("Polynomial Underestimator K2 sample set Monte Carlo"))
			param->add("Polynomial Underestimator K2 sample set Monte Carlo", "10");
		if (!param->get("Polynomial Underestimator K2 sample set mid point"))
			param->add("Polynomial Underestimator K2 sample set mid point", "0");
		if (!param->get("Polynomial Underestimator K2 sample set box ends"))
			param->add("Polynomial Underestimator K2 sample set box ends", "0");
		if (!param->get("Polynomial Underestimator K2 sample set vertices"))
			param->add("Polynomial Underestimator K2 sample set vertices", "0");
			
		PolynomialUnderestimator2 polyuest(param);
		polyuest.polynomial_underestimator(*quad_prob, *minlpdata, ineq_index, quad_obj_c_add, quad_con_c_add);
	
		if (param->get_i("Check polynomial underestimator", 0))
			polyuest.check(*split_prob, *quad_prob, ineq_index);
	}
	out_log << "Time for quadratic underestimators: " << t.stop() << endl;

	conv_obj_c_add=quad_obj_c_add;
	conv_con_c_add=quad_con_c_add;


	for (int c=0; c<ineq_index.size(); c++) {
		if (ineq_index[c]) {
			assert(ineq_index[c]==split_prob->con.size());
			split_prob->add_con(new SepQcFunc(*split_prob->con[c], true), false, quad_prob->con_names[ineq_index[c]]);
			split_prob->con_eq[c]=false;
			min_eigval.push_back(max_eigval[c]);
			max_eigval.push_back(min_eigval[c]);
			for (int k=0; k<min_eigval.back().size(); k++) {
			  min_eigval.back()[k]*=-1.;
			  max_eigval.back()[k]*=-1.;
			}
		}
	}

	prob_is_convex=check_convex(quad_prob ? *quad_prob : *split_prob);
	out_out << "Quadratic problem is convex: " << prob_is_convex << endl;
}

void MinlpOpt::convex_relax() {
	// The problem, we convexify. Taking the best existing.
	MinlpProblem& prob(quad_prob ? *quad_prob : *split_prob);

	convex_prob=new MinlpProblem(prob);
	for (vector<int>::iterator it(convex_prob->i_discr.begin()); it!=convex_prob->i_discr.end(); it++) {
		convex_prob->i_cont.push_back(*it);
		convex_prob->discr[*it]=false;
	}
	convex_prob->i_discr.clear();

	conv_con_c_add.clear();
	if (!quad_con_c_add.size())
		for (int c=0; c<prob.con.size(); c++) conv_con_c_add.push_back(SparseVector<double>(prob.block.size()));
	else
		conv_con_c_add=quad_con_c_add;
	if (!quad_obj_c_add.size()) conv_obj_c_add.resize(prob.block.size());
	else conv_obj_c_add=quad_obj_c_add;

	int convex_nr=0; // the number of convexifications, we computed

	convex_prob->obj=new SepQcFunc(*convex_prob->obj);
	for (int c=0; c<convex_prob->con.size(); c++) convex_prob->con[c]=new SepQcFunc(*convex_prob->con[c]);
	if (!prob_is_convex) {
		Convexify convexify(param);
		if (!quad_prob) ineq_index.resize(orig_prob->con.size());

		out_out << "Convexifying objective and constraints ";
		double max_alpha_norm=0.;
		int orig_con_size=prob.con.size();
		ivector ineq_index_inverse(orig_con_size, -1);
		for (int c=0; c<=orig_con_size; c++) {
			if (c>0 && c-1<minlpdata->con.size() && ineq_index[c-1]>=0) ineq_index_inverse[ineq_index[c-1]]=c-1;
			Pointer<SepQcFunc> f(c ? convex_prob->con[c-1] : convex_prob->obj);
			if ((f->get_curvature()&Func::CONVEX) && ((!c) || (!convex_prob->con_eq[c-1]) || (f->get_curvature()&Func::CONCAVE))) continue;

/*			if (c) { out_log << "Convexifying " << prob.con_names[c-1] << ": " << endl; }
			else out_log << "Convexifying objective: " << endl;*/
 			out_log << ".";

			convexify.new_sampleset(prob.lower, prob.upper, *f);
			double alpha_norm=convexify.convexify(*f, c ? convex_prob->con_eq[c-1] : false, min_eigval[c], max_eigval[c], prob.lower, prob.upper);
			pair<Pointer<SepQcFunc>, Pointer<SepQcFunc> > p(f, NULL);
			convexify.get_decomposed_functions(p);
			max_alpha_norm=MAX(max_alpha_norm, alpha_norm);

			if (c) {
				if (c-1>=minlpdata->con.size()) {
					minlpdata->con[ineq_index_inverse[c-1]].convexification_characteristica_upper=convexify.characteristica1;
					minlpdata->con[ineq_index_inverse[c-1]].concave_overestimator=p.first;
//					out_log << "store convexification information as upper part in " << minlpdata->con[ineq_index_inverse[c-1]].name << endl;
				} else {
					minlpdata->con[c-1].convexification_characteristica_lower=convexify.characteristica1;
					minlpdata->con[c-1].convex_underestimator=p.first;
				}
			} else {
				minlpdata->obj.convexification_characteristica_lower=convexify.characteristica1;
				minlpdata->obj.convex_underestimator=p.first;
			}
/*
			if (c) {
				prob.con[c-1]->print(*out_log_p, prob.var_names);
				p.first->print(*out_log_p, prob.var_names);
				if (p.second) p.second->print(*out_log_p, prob.var_names);
			}
*/

			for (int k=0; k<prob.block.size(); k++)
				if (convexify.convexify_c[k].first)
					(c ? conv_con_c_add[c-1][k] : conv_obj_c_add[k])+=convexify.convexify_c[k].first;

			if (c && convex_prob->con_eq[c-1]) { // equality constraint
				char* name=NULL;
				if (convex_prob->con_names[c-1])
					sprintf(name=new char[strlen(convex_prob->con_names[c-1])+2], "-%s", (char*)convex_prob->con_names[c-1]);
				ineq_index[c-1]=convex_prob->con.size();
//				minlpdata->con[c-1].ineq_index=ineq_index[c-1];
				conv_con_c_add.push_back(-quad_con_c_add[c-1]);
				if (!p.second) p.second=new SepQcFunc(*convex_prob->con[c-1], true); // concave function
				else {
					for (int k=0; k<prob.block.size(); k++)
						if (convexify.convexify_c[k].second)
							conv_con_c_add.back()[k]+=convexify.convexify_c[k].second;
					convex_nr++;
					assert(c-1<minlpdata->con.size());
					minlpdata->con[c-1].convexification_characteristica_upper=convexify.characteristica2;
					minlpdata->con[c-1].concave_overestimator=p.second;
				}

				convex_prob->add_con(p.second, false, name);
				convex_prob->con_eq[c-1]=false;
				split_prob->add_con(new SepQcFunc(*split_prob->con[c-1], true), false, name);
				split_prob->con_eq[c-1]=false;
				if (quad_prob) {
					quad_prob->add_con(new SepQcFunc(*quad_prob->con[c-1], true), false, name);
					quad_prob->con_eq[c-1]=false;
					quad_con_c_add.push_back(-quad_con_c_add[c-1]);
				}
			}
			(c ? convex_prob->con[c-1] : convex_prob->obj)=p.first;

			convex_nr++;
		}
		out_log << endl << "max |\\alpha| = " << max_alpha_norm << endl;
	}
	out_out << "Computed " << convex_nr << " convexifications." << endl;
	out_log << "Initial point violates " << convex_prob->feasible(convex_prob->primal_point, 1E-4, NULL) << " constraints in convexified problem; Objective value: " << convex_prob->obj->eval(convex_prob->primal_point) << endl;

	if (out_log_p && param->get_i("Check convexification", 0)) {
		out_log << "Checking convexification: " << endl;
		Convexify(param).check_convexification(*convex_prob, *split_prob, ineq_index, 20);
	}
}

double MinlpOpt::print_box_reduce_quality(dvector& oldlow, dvector& oldup, Pointer<MinlpProblem> prob, char* prefix) {
	double avg_rel_impr=0.; // average relative improvement
	double min_rel_impr=1.; // minimal relative improvement
	int nr_real_impr=0; // number of improvements, better than 20%
	int nr_bin_fixed=0; // number of fixed binary variables

	for (int i=0; i<oldlow.dim(); i++) {
		if (oldup(i)>=INFINITY || oldlow(i)<=-INFINITY) {
			if (prob->upper(i)<INFINITY && prob->lower(i)>-INFINITY) {
				min_rel_impr=0.;
				nr_real_impr++;
			} else
				avg_rel_impr+=1.;
			continue;
		}
		double rel_impr;
		if (oldup(i)-oldlow(i)<rtol) rel_impr=1.;
		else rel_impr=(prob->upper(i)-prob->lower(i)) / (oldup(i)-oldlow(i));
		avg_rel_impr+=rel_impr;
		if (rel_impr<min_rel_impr) min_rel_impr=rel_impr;
		if (rel_impr<=.8) nr_real_impr++;
		if (prob->discr[i] && (rel_impr<rtol)) nr_bin_fixed++;
	}
	avg_rel_impr/=oldlow.dim();

	out_out.precision(3);
	out_out << prefix << " average relative box diameter improvement: " << avg_rel_impr*100 << endl;
	out_out << prefix << " minimal relative box diameter improvement: " << min_rel_impr*100 << endl;
	out_out << prefix << " percentage number of real box improvements: " << (100*nr_real_impr)/oldlow.dim() << endl;
	if (prob->i_discr.size()) out_out << prefix << " percentage number of fixed binaries: " << (100*nr_bin_fixed/prob->i_discr.size()) << endl;
	out_out.precision(6);
/*
	out_log << "new box: " << endl;
	for (int i=0; i<prob->dim(); i++)
		out_log << prob->var_names[i] << ": \t" << prob->lower[i] << "\t " << prob->upper[i] << endl;
*/
	return avg_rel_impr;
}

Pointer<MinlpProblem> MinlpOpt::get_convex_prob(Pointer<MinlpProblem> prob) {
	Pointer<MinlpProblem> conv_prob(new MinlpProblem(*prob));
	for (int i=0; i<conv_prob->i_discr.size(); i++) {
		conv_prob->i_cont.push_back(conv_prob->i_discr[i]);
		conv_prob->discr[conv_prob->i_discr[i]]=false;
	} conv_prob->i_discr.clear();

	conv_prob->con.clear();
	conv_prob->con_eq.clear();
	conv_prob->con_names.clear();
	for (int c=0; c<prob->con.size(); c++) {
		// convex inequality constraint -> take it
		if (!prob->con_eq[c]) {
			if (prob->con[c]->get_curvature()&Func::CONVEX)
				conv_prob->add_con(prob->con[c], false, prob->con_names[c]);
		}	else {
			// convex, concave equality constraint -> take it
			if (prob->con[c]->get_curvature()==Func::LINEAR)
				conv_prob->add_con(prob->con[c], true, prob->con_names[c]);
			// convex, nonconcave equality constraint -> keep it as inequality
			else if (prob->con[c]->get_curvature()==Func::CONVEX)
				conv_prob->add_con(prob->con[c], false, prob->con_names[c]);
			// nonconvex, concave equality constraint -> change sign and keep it as inequality
			else if (prob->con[c]->get_curvature()==Func::CONCAVE)
				conv_prob->add_con(new SepQcFunc(*prob->con[c], true), false, prob->con_names[c]);
		}
	}

	return conv_prob;
}

Pointer<MinlpProblem> MinlpOpt::get_convex_prob(Pointer<MinlpProblem> prob, Pointer<MinlpProblem> prob_curv_ref) {
	Pointer<MinlpProblem> conv_prob(new MinlpProblem(*prob));
	for (int i=0; i<conv_prob->i_discr.size(); i++) {
		conv_prob->i_cont.push_back(conv_prob->i_discr[i]);
		conv_prob->discr[conv_prob->i_discr[i]]=false;
	}
	conv_prob->i_discr.clear();

	conv_prob->con.clear();
	conv_prob->con_eq.clear();
	conv_prob->con_names.clear();
	for (int c=0; c<prob->con.size(); c++) {
		// convex inequality constraint -> take it
		if (!prob->con_eq[c]) {
			if (prob_curv_ref->con[c]->get_curvature()&Func::CONVEX)
				conv_prob->add_con(prob->con[c], false, prob->con_names[c]);
		}	else {
			// convex, concave equality constraint -> take it
			if (prob_curv_ref->con[c]->get_curvature()==Func::LINEAR)
				conv_prob->add_con(prob->con[c], true, prob->con_names[c]);
			// convex, nonconcave equality constraint -> keep it as inequality
			else if (prob_curv_ref->con[c]->get_curvature()==Func::CONVEX)
				conv_prob->add_con(prob->con[c], false, prob->con_names[c]);
			// nonconvex, concave equality constraint -> change sign and keep it as inequality
			else if (prob_curv_ref->con[c]->get_curvature()==Func::CONCAVE)
				conv_prob->add_con(new SepQcFunc(*prob->con[c], true), false, prob->con_names[c]);
		}
	}

	return conv_prob;
}

void MinlpOpt::box_reduce0() {
	for (int i=0; i<orig_prob->dim(); ++i)
		if (orig_prob->lower(i)<=-INFINITY || orig_prob->upper(i)>=INFINITY) unbounded_var.push_back(i);
	out_out << "Unbounded variables: " << (100*unbounded_var.size())/orig_prob->dim() << "\\%" << endl;
	out_out << "Original box diameter: ";
	if (unbounded_var.empty()) { out_out << sqrt((orig_prob->upper-orig_prob->lower).sq_norm2()) << endl; }
	else out_out << "$\\infty$" << endl;

	out_out << "Boxreduction phase 0:";
	if (!boxreduce_effort) { out_out << " Only unknown bounds." << endl; }
	else out_out << endl;

	Timer t;
//	intervalreduction.do_print=true;
	intervalreduction.compute(split_prob->lower, split_prob->upper, split_prob->lower, split_prob->upper);
	for (list<int>::iterator it(unbounded_var.begin()); it!=unbounded_var.end();)
		if (split_prob->lower(*it)>-INFINITY && split_prob->upper(*it)<INFINITY) it=unbounded_var.erase(it);
		else ++it;

	set<pair<int, IntervalReduction::which_bound_type> > changed_var; // variable bounds, changed by BoundsFinderLinear
	if (boxreduce_effort || !unbounded_var.empty()) {
		BoundsFinderLinear bfl(split_prob, param);
		bfl.known=boxreduce_effort;
		dvector newlow(split_prob->dim()), newup(split_prob->dim());
		int ret=bfl.compute(newlow, newup, &changed_var);
		split_prob->lower=newlow; split_prob->upper=newup;
		if (ret) out_out << "BoundsFinder Linear Constraints problems: " << (100*ret)/(2*split_prob->dim()) << "\\%" << endl;
		for (list<int>::iterator it(unbounded_var.begin()); it!=unbounded_var.end();)
			if (split_prob->lower(*it)>-INFINITY && split_prob->upper(*it)<INFINITY) it=unbounded_var.erase(it);
			else ++it;
	}

	if (!changed_var.empty()) {
		intervalreduction.compute(split_prob->lower, split_prob->upper, split_prob->lower, split_prob->upper, changed_var);
		for (list<int>::iterator it(unbounded_var.begin()); it!=unbounded_var.end();)
			if (split_prob->lower(*it)>-INFINITY && split_prob->upper(*it)<INFINITY) it=unbounded_var.erase(it);
			else ++it;
	}

	boxreduce_time+=t.stop();
	print_box_reduce_quality(orig_prob->lower, orig_prob->upper, split_prob, "Boxreduction phase 0");
	out_log << "Boxreduction phase 0 time: " << t << endl;
	out_log << "Boxreduction phase 0 unbounded variables: " << (100*unbounded_var.size())/orig_prob->dim() << "\\%" << endl;
}

void MinlpOpt::box_reduce1() {
	out_out << "Boxreduction phase 1:";
	if (!boxreduce_effort) { out_out << " Only unknown bounds." << endl; }
	else out_out << endl;

	Timer t;
	set<pair<int, IntervalReduction::which_bound_type> > changed_var;
	if (boxreduce_effort || !unbounded_var.empty()) {
		Pointer<MinlpProblem> conv;
		if (prob_is_convex) {
			conv=new MinlpProblem(*orig_prob);
			for (int i=0; i<conv->i_discr.size(); i++) {
				conv->i_cont.push_back(conv->i_discr[i]);
				conv->discr[conv->i_discr[i]]=false;
			}
			conv->i_discr.clear();
		} else conv=get_convex_prob(orig_prob, split_prob);
		conv->lower=split_prob->lower;
		conv->upper=split_prob->upper;
		
//		out_log << *conv;

		assert(orig_prob->block.size()==1);
		LinearRelax linrelax(param);
		Pointer<dvector> ref(new dvector(conv->primal_point));
		linrelax.init(conv, split_prob->i_discr, ref, false);

		IntervalGradientCutGenerator cutgen(orig_prob);
		LinearizedConCutGenerator cutgen2(conv);
		ref->set_random(conv->lower, conv->upper);
		linrelax.add_cut(cutgen.get_cuts(*ref, 0, conv->lower, conv->upper), 0);
		for (int c=0; c<conv->con.size(); ++c)
			if (conv->con[c]->A[0] || conv->con[c]->s[0]) {
				LinearizationCut cut(cutgen2.get_cut(*ref, c, 0, INFINITY));
				if (cut.coeff) linrelax.add_cut(Pointer<LinearizationCut>(new LinearizationCut(cut)), 0);
			}

		double red=linrelax.box_reduce(split_prob->lower, split_prob->upper, split_prob->discr, !boxreduce_effort, &changed_var);
		out_log << "LinearRelax boxreduction reduction: " << red << endl;

		for (list<int>::iterator it(unbounded_var.begin()); it!=unbounded_var.end();)
			if (split_prob->lower(*it)>-INFINITY && split_prob->upper(*it)<INFINITY) it=unbounded_var.erase(it);
			else ++it;
	}

	if (!changed_var.empty()) {
		intervalreduction.compute(split_prob->lower, split_prob->upper, split_prob->lower, split_prob->upper, changed_var);
		for (list<int>::iterator it(unbounded_var.begin()); it!=unbounded_var.end();)
			if (split_prob->lower(*it)>-INFINITY && split_prob->upper(*it)<INFINITY) it=unbounded_var.erase(it);
			else ++it;
	}

	if (boxreduce_effort>1 || !unbounded_var.empty()) {
		Pointer<MinlpProblem> conv;
		if (prob_is_convex) {
			if (!split_prob->i_discr.size()) conv=split_prob;
			else {
				conv=new MinlpProblem(*split_prob);
				for (int i=0; i<conv->i_discr.size(); i++) {
					conv->i_cont.push_back(conv->i_discr[i]);
					conv->discr[conv->i_discr[i]]=false;
				}
				conv->i_discr.clear();
			}
		} else conv=get_convex_prob(split_prob);

		BoundsFinder boundsfinder(param);
		boundsfinder.known=boxreduce_effort>1;
		pair<int,int> ret=boundsfinder.compute_bounds(*conv, split_prob->discr);
		split_prob->lower=conv->lower; split_prob->upper=conv->upper;
		if (ret.first) out_out << "BoundsFinder boxreduction problems: " << (100*ret.first)/(2*conv->dim()) << "\\%" << endl;
		if (ret.second) out_out << "Guessed bounds: " << (100*ret.second)/conv->dim() << "\\%" << endl;

		for (list<int>::iterator it(unbounded_var.begin()); it!=unbounded_var.end();)
			if (split_prob->lower(*it)>-INFINITY && split_prob->upper(*it)<INFINITY) it=unbounded_var.erase(it);
			else ++it;
		for (list<int>::iterator it(unbounded_var.begin()); it!=unbounded_var.end(); ++it)
			clog << split_prob->var_names[*it] << endl;
		assert(unbounded_var.empty());
	}

	boxreduce_time+=t.stop();
	print_box_reduce_quality(orig_prob->lower, orig_prob->upper, split_prob, "Boxreduction phase 1");
	out_log << "Boxreduction phase 1 time: " << t << endl;

	intervalreduction.print_small_boxes(split_prob->lower, split_prob->upper);
}

double MinlpOpt::box_reduce2() {
	if (boxreduce_effort<=1) return 1.;

	out_out << "Boxreduction phase 2:" << endl;

	dvector oldlow(split_prob->lower);
	dvector oldup(split_prob->upper);

	switch (boxreduce_type) {
		case box_C: {
			assert(convex_prob);
			Timer t;
			pair<int,int> ret=BoundsFinder(param).compute_bounds(*convex_prob, split_prob->discr);
			boxreduce_time+=t.stop();
			out_out << "Boxreduction phase 2 time: " << t << endl;
			if (ret.first) out_log << "BoundsFinder boxreduction problems: " << (100*ret.first)/(2*convex_prob->dim()) << "\\%" << endl;
			split_prob->lower=convex_prob->lower;
			split_prob->upper=convex_prob->upper;

			// lower bound
			t.start();
			Pointer<LocOpt> locopt=LocOpt::get_solver(convex_prob, param, "ConvexSolve", NULL, NULL);
			int ret2=locopt->solve(convex_prob->primal_point);
			out_out << "Lower bound time: " << t.stop() << endl;
			out_out << "Solving convex problem (C): return " << ret2 << "\t value: " << locopt->opt_val() << endl;
			if (ret2==0) low_bound=locopt->opt_val();
		} break;
		case box_Cext:
		case box_Pext: {
			if (!reform) {
				reform=new Reformulation(*this);
				reform->reformulate();
			}
			Timer t;
			DualBounds dualbounds(reform->ext_prob, reform->ext_convex_prob, param, boxreduce_type==box_Cext ? 1 : 2);
			int ret=dualbounds.update_box(split_prob->dim());
			boxreduce_time+=t.stop();
			out_out << "Boxreduction phase 2 time: " << t << endl;
			out_log << "Second boxreduction problems: " << (100*ret)/(2*split_prob->dim()) << "\\%" << endl;
			t.start();
			low_bound=dualbounds.obj_bound();
			out_out << "Lower bound time: " << t.stop() << endl;

			split_prob->upper=convex_prob->upper=reform->get_short_vector(reform->ext_prob->upper);
			split_prob->lower=convex_prob->lower=reform->get_short_vector(reform->ext_prob->lower);
			if (reform->ext_quad_prob) {
				reform->ext_quad_prob->upper=reform->ext_prob->upper;
				reform->ext_quad_prob->lower=reform->ext_prob->lower;
			}
		} break;
	}
	if (quad_prob) {
		quad_prob->lower=convex_prob->lower;
		quad_prob->upper=convex_prob->upper;
	}

	out_out.precision(20);
	out_out.setf(ios::fixed);
	if (low_bound>-INFINITY) out_out << "Lower bound: " << low_bound << endl;
	out_out.unsetf(ios::fixed);
	out_out.precision(6);
	print_box_reduce_quality(orig_prob->lower, orig_prob->upper, split_prob, "Boxreduction phase 2");

	double avg_rel_impr=0.;
	for (int i=0; i<oldlow.dim(); i++)
		if (oldup(i)-oldlow(i)<rtol) avg_rel_impr+=1.;
		else avg_rel_impr+=(split_prob->upper(i)-split_prob->lower(i)) / (oldup(i)-oldlow(i));
	return avg_rel_impr/oldlow.dim();
}

void MinlpOpt::box_reduce3() {
	if (!boxreduce_effort) {  // no more boxreduction
		sol_cand_diam=new dvector(reform->ext_convex_prob->upper);
		*sol_cand_diam-=reform->ext_convex_prob->lower;
		*sol_cand_diam*=param->get_d("heu close points tolerance", .001);
		return;
	}

	set<pair<int, IntervalReduction::which_bound_type> > changed_var;
	Timer t;
	double red=linear_relax->box_reduce(reform->ext_prob->lower, reform->ext_prob->upper, reform->ext_prob->discr, false, &changed_var);
	out_log << "LinearRelax boxreduction reduction: " << red << endl;
	split_prob->lower=reform->ext_prob->lower(0, split_prob->dim()-1);
	split_prob->upper=reform->ext_prob->upper(0, split_prob->dim()-1);

	if (!changed_var.empty()) {
		t.start();
		intervalreduction.compute(reform->ext_prob->lower, reform->ext_prob->upper, reform->ext_prob->lower, reform->ext_prob->upper, changed_var);
		split_prob->lower=reform->ext_prob->lower(0, split_prob->dim()-1);
		split_prob->upper=reform->ext_prob->upper(0, split_prob->dim()-1);
	}

	boxreduce_time+=t.stop();
	out_log << "Boxreduction phase 3 time: " << t << endl;

	// updating box
	convex_prob->lower=split_prob->lower;
	convex_prob->upper=split_prob->upper;
	if (quad_prob) {
		quad_prob->lower=convex_prob->lower;
		quad_prob->upper=convex_prob->upper;
	}
	reform->ext_convex_prob->lower=reform->ext_prob->lower;
	reform->ext_convex_prob->upper=reform->ext_prob->upper;
	if (reform->ext_quad_prob) {
		reform->ext_quad_prob->lower=reform->ext_prob->lower;
		reform->ext_quad_prob->upper=reform->ext_prob->upper;
	}
	linear_relax->set_box(reform->ext_prob->lower, reform->ext_prob->upper);
	print_box_reduce_quality(orig_prob->lower, orig_prob->upper, split_prob, "Boxreduction phase 3");

	// check, if box is complete now:
	for (int i=0; i<reform->ext_convex_prob->dim(); i++) {
		if (reform->ext_convex_prob->lower(i)==-INFINITY) {
			out_err << "Lower bound of variable " << i << " " << reform->ext_convex_prob->var_names[i] << " still not known, arborting." << endl;
			exit(-1);
		}
		if (reform->ext_convex_prob->upper(i)==INFINITY) {
			out_err << "Upper bound of variable " << i << " " << reform->ext_convex_prob->var_names[i] << " still not known, arborting." << endl;
			exit(-1);
		}
	}

	sol_cand_diam=new dvector(reform->ext_convex_prob->upper);
	*sol_cand_diam-=reform->ext_convex_prob->lower;
	*sol_cand_diam*=param->get_d("heu close points tolerance", .001);
}

void MinlpOpt::init() {
	check_initial_point();

	decompose();

	box_reduce0();

	prob_is_convex=check_convex(*split_prob); // checks the splitted problem for convexity.
	out_out << "Splitted problem convex: " << (prob_is_convex ? "y" : "n") << endl;
	for (int c=0; c<orig_prob->con.size(); ++c)
		orig_prob->con[c]->set_curvature(split_prob->con[c]->get_curvature());

	box_reduce1();

	minlpdata=new MINLPData(*split_prob);
	MINLP minlp(*minlpdata);
//	out_log << minlp;

/*	for (int i=0; i<split_prob->dim(); ++i)
		out_log << split_prob->var_names[i] << ": \t" << split_prob->lower[i] << '\t' << split_prob->upper[i] << endl;*/
	
	quad_relax();

	Pointer<MinlpProblem> split_prob_orig(new MinlpProblem(*split_prob));
	Pointer<MinlpProblem> quad_prob_orig(new MinlpProblem(*quad_prob));

	convex_relax();

//	out_log << minlpdata->obj;
//	for (int c=0; c<minlpdata->con.size(); ++c)
//		out_log << "Constraint " << minlpdata->con[c].name << endl << minlpdata->con[c];

 	if (boxreduce_type!=box_off)	{ // updating box 2nd time und compute lower bound
		double quality=box_reduce2();
		split_prob_orig->lower=split_prob->lower;
		split_prob_orig->upper=split_prob->upper;
		if (quad_prob) {
			quad_prob_orig->lower=quad_prob->lower;
			quad_prob_orig->upper=quad_prob->upper;
		}

		if (quality<param->get_d("Boxreduction limit for reconvexification", 0.8)) { // compute new convex relaxation
			out_out << "Computing new convex relaxation" << endl;
			convex_prob=NULL;
			split_prob=split_prob_orig;
			quad_prob=quad_prob_orig;
			check_convex(quad_prob ? *quad_prob : *split_prob);
			convex_relax();
			reform=NULL;
		}
	}

	// setting tolerance parameters for SolCandidate-set
	sol_cand_closeval_tol=param->get_d("heu close value tolerance", 10e-8);
	sol_cand_diam=new dvector(convex_prob->upper); *sol_cand_diam-=convex_prob->lower; // these bounds should exist
	*sol_cand_diam*=param->get_d("heu close points tolerance", .0001);

	out_out << "Relaxation finished." << endl << endl;
}

void MinlpOpt::init2() {
	if (!reform) {
		reform=new Reformulation(*this);
		reform->reformulate();

		intervalreduction.set_problem(reform->ext_prob);
		Timer t;
		intervalreduction.compute(reform->ext_prob->lower, reform->ext_prob->upper, reform->ext_prob->lower, reform->ext_prob->upper);
		boxreduce_time+=t.stop();
		out_log << "IntervalReduction time: " << t << endl;
		convex_prob->lower=split_prob->lower=reform->ext_prob->lower(0, split_prob->dim()-1);
		convex_prob->upper=split_prob->upper=reform->ext_prob->upper(0, split_prob->dim()-1);

		if (quad_prob) {
			quad_prob->lower=convex_prob->lower;
			quad_prob->upper=convex_prob->upper;
		}
		reform->ext_convex_prob->lower=reform->ext_prob->lower;
		reform->ext_convex_prob->upper=reform->ext_prob->upper;
		if (reform->ext_quad_prob) {
			reform->ext_quad_prob->lower=reform->ext_prob->lower;
			reform->ext_quad_prob->upper=reform->ext_prob->upper;
		}
	}

//	out_log << *reform->ext_convex_prob;
//	out_log << *reform->ext_quad_prob;
//	out_log << minlpdata->obj;
//	for (int c=0; c<minlpdata->con.size(); ++c)
//		out_log << "Constraint " << minlpdata->con[c].name << endl << minlpdata->con[c];

	if (out_log_p && param->get_i("Check convexification", 0)) {
		out_log << "Checking extended convexification: " << endl;
		Convexify(param).check_convexification(*reform->ext_convex_prob, *reform->ext_prob, ineq_index, 20);
	}

	if (LocOpt::nlp_solver_available()) {
		Pointer<LocOpt> locopt=LocOpt::get_solver(reform->ext_convex_prob, param, "ConvexSolve", NULL, NULL);
		int ret=locopt->solve(reform->ext_convex_prob->primal_point);

		out_out.setf(ios::fixed);
		out_out.precision(20);
		out_out << "Solving extended convex problem (Cext): return " << ret << "\t time: " << locopt->time() << "\t value: " << locopt->opt_val() << endl;
		out_out.unsetf(ios::fixed);
		out_out.precision(6);

		if ((!ret) || (!reform->ext_convex_prob->feasible(locopt->sol_point, tol, out_log_p))) // feasible solution
			sol_Cext=new dvector(locopt->sol_point);
		sol_Cext_is_solution=(ret==0);
		if (sol_Cext_is_solution && low_bound<locopt->opt_val()) low_bound=locopt->opt_val();


		locopt=NULL;		
		locopt=LocOpt::get_solver(convex_prob, param, "ConvexSolve", NULL, NULL);
		ret=locopt->solve(convex_prob->primal_point);

		out_out.setf(ios::fixed);
		out_out.precision(20);
		out_out << "Solving convex problem (C): return " << ret << "\t time: " << locopt->time() << "\t value: " << locopt->opt_val() << endl;
		out_out.unsetf(ios::fixed);
		out_out.precision(6);
	}

	linear_relax=new LinearRelax(param);
	linear_relax->init(reform->ext_convex_prob, reform->ext_prob->i_discr, sol_Cext, sol_Cext_is_solution);

	box_reduce3();
	out_out << "Reduced box diameter: " << sqrt((split_prob->upper-split_prob->lower).sq_norm2()) << endl;
	out_out << "Complete boxreduction time: " << boxreduce_time << endl;

	if (param->get_i("Level Cuts", 1)) { // Level cut
		levelcuts=new LevelCutHandler(linear_relax->get_upper_bound());
//		levelcuts->add_problem(split_prob);
		levelcuts->add_problem(convex_prob);
//		levelcuts->add_problem(reform->ext_prob);
		levelcuts->add_problem(reform->ext_convex_prob);
		levelcuts->add_problem(linear_relax);
	}

	out_out << "Preprocessing time: " << timer->stop() << endl;
}

int MinlpOpt::solve() {
	init();
	int ret;

	if (!strcmp(param->get("MinlpOpt mode", "BCP"), "BCP")) {
		out_out << "Starting Branch Cut and Price." << endl;
		ret=start_bb();
	} else if (!strcmp(param->get("MinlpOpt mode"), "off")) {
		out_out << "Preprocessing only." << endl;
		init2();
		ret=1;
	} else {
		out_err << "Error: MinlpOpt mode not known." << endl;
		exit(-1);
	}

	return ret;
}

int MinlpOpt::start_bb() {
	init2();

	Timer t;
	MinlpBCP bcp(orig_prob, split_prob, linear_relax, is_gams_prob, sol_cand_closeval_tol, sol_cand_diam, param, out_solver_p, out_solver_log_p);
	bcp.set_convex_prob(convex_prob);
	bcp.set_reform(reform, sol_Cext, sol_Cext_is_solution);
	bcp.set_levelcut_handler(levelcuts);
	bcp.set_problem_is_convex(prob_is_convex);
	bcp.set_MINLPData(minlpdata);
	bcp.set_timer(timer);
	int ret=bcp.solve(reform->ext_prob->primal_point);
	out_out << "BCP time: " << t.stop() << endl;

	if (bcp.sol_cand.size()) {
		opt_val_=bcp.sol_cand.begin()->first;
		sol_point=bcp.sol_cand.begin()->second;
	}
	sol_cand=bcp.sol_cand;

	return ret;
}


void MinlpOpt::check_initial_point() {
	bool out_of_bounds = 0;
	for (int j=0; j<orig_prob->dim(); j++)
		if (orig_prob->primal_point[j] < orig_prob->lower[j]-rtol || orig_prob->primal_point[j] > orig_prob->upper[j]+rtol ) {
			out_log << "Variable " << j << " named " << orig_prob->var_names[j]
			<< " = " << orig_prob->primal_point[j]
			<< " Out of bounds [" << orig_prob->lower[j] << "," <<  orig_prob->upper[j] << "];" << endl;
			out_of_bounds = 1;
			if (orig_prob->primal_point[j] < orig_prob->lower[j]-rtol) orig_prob->primal_point[j] = orig_prob->lower[j];
			if (orig_prob->primal_point[j] > orig_prob->upper[j]+rtol) orig_prob->primal_point[j] = orig_prob->upper[j];
    }
	//	if(out_of_bounds) exit(-1);
	// this should be uncommented for super2, ...
}
