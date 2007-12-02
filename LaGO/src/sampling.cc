// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "sampling.h"

// ----------------------------------------------- Sampling -------------------------------------------

Sampling::Sampling(Pointer<Param> param, const char* param_prefix)
: midpoint(false), bounds(false), montecarlo(5), vertices(5)
{	if (!param) return;

	int prefix_end=param_prefix ? strlen(param_prefix)+1 : 0;
	char* name=new char[prefix_end+40];
	if (param_prefix) {
		strcpy(name, param_prefix);
		name[prefix_end-1]=' ';
	}

	strcpy(name+prefix_end, "sample set mid point");
	midpoint=param->get_i(name, 0);

	strcpy(name+prefix_end, "sample set box ends");
	bounds=param->get_i(name, 0);

	strcpy(name+prefix_end, "sample set Monte Carlo");
	montecarlo=param->get_i(name, 5);

	strcpy(name+prefix_end, "sample set vertices");
	vertices=param->get_i(name, 5);

	delete[] name;
}


int Sampling::get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper) {
	int size=montecarlo;
	for (int i=0; i<montecarlo; i++) {
		sample_set.push_back(dvector(lower.dim()));
		sample_set.back().set_random(lower, upper);
	}

	if (midpoint || bounds || vertices) {
		bool bounds_missing=false;
		dvector dummy(lower);
		for (int j=0; j<lower.dim(); j++)
			if (lower(j)<=-INFINITY || upper(j)>=INFINITY) {
				if (lower(j)<=-INFINITY) dummy[j]=random(lower(j), upper(j));
//				if (lower(j)<=-INFINITY) dummy[j]=0;
				bounds_missing=true;
			}

		if (!bounds_missing) {
			if (midpoint) {
				sample_set.push_back(dvector(upper));
				sample_set.back()+=lower;
				sample_set.back()*=.5;
				size++;
			}

			if (bounds) {
				sample_set.push_back(lower);
				sample_set.push_back(upper);
				size+=2;
			}

			for (int i=0; i<vertices; i++) {
				sample_set.push_back(dvector(lower));
				for (int j=0; j<lower.dim(); j++)
					if (random(0.,1.)>=.5) sample_set.back()[j]=upper[j];
			}
			size+=vertices;
		} else { // add zero-vectors to preserve correct size
			sample_set.resize(sample_set.size()+(midpoint ? 1 : 0)+(bounds ? 2 : 0)+vertices, dummy);
		}
	}

	return size;
}

// ----------------------------------------- Sampling_Vertices -------------------------------

Sampling_Vertices::Sampling_Vertices(Pointer<Param> param, const char* param_prefix)
: vertices(0)
{	if (!param) return;
	int prefix_end=param_prefix ? strlen(param_prefix)+1 : 0;
	char* name=new char[prefix_end+40];
	if (prefix_end) {
		strcpy(name, param_prefix);
		name[prefix_end-1]=' ';
	}
	strcpy(name+prefix_end, "sample set vertices2");

	vertices=param->get_i(name, 0);
	delete[] name;
}

void Sampling_Vertices::get_points(vector<vector<dvector> >& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, const SepQcFunc& f) {
	if (!vertices) return;
	assert(sample_set.size()==f.block.size());
	int max_size=0;
	for (int k=0; k<f.block.size(); k++) {
		if (!f.s[k]) continue;
		assert(f.s[k]->sparsity_available());
		const SparsityInfo& si(((const Func*)(Func*)f.s[k])->get_sparsity());
		if (si.nonlinear->empty()) continue;

		get_points(sample_set[k], lower(f.block[k]), upper(f.block[k]), si);
		for (int l=0; l<f.block.size(); l++) {
			if (k==l) continue;
			dvector low(lower, f.block[l]);
			dvector up(upper, f.block[l]);
			while (sample_set[l].size()<sample_set[k].size()) {
				sample_set[l].push_back(dvector(f.block[l].size()));
				sample_set[l].back().set_random(low, up);
			}
		}
	}
}

void Sampling_Vertices::get_points(vector<vector<dvector> >& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, const vector<ivector>& block, const vector<Pointer<set<int> > >& i_quad, const vector<Pointer<set<int> > > &i_nonquadlin) {
	if (!vertices) return;
	assert(sample_set.size()==block.size());
	int max_size=0;
	for (int k=0; k<block.size(); k++) {
		if ((!i_quad[k]) || (!i_nonquadlin[k])) continue;

		get_points(sample_set[k], lower(block[k]), upper(block[k]), *i_quad[k], *i_nonquadlin[k]);
//		out_log << "added " << added << "; size was " << i_quad[k]->size()+i_nonquadlin[k]->size() << endl;
//		max_size=MAX(max_size, sample_set[k].size());
		for (int l=0; l<block.size(); l++) {
			if (k==l) continue;
			dvector low(lower, block[l]);
			dvector up(upper, block[l]);
			while (sample_set[l].size()<sample_set[k].size()) {
				sample_set[l].push_back(dvector(block[l].size()));
				sample_set[l].back().set_random(low, up);
			}
		}
	}

/*
	for (int k=0; k<block.size(); k++) {
		if (i_quad[k] && i_nonquadlin[k]) {
			get_points(sample_set[k], lower(block[k]), upper(block[k]), *i_quad[k], *i_nonquadlin[k]);
			max_size=MAX(max_size, sample_set[k].size());
		}
	}
	for (int k=0; k<block.size(); k++) // fill other blocks with random vectors
		if (sample_set[k].size()<max_size) {
			dvector rnd(block[k].size()); rnd.set_random(lower(block[k]), upper(block[k]));
			sample_set[k].resize(max_size, rnd);
		}
*/
}

int Sampling_Vertices::get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, set<int>& i_quad, set<int>& i_nonquadlin) {
	if (!vertices) return 0;

	dvector x(lower); x+=upper; x*=.5;
	long maxnum=0;

	set<int>::iterator it(i_quad.size() ? i_quad.begin() : i_nonquadlin.begin());

	while (it!=i_nonquadlin.end()) {
		if (lower(*it)>-INFINITY && upper(*it)<INFINITY) maxnum++;
		x[*it]=lower(*it);
		if (++it==i_quad.end()) it=i_nonquadlin.begin();
	}

	if (!maxnum) return 0; // no bounded variables regarding given index sets

	for (int i=0; i<x.dim(); i++) // checking for unboundness
		if ((lower(i)<=-INFINITY) || (upper(i)>=INFINITY))
			x[i]=random(lower(i), upper(i));

	maxnum=(long)pow(2., (double)maxnum);
	assert(maxnum>0);
	long dist=MAX(maxnum/vertices, 1);

	long switch_mask=0;
	int added=0;
	while (added<=vertices) {
		switch_mask+=dist;
		long sm=switch_mask^(switch_mask-dist);
		it=i_quad.size() ? i_quad.begin() : i_nonquadlin.begin();
		while (it!=i_nonquadlin.end()) {
			if (lower(*it)>-INFINITY && upper(*it)<INFINITY) {
				if (sm%2) x[*it]=lower[*it]+upper[*it]-x(*it);
				sm/=2;
			}
			if (++it==i_quad.end()) it=i_nonquadlin.begin();
		}

		sample_set.push_back(x);
		added++;
	}

	return added; // should be nearly the same as vertices
}

int Sampling_Vertices::get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, const SparsityInfo& si) {
	if (!vertices) return 0;

	dvector x(lower); x+=upper; x*=.5;
	long maxnum=0;

	for (VariableIterator it(si, false, true, true); it; ++it) {
		if (lower(it())>-INFINITY && upper(it())<INFINITY && lower(it())!=upper(it())) maxnum++;
		x[it()]=lower(it());
	}

	if (!maxnum) return 0; // no bounded variables regarding given index sets

	for (int i=0; i<x.dim(); i++) // checking for unboundness
		if ((lower(i)<=-INFINITY) || (upper(i)>=INFINITY))
			x[i]=random(lower(i), upper(i));

	int added=0;
	if (maxnum<=20) {
		maxnum=(long)pow(2., (double)maxnum);
		long dist=MAX(maxnum/vertices, 1);
	
		long switch_mask=0;
		while (added<=vertices) {
			switch_mask+=dist;
			long sm=switch_mask^(switch_mask-dist);
			for (VariableIterator it(si, false, true, true); it; ++it) {
				if (lower(it())>-INFINITY && upper(it())<INFINITY && lower(it())!=upper(it())) {
					if (sm%2) x[it()]=lower[it()]+upper[it()]-x(it());
					sm/=2;
				}
			}
	
			sample_set.push_back(x);
			added++;
		}
	} else while (added<vertices) {
			for (VariableIterator it(si, false, true, true); it; ++it)
				if (lower(it())>-INFINITY && upper(it())<INFINITY && lower(it())!=upper(it()))
					if (random(0.,1.)>0.5) x[it()]=upper(it());
					else x[it()]=lower(it());
			sample_set.push_back(x);
			++added;
		}		

	return added; // should be nearly the same as vertices
}


int Sampling_Vertices::get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper) {
	if (!vertices) return 0;

	dvector x(lower);

	long maxnum=0;
	for (int i=0; i<x.dim(); ++i) {
		if (lower(i)>-INFINITY && upper(i)<INFINITY && lower(i)!=upper(i)) maxnum++;
		else x[i]=random(lower(i), upper(i));
	}
	if (!maxnum) return 0; // no bounded variables regarding given index sets

	int added=0;
	if (maxnum<=20) {
		maxnum=(long)pow(2., (double)maxnum);
		
		long dist=MAX(maxnum/vertices, 1);
		long switch_mask=0;
		while (added<maxnum && added<vertices) {
			switch_mask+=dist;
			long sm=switch_mask^(switch_mask-dist);
			for (int i=0; i<x.dim(); ++i) {
				if (lower(i)>-INFINITY && upper(i)<INFINITY && lower(i)!=upper(i)) {
					if (sm%2) x[i]=lower(i)+upper(i)-x(i);
					sm/=2;
				}
			}
			sample_set.push_back(x);
			++added;
		}
	} else while (added<vertices) {
			for (int i=0; i<x.dim(); ++i)
				if (lower(i)>-INFINITY && upper(i)<INFINITY && lower(i)!=upper(i))
					if (random(0.,1.)>0.5) x[i]=upper(i);
					else x[i]=lower(i);
			sample_set.push_back(x);
			++added;
		}		
	
	out_log << 'V' << added << ' ';
	return added; // should be vertices, approximately
}

// ------------------------------------ SimpleEquationSolver -----------------------------------------

int SimpleEquationSolver::solve(dvector& start) {
	double val1=f.eval(start);
	assert(fabs(val1)>rtol);

	if (fabs(val1)<tol) {
		sol_point=start;
		opt_val_=val1;
		return 0;
	}

	dvector grad(start.dim());
	f.grad(grad, start); if (val1<0) grad*=-1.;

	dvector bestpoint_so_far(start);
	opt_val_=fabs(val1);

	int ret=run(bestpoint_so_far, val1, grad);
	iter_=0;
	while ((ret==1 || ret==3) && iter_++<iter_max) { // return's 2 and 4 means, that we were not able to change the point
		val1=f.eval(sol_point);
		if ((val1!=val1) || (opt_val_<fabs(val1)+rtol)) {
			sol_point=bestpoint_so_far;
			break; // couldn't improve point
		}

		bestpoint_so_far=sol_point;
		opt_val_=fabs(val1);
		f.grad(grad, bestpoint_so_far); if (val1<0) grad*=-1;
		ret=run(bestpoint_so_far, val1, grad); // this updates sol_point
	}

//	if (ret) out_log << "SimpleEquationSolver: couldn't make point feasible after " << iter_ << " Iterations. last return: " << ret << "\t |val|: " << opt_val_ << endl;
	return ret;
}

int SimpleEquationSolver::run(const dvector& start, double val1, dvector& dir) {
	double t=INFINITY;
	for (int i=0; i<f.dim(); i++) {
		if (dir(i)!=dir(i)) {
//			out_log << "SimpleEquationSolver: nan in direction." << endl;
			return 2;
		}
		if (fabs(dir(i))<rtol) { dir[i]=0.; continue; }

		double t_new=0.;
		t_new = ((dir(i)<0) ? upper(i) : lower(i))-start(i);
		if (fabs(t_new)<rtol) { dir[i]=0.; continue; } // close to bounds
		t_new/=-dir(i);
		if (t_new>0) t=MIN(t, t_new);
		else dir[i]=0.; // points outside of box -> project on box
	}
//	out_log << "SimpleEquationSolver initial t: " << t << endl;
	if (t==INFINITY) {
//		out_log << "SimpleEquationSolver: couldn't find initial steplength (or steplength==0)." << endl;
		return 4;
	}

	dir*=-t;

	sol_point=start; sol_point+=dir;
	double val2=f.eval(sol_point);
	if (val2!=val2) {
//		out_log << "SimpleEquationSolver: point at other extreme of the box not evaluable (val=nan)" << endl;
		sol_point=start; // restore start-point
		return 2;
	}

	if (fabs(val2)<tol) {
		opt_val_=val2;
		return 0;
	}

	if (val1*val2>0) {
//		out_log << "SimpleEquationSolver: maybe no zero between starting point and other edge, val1=" << val1 << "\t val2=" << val2 << endl;
		return 1;
	}

	dvector x1(start);

	int bisection_iter=0;
	while (bisection_iter++<bisection_iter_max) { // bisection algorithm
//		out_log << "SimpleEquationSolver iter " << iter_ << " val1: " << val1 << "\t val2: " << val2 << endl;
		dir*=.5;

		dvector x_tmp(x1); x_tmp+=dir;
		double val_tmp=f.eval(x_tmp);

		if (val_tmp!=val_tmp) {
//			out_log << "SimpleEquationSolver: nonevaluable (val=nan) point" << endl;
			return 2;
		}

		if (fabs(val_tmp)<tol) {
			sol_point=x_tmp;
			opt_val_=val_tmp;
			return 0;
		}

		if (val1*val_tmp<0) { // sign changed between x1 and x_tmp, so use this as new interval
			sol_point=x_tmp;
			val2=val_tmp;
		} else { // sign changed between x_tmp and sol_point, so use this as new interval
			x1=x_tmp;
			val1=val_tmp;
		}
	}

//	out_log << "SimpleEquationSolver: iteration limit exceeded " << val1 << "\t " << val2 << endl;
	return 3;
}

// ---------------------------------------- Sampling_check -------------------------------------------

void Sampling_check::check(vector<vector<dvector> >& sample_set, const SepQcFunc& f, bool eq, const vector<Pointer<set<int> > >& i_lin, const UserVector<double>& lower, const UserVector<double>& upper, const vector<int>& start, int min_keep) {
	assert(sample_set.size()==f.block.size());
//	if (!start.size()) start.resize(f.block.size(), 0);

	list<pair<pair<int,int>, double> > lin; // linear variables and their coefficients
	// ((block-nr, index-in-block), coefficient)

	int minsize=INF;
	int maxstart=0;

	double lin_coeff;
	for (int k=0; k<f.block.size(); k++) {
		minsize=MIN(sample_set[k].size(), minsize);
		if (start.size()) maxstart=MAX(start[k], maxstart);
		if (!minsize) return;
		if (i_lin[k] && i_lin[k]->size()) {
			SparseVector<double> grad(f.block[k].size());
			if (f.s[k] && i_lin[k]->size()) f.s[k]->grad(grad, sample_set[k][0]);
			for (set<int>::iterator it(i_lin[k]->begin()); it!=i_lin[k]->end(); it++) {
				// use only bounded variables
				if (lower(f.block[k][*it])<=-INFINITY || upper(f.block[k][*it])>=INFINITY) continue;

				lin_coeff=0;
				if (f.s[k]) lin_coeff+=grad(*it);
				if (f.b[k]) lin_coeff+=2*(*f.b[k])(*it);
				if (lin_coeff) lin.push_back(pair<pair<int,int>, double>(pair<int,int>(k, *it), lin_coeff));
			}
		}
	}

	list<pair<double, int> > bad_points; // list of bad points, we want to remove at the end, sorted by function value

	for (int i=maxstart; i<minsize; i++) {
		double val=f.c;
		for (int k=0; k<f.block.size(); k++)
			val+=f.eval(sample_set[k][i], k);

		bool accept=false;
		if (val==val && val>-INFINITY && val<INFINITY) {
			if (val<rtol && ((!eq) || val>-rtol)) accept=true; // already feasible

			// setting linear variables
			for (list<pair<pair<int,int>, double> >::iterator it(lin.begin()); (!accept) && (it!=lin.end()); it++) {
				double x_new=sample_set[it->first.first][i][it->first.second]-val/it->second;
				double low=lower[f.block[it->first.first][it->first.second]];
				double up=upper[f.block[it->first.first][it->first.second]];

				if (eq) {
					if (x_new<low-rtol) x_new=low;
					else if (x_new>up+rtol) x_new=up;
					else accept=true;
				} else {
					if (it->second>0) {
						if (x_new<low-rtol) x_new=low;
						else {
							x_new=random(low, MIN(x_new, up));
							accept=true;
						}
					} else {
						if (x_new>up+rtol) x_new=up;
						else {
							x_new=random(MAX(x_new, low), up);
							accept=true;
						}
					}
				}

				val+=it->second*(x_new-sample_set[it->first.first][i][it->first.second]); // this is the new value
				sample_set[it->first.first][i][it->first.second]=x_new;
			}

			if (!accept) {
				dvector x(f.dim());
				for (int k=0; k<f.block.size(); k++) x.set_block(sample_set[k][i], f.block[k]);

				SimpleEquationSolver ses(f, lower, upper, param);
				accept=ses.solve(x)==0;
				if (accept)
					for (int k=0; k<f.block.size(); k++) sample_set[k][i]=ses.sol_point(f.block[k]);
			}

			if (!accept) bad_points.push_back(pair<double,int>(fabs(val), i));

		} else {
//			out_log << "Sampling: sample points function value is nan, will remove point" << endl;
			bad_points.push_back(pair<double, int>(INFINITY, i));
		}
	}

	bad_points.sort();
	while (bad_points.size() && (minsize>min_keep || bad_points.back().first==INFINITY)) { // remove as much bad points as allowed or needed
		int index=bad_points.back().second;
		bad_points.pop_back();
		for (int k=0; k<f.block.size(); k++) {
			vector<dvector>::iterator it(sample_set[k].begin());
			for (int j=0; j<index; j++) it++;
			sample_set[k].erase(it);
		}
		minsize--;
		for (list<pair<double, int> >::iterator it(bad_points.begin()); it!=bad_points.end(); it++)
			if (it->second>index) it->second--;
	}

	if (bad_points.size())
		out_log << "Sampling: Kept " << bad_points.size() << " bad points, max. function value: " << bad_points.back().first << endl;
	if (minsize<min_keep)
		out_log << "Sampling: Couldn't reach minimium set size " << min_keep << ", generated only " << minsize << " points, max. function value: " << bad_points.back().first<< endl;

	if (maxstart==minsize)
		out_log << "Sampling: Warning: all points couldn't be made feasible" << endl;
}

void Sampling_check::check(vector<vector<dvector> >& sample_set, const SepQcFunc& f, bool eq, const UserVector<double>& lower, const UserVector<double>& upper, const vector<int>& start, int min_keep) {
	assert(sample_set.size()==f.block.size());
//	if (!start.size()) start.resize(f.block.size(), 0);

//	list<pair<pair<int,int>, double> > lin; // linear variables and their coefficients
	// ((block-nr, index-in-block), coefficient)

	int minsize=INF;
	int maxstart=0;

	double lin_coeff;
	for (int k=0; k<f.block.size(); k++) {
		minsize=MIN(sample_set[k].size(), minsize);
		if (start.size()) maxstart=MAX(start[k], maxstart);
		if (!minsize) return;
/*		if (i_lin[k] && i_lin[k]->size()) {
			SparseVector<double> grad(f.block[k].size());
			if (f.s[k] && i_lin[k]->size()) f.s[k]->grad(grad, sample_set[k][0]);
			for (set<int>::iterator it(i_lin[k]->begin()); it!=i_lin[k]->end(); it++) {
				// use only bounded variables
				if (lower(f.block[k][*it])<=-INFINITY || upper(f.block[k][*it])>=INFINITY) continue;

				lin_coeff=0;
				if (f.s[k]) lin_coeff+=grad(*it);
				if (f.b[k]) lin_coeff+=2*(*f.b[k])(*it);
				if (lin_coeff) lin.push_back(pair<pair<int,int>, double>(pair<int,int>(k, *it), lin_coeff));
			}
		}
*/	}

	list<pair<double, int> > bad_points; // list of bad points, we want to remove at the end, sorted by function value

	for (int i=maxstart; i<minsize; i++) {
		double val=f.c;
		for (int k=0; k<f.block.size(); k++)
			val+=f.eval(sample_set[k][i], k);

		bool accept=false;
		if (val==val && val>-INFINITY && val<INFINITY) {
			if (val<rtol && ((!eq) || val>-rtol)) accept=true; // already feasible

			for (SepQcFunc::VariableIterator it(f, true, false, false); it; ++it) {
			// setting linear variables
//			for (list<pair<pair<int,int>, double> >::iterator it(lin.begin()); (!accept) && (it!=lin.end()); it++) {
				double x_new=sample_set[it.bnr()][i][it.bindex()]-val/it.coeff_lin();
				double low=lower[it()];
				double up=upper[it()];

				if (eq) {
					if (x_new<low-rtol) x_new=low;
					else if (x_new>up+rtol) x_new=up;
					else accept=true;
				} else {
					if (it.coeff_lin()>0) {
						if (x_new<low-rtol) x_new=low;
						else {
							x_new=random(low, MIN(x_new, up));
							accept=true;
						}
					} else {
						if (x_new>up+rtol) x_new=up;
						else {
							x_new=random(MAX(x_new, low), up);
							accept=true;
						}
					}
				}

				val+=it.coeff_lin()*(x_new-sample_set[it.bnr()][i][it.bindex()]); // this is the new value
				sample_set[it.bnr()][i][it.bindex()]=x_new;
			}

			if (!accept) {
				dvector x(f.dim());
				for (int k=0; k<f.block.size(); k++) x.set_block(sample_set[k][i], f.block[k]);

				SimpleEquationSolver ses(f, lower, upper, param);
				accept=ses.solve(x)==0;
				if (accept)
					for (int k=0; k<f.block.size(); k++) sample_set[k][i]=ses.sol_point(f.block[k]);
			}

			if (!accept) bad_points.push_back(pair<double,int>(fabs(val), i));

		} else {
//			out_log << "Sampling: sample points function value is nan, will remove point" << endl;
			bad_points.push_back(pair<double, int>(INFINITY, i));
		}
	}

	bad_points.sort();
	while (bad_points.size() && (minsize>min_keep || bad_points.back().first==INFINITY)) { // remove as much bad points as allowed or needed
		int index=bad_points.back().second;
		bad_points.pop_back();
		for (int k=0; k<f.block.size(); k++) {
			vector<dvector>::iterator it(sample_set[k].begin());
			for (int j=0; j<index; j++) it++;
			sample_set[k].erase(it);
		}
		minsize--;
		for (list<pair<double, int> >::iterator it(bad_points.begin()); it!=bad_points.end(); it++)
			if (it->second>index) it->second--;
	}

	if (bad_points.size())
		out_log << "Sampling: Kept " << bad_points.size() << " bad points, max. function value: " << bad_points.back().first << endl;
	if (minsize<min_keep)
		out_log << "Sampling: Couldn't reach minimium set size " << min_keep << ", generated only " << minsize << " points, max. function value: " << bad_points.back().first<< endl;

	if (maxstart==minsize)
		out_log << "Sampling: Warning: all points couldn't be made feasible" << endl;
}

// ---------------------------------- Sampling_Minimizer ---------------------------------------------------

Sampling_Minimizer::Sampling_Minimizer(Pointer<Param> param_, const char* param_prefix)
: param(param_)
{	if (!param) return;
	int prefix_end=param_prefix ? strlen(param_prefix)+1 : 0;
	char* name=new char[prefix_end+30];
	if (prefix_end) {
		strcpy(name, param_prefix);
		name[prefix_end-1]=' ';
	}
	strcpy(name+prefix_end, "sample set minimizer");

	minimizer=param->get_i(name, 0);
	delete[] name;
}

bool Sampling_Minimizer::add_minimizer(vector<vector<dvector> >& sample_set, Pointer<SepQcFunc> f, const UserVector<double>& lower, const UserVector<double>& upper) {
	if (!minimizer) return false;

	Pointer<MinlpProblem> prob=new MinlpProblem(f, lower, upper);
	Pointer<LocOpt> solver(LocOpt::get_solver(prob, param, NULL, NULL, NULL));

	double min_val=INFINITY;
	int min_index=-1;
	for (int i=0; i<sample_set[0].size(); i++) {
		double val=0.;
		for (int k=0; k<sample_set.size(); k++) {
			assert(i<sample_set[k].size());
			val+=f->eval(sample_set[k][i], k);
		}

		if (val<min_val) {
			min_val=val;
			min_index=i;
		}
	}

	int ret;
	if (min_index>=0) {
		dvector x(f->dim());
		for (int k=0; k<sample_set.size(); k++)
			x.set_block(sample_set[k][min_index], f->block[k]);
		ret=solver->solve(x);
	} else
		ret=solver->solve();

	if (ret==0 || ret==3 || ret==9) {
		out_log << "Adding minimizer to sample set, SnOpt return = " << ret << " value= " << solver->opt_val() << endl;
		for (int k=0; k<f->block.size(); k++)
			sample_set[k].push_back(solver->sol_point(f->block[k]));
		return true;
	}

	out_log << "Cannot add minimizer. SnOpt return = " << ret << endl;
	return false;
}

bool Sampling_Minimizer::add_minimizer(vector<dvector>& sample_set, Pointer<Func> f, const UserVector<double>& lower, const UserVector<double>& upper, dvector& start) {
	if (!minimizer) return false;

	Pointer<SepQcFunc> obj(new SepQcFunc(f->dim()));
	obj->s[0]=f;
	Pointer<MinlpProblem> prob=new MinlpProblem(obj, lower, upper);
	Pointer<LocOpt> solver(LocOpt::get_solver(prob, param, NULL, NULL, NULL));

	int ret=solver->solve(start);

	if (!finite(solver->opt_val()) || solver->opt_val()>=INFINITY) ret=-1;

	out_log << 'M' << ret << ' ';
	if (ret==0 || ret==2 || ret==3 || ret==4 || ret==9) {
//		out_log << "Adding minimizer to sample set, LocOpt return = " << ret << " value= " << solver->opt_val() << endl;
		sample_set.push_back(solver->sol_point);
		return true;
	}

//	out_log << "Cannot add minimizer. LocOpt return = " << ret << endl;
	return false;
}
