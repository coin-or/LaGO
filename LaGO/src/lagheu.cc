// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "lagheu.h"

void LagHeu::make_project_LP() {
	project_LP=NULL;

	int coresize=linear_relax->core_size();
	x_con_start=coresize;
	int R_dim=linear_relax->obj->dim();
	// variables x,s,t and constraints from core (R), x-s=\hat x, s-t<=0, -s-t<=0
	MipProblem baselp(3*R_dim, coresize+3*R_dim);

//	for (int i=0; i<R_dim; ++i)
//		baselp.setColBounds(i, node->lower(i), node->upper(i));

	int c=0;
	for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->couple_con.begin()); it!=linear_relax->couple_con.end(); ++it, ++c) {
		for (int k=0; k<it->b.size(); ++k)
			if (it->b[k]) {
				SparseVector<double> row(*it->b[k]); row*=2.;
				baselp.setRow(c, row, linear_relax->obj->block[k], it->eq ? -it->c : -INFINITY, -it->c);
			}
	}
	for (int k=0; k<linear_relax->block_con.size(); ++k)
		for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->block_con[k].begin()); it!=linear_relax->block_con[k].end(); ++it, ++c) {
			SparseVector<double> row(*it->b[0]); row*=2.;
			baselp.setRow(c, row, linear_relax->obj->block[k], it->eq ? -it->c : -INFINITY, -it->c);
		}
	assert(c==coresize);

	dvector row(2); row[0]=1.; row[1]=-1.;
	ivector indices(2);
	for (int i=0; i<R_dim; ++i) { // x-s = \hat x
		indices[0]=i;
		indices[1]=R_dim+i;
		baselp.setRow(c++, row, indices, 0., 0.); // correct bounds will be set in update_project_LP
	}
	for (int i=0; i<R_dim; ++i) { // s-t <= 0
		indices[0]=R_dim+i;
		indices[1]=2*R_dim+i;
		baselp.setRow(c++, row, indices, -INFINITY, 0.);
	}
	row[0]=-1.;
	for (int i=0; i<R_dim; ++i) { // -s-t < =0
		indices[0]=R_dim+i;
		indices[1]=2*R_dim+i;
		baselp.setRow(c++, row, indices, -INFINITY, 0.);
	}

	project_LP_obj=new SparseVector<double>(linear_relax->obj->b, linear_relax->obj->block);
	project_LP_obj->resize(baselp.dim());
	baselp.setObj((Pointer<UserVector<double> >)project_LP_obj, linear_relax->obj->c);

	baselp.finish();

	project_LP=MIPSolver::get_solver(baselp, param);
}

void LagHeu::update_project_LP() {
	if (!project_LP) make_project_LP();

	assert(node);

	int R_dim=linear_relax->obj->dim();

	// set box and scaling of t-term in objective
	double scal;
	for (int i=0; i<R_dim; ++i) {
		project_LP->modify_col(i, node->lower(i), node->upper(i), MipProblem::CONTINUOUS);
		scal=node->upper[i]-node->lower[i]; if (scal<1.) scal=1.;
		project_LP->modify_obj(2*R_dim+i, .5*delta/scal); // t
	}

	// add cuts for this node
	project_LP->delete_rows(project_LP_cuts); project_LP_cuts.clear();
	list<CutPool::CutInfo> cutinfos;
	linear_relax->get_cuts(cutinfos, node);
	for (list<CutPool::CutInfo>::iterator it(cutinfos.begin()); it!=cutinfos.end(); ++it) {
		if (it->type!=CutPool::SIMPLE && it->type!=CutPool::LINEARIZATION) continue;
		const SimpleCut& cut(it->type==CutPool::SIMPLE ? it->it_simplecuts->get_cut() : it->it_linearizationcuts->get_cut());
		Pointer<UserVector<double> > b(cut.coeff->getcopy()); *b*=2.;
		project_LP_cuts.push_back(project_LP->add_row(*b, linear_relax->obj->block[it->block_nr], -INFINITY, -cut.constant));
	}

	// scaling of c^Tx-term in objective
	double max_cxprod=0.; double prod;
	for (int i=0; i<R_dim; ++i) {
		prod=fabs((*project_LP_obj)(i))*MAX(fabs(node->lower(i)), fabs(node->upper(i)));
		if (prod>max_cxprod) max_cxprod=prod;
	}
	for (int i=0; i<R_dim; ++i)
		if ((*project_LP_obj)(i)) project_LP->modify_obj(i, (*project_LP_obj)(i)/max_cxprod);
}

void LagHeu::update_project_LP_x(const dvector& x) {
	// set \hat x in x-\hat x = s constraints
	for (int i=0; i<x.dim(); ++i) project_LP->modify_row(x_con_start+i, x(i), x(i));
}

MIPSolver::SolutionStatus LagHeu::solve_project_LP(dvector& sol) {
	MIPSolver::SolutionStatus ret=project_LP->solve();
	if (ret==MIPSolver::SOLVED || ret==MIPSolver::FEASIBLE) project_LP->get_primal(sol);
	return ret;
}

MIPSolver::SolutionStatus LagHeu::solve_project_LP(double& val) {
	MIPSolver::SolutionStatus ret=project_LP->solve();
	if (ret==MIPSolver::SOLVED || ret==MIPSolver::FEASIBLE) val=project_LP->get_optval();
	return ret;
}

void LagHeu::project_LP_fix_variable(int i, double val) {
	project_LP->modify_col(i, val, val, MipProblem::CONTINUOUS);
	// maybe fix also s and t ?
//	int R_dim=linear_relax->obj->dim();
}

void LagHeu::project_LP_unfix_variable(int i) {
	project_LP->modify_col(i, node->lower(i), node->upper(i), MipProblem::CONTINUOUS);
}


//-\min && e^Tt\\
//-\st &&AWz-b=s\\
//-&& s\leq t\\
//-&& -t\leq s\\
//-&&e^Tz_{I_k}=1 & k=1,\dots,p\\
//-&& z_{I_k}\geq 0 & k=1,\dots,p
void LagHeu::make_lagheu_LP(double delta) {
	lagheu_LP=NULL;
	lagheu_LP_columns.clear();

	int ccsize=linear_relax->couple_con.size();
	int blocknr=node->i_ExtremePoints.size();
	// coupling constraints from (R), s-t,s+t constraints, convex-combination constraints
	// add dummy variable in the case of no coupling constraints (-> no s,t)
	MipProblem baselp(MAX(2*ccsize,1), 3*ccsize+blocknr);

	int c=0;
	// coupling constraints
	dvector row(1); row[0]=-1.;
	ivector indices(1);
	for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->couple_con.begin()); it!=linear_relax->couple_con.end(); ++it, ++c) {
		indices[0]=c;
		baselp.setRow(c, row, indices, it->eq ? -it->c : -INFINITY, -it->c);
	}

	// s-t, s+t constraints
	row.resize(2); row[0]=1.; row[1]=-1.;
	indices.resize(2);
	while (c<2*ccsize) { // s-t <= 0
		indices[0]=c-ccsize;
		indices[1]=c;
		baselp.setRow(c++, row, indices, -INFINITY, 0.);
	}
	row[0]=-1.;
	while (c<3*ccsize) { // -s-t <= 0
		indices[0]=c-2*ccsize;
		indices[1]=c-ccsize;
		baselp.setRow(c++, row, indices, -INFINITY, 0.);
	}

	// convex-combination constraints
	for (int k=0; k<blocknr; ++k, ++c)
		baselp.setRowBounds(c, 1., 1.);

	// objective: delta * sum of t's
	Pointer<UserVector<double> > obj=new dvector(baselp.dim());
	for (int i=ccsize; i<2*ccsize; ++i) (*obj)[i]=delta;
	baselp.setObj(obj, linear_relax->obj->c);

	baselp.finish();
	z_start=baselp.dim();

	lagheu_LP=MIPSolver::get_solver(baselp, param);

	lagheu_LP_columns.resize(blocknr);

	// compute and add columns
	for (int k=0; k<blocknr; ++k) {
		vector<Pointer<UserVector<double> > > columns;
		vector<double> objcoeff;
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it) {
			Pointer<SparseVector<double> > col(new SparseVector<double>((**it).rmpcolumn));
			col->resize(ccsize); col->resize(lagheu_LP->nr_row()); // keep only coupling constraints
			col->SetElement(3*ccsize+k, 1.); // convex combination constraint
			columns.push_back(Pointer<UserVector<double> >(col));
			objcoeff.push_back((**it).rmpobjcoeff);
		}
		dvector low(columns.size()), up(columns.size(), 1.); // lower and upper bounds for z
		lagheu_LP_columns[k].clear();
		lagheu_LP->add_cols(lagheu_LP_columns[k], columns, objcoeff, low, up);
	}
}

void LagHeu::lagheu_LP_fixblock(int k, const MIPSolver::ColItem* colitem) {
	for (list<const MIPSolver::ColItem*>::iterator it(lagheu_LP_columns[k].begin()); it!=lagheu_LP_columns[k].end(); ++it) {
		if (*it==colitem) lagheu_LP->modify_col(*colitem, 1., 1., MipProblem::BINARY);
		else lagheu_LP->modify_col(**it, 0., 0., MipProblem::BINARY);
	}
}

void LagHeu::lagheu_LP_unfixblock(int k) {
	for (list<const MIPSolver::ColItem*>::iterator it(lagheu_LP_columns[k].begin()); it!=lagheu_LP_columns[k].end(); ++it)
		lagheu_LP->modify_col(**it, 0., 1., MipProblem::BINARY);
}

MIPSolver::SolutionStatus LagHeu::solve_lagheu_LP() {
	MIPSolver::SolutionStatus ret=lagheu_LP->solve();
//	out_solver_log << "Solving lagheu LP returns " << ret << "\t value: " << lagheu_LP->get_optval();
	return ret;
}

MIPSolver::SolutionStatus LagHeu::solve_lagheu_LP(dvector &z, dvector& x) {
	MIPSolver::SolutionStatus ret=solve_lagheu_LP();

	z.resize(lagheu_LP->nr_col()-z_start);
	x=0;
	for (int k=0; k<node->i_ExtremePoints.size(); ++k) {
		list<const MIPSolver::ColItem*>::iterator it_col(lagheu_LP_columns[k].begin());
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it, ++it_col) {
			int index=lagheu_LP->get_colindex(**it_col)-z_start;
			z[index]=lagheu_LP->get_primal(**it_col);
			if (z(index))
				for (int j=linear_relax->obj->block[k].size()-1; j>=0; --j)
					x[linear_relax->obj->block[k][j]]+=z(index)*(**it)(j);
		}
	}

//	out_log << "\t conviol: " << lagheu_LP->get_optval()-linear_relax->obj->eval(x) << endl;

	return ret;
}

MIPSolver::SolutionStatus LagHeu::solve_lagheu_LP(vector<dvector> &z, dvector& x) {
	MIPSolver::SolutionStatus ret=solve_lagheu_LP();

	z.resize(node->i_ExtremePoints.size());
	x=0;
	for (int k=0; k<node->i_ExtremePoints.size(); ++k) {
		list<const MIPSolver::ColItem*>::iterator it_col(lagheu_LP_columns[k].begin());
		z[k].resize(node->i_ExtremePoints[k].size());
		int i=0;
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it, ++it_col, ++i) {
			z[k][i]=lagheu_LP->get_primal(**it_col);
			if (z[k](i))
				for (int j=linear_relax->obj->block[k].size()-1; j>=0; --j)
					x[linear_relax->obj->block[k][j]]+=z[k](i)*(**it)(j);
		}
	}

//	out_log << "\t conviol: " << lagheu_LP->get_optval()-linear_relax->obj->eval(x) << endl;

	return ret;
}

double LagHeu::couple_con_violation(const dvector& x) {
	double viol=0.;
	for (list<LinearRelax::LinConstraint>::iterator it_con(linear_relax->couple_con.begin()); it_con!=linear_relax->couple_con.end(); it_con++) {
		double valc=it_con->c;
		for (int k_=0; k_<it_con->b.size(); k_++)
			if (it_con->b[k_]) valc+=2*(*it_con->b[k_] * x(linear_relax->obj->block[k_]));
		if (valc>0 || it_con->eq) viol+=fabs(valc);
	}
	return viol;
}


bool LagHeu::project_and_round(dvector& x) {
	update_project_LP_x(x);
	MIPSolver::SolutionStatus ret=solve_project_LP(x);
	if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) return false;

	// collect binary variables
	map<int, double> unfixed;
	for (int i=0; i<orig_prob->i_discr.size(); ++i) {
		int i0=orig_prob->i_discr[i];
		if (node->lower(i0)==node->upper(i0)) continue;

		double key=MIN(fabs(node->lower(i0)-x(i0)), fabs(node->upper(i0)-x(i0)))/(node->upper(i0)-node->lower(i0));
		unfixed.insert(pair<int,double>(i0, key));
	}
	if (unfixed.empty()) return true; // continuous problem or all binaries are already fixed

	int switch_count=0;
	bool ret2=project_and_round_rek(x, unfixed, switch_count);
	out_log << (ret2 ? '.' : (switch_count==switch_count_limit ? 's' : '!'));// << endl;

	return ret2;
}

bool LagHeu::project_and_round_rek(dvector& x, map<int, double>& unfixed, int& switch_count) {
	if (unfixed.empty()) return true;

	// fix a variable with high key
	map<int,double>::iterator index_it(unfixed.begin());
	for (map<int,double>::iterator it(++unfixed.begin()); it!=unfixed.end(); ++it)
		if (it->second<index_it->second) index_it=it;
	int index=index_it->first;
	double key=index_it->second;
	unfixed.erase(index_it);

	double rounded=2.*x(index)>node->lower(index)+node->upper(index) ? node->upper(index) : node->lower(index);
	project_LP_fix_variable(index, rounded); // fix variable

//	out_log << 'f' << index << ' ';

	bool needswitch=false;

	if (key>.1) { // if variable was very unsure, resolve project_LP
		MIPSolver::SolutionStatus ret=solve_project_LP(x);
		if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) needswitch=true;
		else { // update keys and call recursion
			for (map<int,double>::iterator it(unfixed.begin()); it!=unfixed.end(); ++it)
				it->second=MIN(fabs(node->lower(it->first)-x(it->first)), fabs(node->upper(it->first)-x(it->first)))/(node->upper(it->first)-node->lower(it->first));
			needswitch=!project_and_round_rek(x, unfixed, switch_count);
		}
	} else {
		if (unfixed.empty()) { // need to be sure, that fixation is feasible
			if (key==0) needswitch=false; // variable was implicitly fixed by LP before
			else {
				MIPSolver::SolutionStatus ret=solve_project_LP(x);
				if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) needswitch=true;
				else needswitch=false;
			}
		} else needswitch=!project_and_round_rek(x, unfixed, switch_count);
	}

	bool win=true;

	if (needswitch && switch_count<switch_count_limit) {
		++switch_count;
		project_LP_fix_variable(index, node->lower(index)+node->upper(index)-rounded); // switch variable
//		out_log << 's' << index << ' ';

		MIPSolver::SolutionStatus ret=solve_project_LP(x);
		if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) win=false;
		else { // update keys and call recursion
			for (map<int,double>::iterator it(unfixed.begin()); it!=unfixed.end(); ++it)
				it->second=MIN(fabs(node->lower(it->first)-x(it->first)), fabs(node->upper(it->first)-x(it->first)))/(node->upper(it->first)-node->lower(it->first));
			win=project_and_round_rek(x, unfixed, switch_count);
		}
	} else if (needswitch) win=false; // switch count limit exceeded

	project_LP_unfix_variable(index);
	unfixed.insert(pair<int,double>(index,key));
//	out_log << 'u' << index << ' ';
	return win;
}

// -------------------------------------------- LagHeu1 ------------------------------------------------------

// find solution candidates via Lagrange heuristic
int LagHeu1::solve(Pointer<MinlpNode> node_) {
	node=node_;
	update_project_LP();

	int max_sol=20;
	dvector x(linear_relax->obj->dim());
	dvector x2(x.dim()), x_star(x.dim());

	// compute a solution z of lagheu-LP
	make_lagheu_LP(delta);
	dvector z;
	MIPSolver::SolutionStatus ret=solve_lagheu_LP(z, x_star);
	if (ret!=MIPSolver::SOLVED) {
		out_solver_log << "Couldn't solve Lagheu LP: " << ret << endl;
		return 1;
	}
// int ret; x_star=node->ref_point;
	vector<vector<pair<double,list<ExtremePoint>::iterator> > > sorted_RMP_points(node->i_ExtremePoints.size());
	for (int k=0; k<node->i_ExtremePoints.size(); ++k) {
		multimap<double, list<ExtremePoint>::iterator> RMPpts;
		x2=x_star;
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it) {
			// compute key for the RMP-point
			x2.set_block(**it, linear_relax->obj->block[k]);
			double val=linear_relax->obj->eval(x2)+delta*couple_con_violation(x2);
			RMPpts.insert(pair<double, list<ExtremePoint>::iterator>(val, *it));
//			out_log << k << "\t " << val << "\t " << (*it)->first;
		}
		sorted_RMP_points[k].reserve(RMPpts.size());
		for (multimap<double, list<ExtremePoint>::iterator>::iterator it(RMPpts.begin()); it!=RMPpts.end(); it++) {
//			out_log << k << ":\t " << it->first << "\t " << it->second->first;
			sorted_RMP_points[k].push_back(*it);
		}
		x.set_block(*RMPpts.begin()->second, linear_relax->obj->block[k]);
	}
	out_log << "RMP points sorted." << endl;

	int locoptret=1;
	if (project_and_round(x)) {
		locoptret=locopt_NLP(x).first;
		out_log << locoptret << ' ';
	}

// $\lb p_k=\min_j p_{k,j}$, $\ub p_k=\max_j p_{k,j}$
// $s_{k,j}=\sum_{i\leq j} 1/((p_{k,i}-\lb p_k)/((\ub p_k-\lb p_k)+0.1)$

	long poss_combinations=1;
	for (int k=0; k<sorted_RMP_points.size(); k++) {
		double pmax=-INFINITY;
		double pmin=INFINITY;
		for(int i=0; i<sorted_RMP_points[k].size(); i++) {
			if(sorted_RMP_points[k][i].first>pmax) pmax=sorted_RMP_points[k][i].first;
			if(sorted_RMP_points[k][i].first<pmin) pmin=sorted_RMP_points[k][i].first;
		}
		double sum=0;
		for(int i=0;i<sorted_RMP_points[k].size(); i++) {
			sum+=1/((sorted_RMP_points[k][i].first-pmin)/(fabs(pmax-pmin)+1)+0.1);
			sorted_RMP_points[k][i].first=sum;
//			out_log << k << ":\t " << 1/((sorted_RMP_points[k][i].first-pmin)/(fabs(pmax-pmin)+1)+0.1) << "\t " << sum << endl;
		}
//		if (sorted_RMP_points[k].size()>9) out_log << " "; out_log << sorted_RMP_points[k].size(); if (sorted_RMP_points[k].size()>9) out_log << " ";
		if (poss_combinations<max_sol) poss_combinations*=sorted_RMP_points[k].size();
	}
//	out_log << endl;
	if (max_sol>poss_combinations) max_sol=poss_combinations;

	// compute solution candidates
	while (max_sol--) {
		for (int k=0; k<sorted_RMP_points.size(); k++) {
			//int l=(int)genexp(sorted_RMP_points[k].size()/2.);
			//if (l>=sorted_RMP_points[k].size()) l=sorted_RMP_points[k].size()-1;
			int l=0;
			if (sorted_RMP_points[k].size()>1) {
				double s=random(0., sorted_RMP_points[k].back().first); // compute random number
//				double s=genexp(sorted_RMP_points[k].back().first/4.); // compute random number
				for(int i=0;i<sorted_RMP_points[k].size(); i++) {
					l=i;
					if (sorted_RMP_points[k][i].first>s) break;
				}
			}
//			if (l>9) out_log << " "; out_log << l; if (l>9) out_log << " ";
//			out_log << k << ":\t " << sorted_RMP_points[k].back().first/2 << "\t " << s << "\t " << l << endl;
			x.set_block(*sorted_RMP_points[k][l].second, linear_relax->obj->block[k]);
		}

//		out_log << endl;
		if (project_and_round(x)) {
			int ret=locopt_NLP(x).first;
			out_log << ret;
			if (locoptret) locoptret=ret;
		}
	}
	out_log << endl;

	return locoptret; // best local optimization return, we obtained
}

// ---------------------------------------- LagHeu_SimAnnealing ------------------------------------

LagHeu_SimAnnealing::LagHeu_SimAnnealing(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob, double closeval_tol_, Pointer<dvector> diam_, Pointer<Param> param_, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: LagHeu(orig_prob_, linrelax_, is_gams_prob, closeval_tol_, diam_, param_, out_solver_p_, out_solver_log_p_),
	Z(linear_relax->obj->block.size())
{ iter_max=param->get_i("Simulated Annealing iter max", MAX(20, Z.size()));
	minor_iter_max=linear_relax->obj->dim()/Z.size(); // average block size
	if (!strcmp(param->get("Simulated Annealing weights", "violation"), "violation")) weight_type=VIOLATION;
	else weight_type=DISTANCE;
}

void LagHeu_SimAnnealing::get_Wz(dvector& Wz, const ivector& z) {
	for (int k=0; k<Z.size(); ++k)
		Wz.set_block(*Z[k][z(k)], linear_relax->obj->block[k]);
}

double LagHeu_SimAnnealing::weight(const ivector& z, const dvector& Wz) {
	switch (weight_type) {
		case VIOLATION: {
			double val=linear_relax->obj->c;
			for (int k=0; k<z.dim(); ++k) val+=Z[k][z(k)]->rmpobjcoeff;
			int c=0;
			for (list<LinearRelax::LinConstraint>::iterator it_con(linear_relax->couple_con.begin()); it_con!=linear_relax->couple_con.end(); ++it_con, ++c) {
				double valc=it_con->c;
				for (int k=0; k<it_con->b.size(); k++) valc+=Z[k][z(k)]->rmpcolumn[c];
				if (valc>0 || it_con->eq) val+=delta*fabs(valc);
			}
			return val;
		} break;
		case DISTANCE: {
			update_project_LP_x(Wz);
			double val;
			MIPSolver::SolutionStatus ret=solve_project_LP(val);
			if (ret!=MIPSolver::SOLVED) out_solver << "Solving Project LP return: " << ret << endl;
			return val;
		} break;
	}
	out_err << "Weight type now known. Aborting." << endl; exit(-1);
	return 0;
}

void LagHeu_SimAnnealing::walk(ivector& z, dvector& Wz, double& weight_z, double T, int max_iter) {
	ivector z2(z);
	dvector Wz2(Wz);
	double weight_z2;

	out_solver_log << "Walk for temperature T=" << T << "\t iteration " << iter() << " / " << iter_max << '\t';

	for (int it=0; it<max_iter; ++it) {
//		out_log << "z: " << z;
 		// pick a block, which contains more than one RMP point
 		int k=nonsingles[random(0, nonsingles.size()-1)];
 		// pick a RMP point, which is different from the current
 		int i=random(0, Z[k].size()-2);
		if (i>=z(k)) ++i;

		z2[k]=i;
		Wz2.set_block(*Z[k][i], linear_relax->obj->block[k]);
		weight_z2=weight(z2, Wz2);

		double p=exp((weight_z-weight_z2)/T);
//		out_solver_log << "w(z) = " << weight_z << "\t w(z') = " << weight_z2 << "\t prob: " << p;
		double r=1.;
		if (p>1 || (r=random(0.,1.))<p) { // accept
			out_solver_log << 'a';
//			out_solver_log << "\t accepted (" << r << ")." << endl;
			z[k]=z2[k];
			Wz.set_block(*Z[k][i], linear_relax->obj->block[k]);
			weight_z=weight_z2;
		}// else out_solver_log << "\t declined (" << r << ")." << endl;
			else out_solver_log << 'd';
	}
	out_solver_log << '\t';
}

int LagHeu_SimAnnealing::solve(Pointer<MinlpNode> node_) {
	node=node_;
	update_project_LP();

	// compute a solution z* of lagheu-LP
	dvector z_star;
	dvector x(linear_relax->obj->dim());

	make_lagheu_LP(delta);
	MIPSolver::SolutionStatus ret=solve_lagheu_LP(z_star, x);
	if (ret!=MIPSolver::SOLVED) {
		out_solver_log << "Couldn't solve Lagheu LP: " << ret << endl;
		return 1;
	}

	ivector z(Z.size()); // starting point

	// copy RMP points to vector-vector structure and compute z0
	int i=0;
	long poss=1;
	nonsingles.clear();
	for (int k=0; k<node->i_ExtremePoints.size(); ++k) {
		Z[k].clear();
		double maxz=0.; int j=0;
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it, ++i, ++j) {
			Z[k].push_back(*it);
			if (maxz<z_star(i)) { maxz=z_star(i); z[k]=j; }
		}
		if (Z[k].size()>1) nonsingles.push_back(k);
		poss*=Z[k].size();
	}

	out_solver_log << "Possible combinations: " << poss << endl;

	if (nonsingles.empty()) {
		out_solver << "Only one RMP point per block. Not many possibilites to walk." << endl;
		if (project_and_round(x)) {
			int locoptret=locopt_NLP(x).first;
			out_solver_log << "LocOpt return: " << locoptret << endl;
			return locoptret;
		}
	}

	out_solver_log << "Nonsingles: " << nonsingles.size() << " \t out of " << Z.size() << endl;

	dvector Wz(x.dim());
	get_Wz(Wz, z);
	double start_weight=weight(z, Wz);
	double T=fabs(10*start_weight);
	double rho=pow(10.,-3./iter_max);
	out_solver_log << "rho: " << rho << endl;

	int locoptret=1;
	double weight_z=start_weight;
	for (iter_=0; iter()<iter_max; ++iter_) {
		walk(z, Wz, weight_z, T, minor_iter_max);

		x=Wz;
		if (project_and_round(x)) {
			int ret=locopt_NLP(x).first;
			out_solver_log << "LocOpt return: " << ret << endl;
			if (locoptret) locoptret=ret;
		}

//		T=fabs(start_weight)*((1./100.-10.)*iter()/(double)iter_max+10.); // linear decreasing
		T*=rho;
	}

	return locoptret;
}

// ------------------------------------- LagHeu2 ------------------------------------------

LagHeu2::LagHeu2(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob, double closeval_tol_, Pointer<dvector> diam_, Pointer<Param> param_, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: LagHeu(orig_prob_, linrelax_, is_gams_prob, closeval_tol_, diam_, param_, out_solver_p_, out_solver_log_p_),
	x(linrelax_->obj->dim()), penalty(new MinlpPenaltyFunc(orig_prob_))
{ max_candidates=param->get_i("LagHeu2 max locopt", 100);
  max_ns=param->get_i("LagHeu2 max combinations", 10000);
	penalty->delta=1.;
}

int LagHeu2::solve(Pointer<MinlpNode> node_) {
	node=node_;
	update_project_LP();

	dvector z_star;

	make_lagheu_LP(delta);
	MIPSolver::SolutionStatus lagheulpret=solve_lagheu_LP(z_star, x);
	if (lagheulpret!=MIPSolver::SOLVED) {
		out_solver_log << "Couldn't solve Lagheu LP: " << lagheulpret << endl;
		return 1;
	}

	// sort ExtremePoints; compute ns
	W.clear();
	candidates.clear();

	double ns=1;
	int i=0;
	for (int k=0; k<node->i_ExtremePoints.size(); ++k) {
		pair<pair<double, int>, multimap<double, list<ExtremePoint>::iterator> > sortedblock;
		sortedblock.first.first=0.;
		sortedblock.first.second=k;

		int size=0;
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it, ++i) {
			if (z_star(i)<1E-4) continue;
			++size;
			sortedblock.second.insert(pair<double, list<ExtremePoint>::iterator>(z_star(i), *it));
			sortedblock.first.first+=z_star(i)*z_star(i);
		}
		sortedblock.first.first-=1./(double)size;
		sortedblock.first.first/=(double)size;

		W.insert(sortedblock);
		if (ns<max_ns) ns*=size;

//		out_log << k << ":" << sortedblock.first.first << ',' << size << ' ';
	}
	if (ns>max_ns) ns=max_ns;
	ns_all=ns;
//	out_log << endl;
	ns_done=0;
	lastprint=0;

	out_solver_log << "Start searching for " << ns << " combinations. Percentage done: ";
	search(W.rbegin(), ns);
	out_solver_log << endl;

	if (candidates.empty()) {
		out_solver << "LagHeu2: No candidates for local optimization found." << endl;
		return 1;
	}

	int ret=1;
	out_solver_log << "Start local minimization for " << candidates.size() << " points: " << endl;
	for (map<double, dvector>::iterator it(candidates.begin()); it!=candidates.end(); ++it) {
//		Round::round(it->second, it->second, split_prob->i_discr, node->lower, node->upper);
		int locret=locopt_NLP(it->second).first;
		out_solver_log << locret;
		if (ret) ret=locret;
	}
	out_solver_log << endl;

	return ret;
}

void LagHeu2::search(map<pair<double, int>, multimap<double, list<ExtremePoint>::iterator> >::reverse_iterator it_k, double ns) {
	if (it_k==W.rend() || ns<1) {
		pair<double, dvector> cand(0, x);
		if (project_and_round(cand.second)) {
			cand.first=key(cand.second);
			candidates.insert(cand);
			if (candidates.size()>max_candidates) candidates.erase(--candidates.end());
		}
		ns_done+=ns;
		while (100.*ns_done/ns_all>=lastprint+10) {
			lastprint+=10;
			out_solver_log << lastprint << ' ';
		}
		//out_solver_log << ns_done << ' ';
		return;
	}

//out_log << it_k->first.second << ',' << ns << ' ';
	map<pair<double, int>, multimap<double, list<ExtremePoint>::iterator> >::reverse_iterator next(it_k); ++next;

	dvector xk(x, linear_relax->obj->block[it_k->first.second]);

	double restns=ns;
	for (multimap<double, list<ExtremePoint>::iterator>::reverse_iterator it(it_k->second.rbegin()); it!=it_k->second.rend(); ++it) {
		x.set_block(*it->second, linear_relax->obj->block[it_k->first.second]);
		search(next, it->first*ns);
		restns-=it->first*ns;
	}
	x.set_block(xk, linear_relax->obj->block[it_k->first.second]);
}

void LagHeu2::set_reform(Pointer<Reformulation> reform_) {
	LagHeu::set_reform(reform_);
	penalty=new MinlpPenaltyFunc(reform->ext_prob);
	penalty->delta=1.;
}

double LagHeu2::key(const dvector& x) {
	double val=0;
	if (reform) val=reform->ext_prob->obj->eval(x);
	else val=orig_prob->obj->eval(x);
	val+=penalty->eval(x);
	return val;
}

bool LagHeu2::project(dvector& res) {
	update_project_LP_x(x);
	MIPSolver::SolutionStatus ret=solve_project_LP(res);
	if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) {
		out_solver_log << "project LP return: " << ret << endl;
		return false;
	}
	return true;
}

// ------------------------------------- LagHeu2b ------------------------------------------

LagHeu2b::LagHeu2b(Pointer<MinlpProblem> orig_prob_, Pointer<LinearRelax> linrelax_, bool is_gams_prob, double closeval_tol_, Pointer<dvector> diam_, Pointer<Param> param_, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: LagHeu(orig_prob_, linrelax_, is_gams_prob, closeval_tol_, diam_, param_, out_solver_p_, out_solver_log_p_),
	x(linrelax_->obj->dim()), penalty(new MinlpPenaltyFunc(orig_prob_)), W_val(linrelax_->obj->block.size())
{ max_candidates=param->get_i("LagHeu2 max locopt", 100);
  max_ns=param->get_i("LagHeu2 max combinations", 10000);
	penalty->delta=1.;
}

int LagHeu2b::block_value(double& val, int k) {
	val=0.;
	int size=0;
	for (int i=0; i<z[k].dim(); ++i) {
		if (z[k](i)<1E-4) continue;
		++size;
		val+=z[k](i)*z[k](i);
	}
	val-=1./(double)size;
	val/=(double)size;

	return size;
}

int LagHeu2b::solve(Pointer<MinlpNode> node_) {
	node=node_;
	update_project_LP();

	make_lagheu_LP(delta);
	MIPSolver::SolutionStatus lagheulpret=solve_lagheu_LP(z, x);
	if (lagheulpret!=MIPSolver::SOLVED) {
		out_solver_log << "Couldn't solve Lagheu LP: " << lagheulpret << endl;
		return 1;
	}

	double ns=1;
	unfixed_blocks.clear();
	for (int k=0; k<W_val.size(); ++k) {
		int size=block_value(W_val[k], k);
		if (ns<max_ns) ns*=size;
		unfixed_blocks.push_back(k);
	}
	if (ns>max_ns) ns=max_ns;
	ns_all=ns;
	ns_done=0;
	lastprint=0;
	candidates.clear();

	out_solver_log << "Start searching for " << ns << " combinations. Percentage done: ";
	search(ns);
	out_solver_log << endl;

	if (candidates.empty()) {
		out_solver << "LagHeu2b: No candidates for local optimization found." << endl;
		return 1;
	}

	int ret=1;
	out_solver_log << "Start local minimization for " << candidates.size() << " points: " << endl;
	for (multimap<double, dvector>::iterator it(candidates.begin()); it!=candidates.end(); ++it) {
//		Round::round(it->second, it->second, split_prob->i_discr, node->lower, node->upper);
		int locret=locopt_NLP(it->second).first;
		out_solver_log << locret;
		if (ret) ret=locret;
	}
	out_solver_log << endl;

	return ret;
}

void LagHeu2b::search(double ns) {
	if (unfixed_blocks.empty() || ns<1) {
		if (project_and_round(x)) {
			candidates.insert(pair<double, dvector>(key(), x));
			if (candidates.size()>max_candidates) candidates.erase(--candidates.end());
		}
		ns_done+=ns;
		while (100.*ns_done/ns_all>=lastprint+10) {
			lastprint+=10;
			out_solver_log << lastprint << ' ';
		}
//		out_solver_log << ns_done << ' ';
//out_solver_log << '.';
		return;
	}

	// get block with highest value
	list<int>::iterator myk_it=unfixed_blocks.begin();
	for (list<int>::iterator it(++unfixed_blocks.begin()); it!=unfixed_blocks.end(); ++it)
		if (W_val[*it]>W_val[*myk_it]) myk_it=it;
	int myk=*myk_it;
	unfixed_blocks.erase(myk_it);
//out_log << 'b' << myk << ": " << z[myk];
	// sort z-values for this block
	multimap<double, pair<list<ExtremePoint>::iterator, const MIPSolver::ColItem*> > z_sort;
	int i=0;
	list<const MIPSolver::ColItem*>::iterator it_colitem(lagheu_LP_columns[myk].begin());
	for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[myk].begin()); it!=node->i_ExtremePoints[myk].end(); ++it, ++it_colitem, ++i) {
		if (z[myk](i)<1E-4) continue;
		z_sort.insert(pair<double, pair<list<ExtremePoint>::iterator, const MIPSolver::ColItem*> >(z[myk](i), pair<list<ExtremePoint>::iterator, const MIPSolver::ColItem*>(*it, *it_colitem)));
	}

	dvector x0(x);

	double restns=ns;
	for (multimap<double, pair<list<ExtremePoint>::iterator, const MIPSolver::ColItem*> >::reverse_iterator it(z_sort.rbegin()); it!=z_sort.rend() && it->first>1E-4; ++it) {
		lagheu_LP_fixblock(myk, it->second.second);

		if (it->first<.8) { // resolve lagheu_LP
			solve_lagheu_LP(z, x);
			MIPSolver::SolutionStatus lagheulpret=solve_lagheu_LP(z, x);
			if (lagheulpret!=MIPSolver::SOLVED) {
				out_solver_log << "Couldn't solve Lagheu LP: " << lagheulpret << endl;
				ns_done+=it->first*ns;
				continue;
			}
			for (list<int>::iterator it_k(unfixed_blocks.begin()); it_k!=unfixed_blocks.end(); ++it_k)
				block_value(W_val[*it_k], *it_k);
		} else {
			x.set_block(*it->second.first, linear_relax->obj->block[myk]);
		}

		search(it->first*ns);
		restns-=it->first*ns;
	}

	lagheu_LP_unfixblock(myk);
	unfixed_blocks.push_back(myk);
	x=x0;
}

void LagHeu2b::set_reform(Pointer<Reformulation> reform_) {
	LagHeu::set_reform(reform_);
	penalty=new MinlpPenaltyFunc(reform->ext_prob);
	penalty->delta=1.;
}

double LagHeu2b::key() {
	double val=0;
	if (reform) val=reform->ext_prob->obj->eval(x);
	else val=orig_prob->obj->eval(x);
	val+=penalty->eval(x);
	return val;
}

bool LagHeu2b::project() {
	update_project_LP_x(x);
	MIPSolver::SolutionStatus ret=solve_project_LP(x);
	if (ret!=MIPSolver::SOLVED && ret!=MIPSolver::FEASIBLE) {
		out_solver_log << "project LP return: " << ret << endl;
		return false;
	}
	if (split_prob->i_discr.size()) {
		Round::round(x,x, split_prob->i_discr, node->lower, node->upper);
		make_binaries_feasible();
	}
	return true;
}

bool LagHeu2b::switch_binaries(double val, bool eq, map<int,double>& coeff) {
	if (eq) out_log << 'e' << coeff.size();
	double oldval=val;
	for (int iter=0; iter<coeff.size(); ++iter) {
		vector<int> poss; poss.reserve(coeff.size()); // possiblities of improvement
		for (map<int,double>::iterator it(coeff.begin()); it!=coeff.end(); ++it) {
			if (val*it->second>0 && x(it->first)==split_prob->upper(it->first)) poss.push_back(it->first);
			if (val*it->second<0 && x(it->first)==split_prob->lower(it->first)) poss.push_back(it->first);
		}
		if (poss.empty()) {
			out_log << 'p';
			return false;
		}

		int index=poss[random(0, poss.size()-1)];
		double newval=val+coeff[index]*(split_prob->lower(index)+split_prob->upper(index)-2*x(index));
		if (fabs(val)>(eq ? fabs(newval) : newval)) { // closer to feasiblity -> accept
			x[index]=split_prob->lower(index)+split_prob->upper(index)-x(index);
			val=newval;
		}

		if (val<1E-4 && (!eq || val>-1E-4)) {
			out_log << '.';
			return true;
		}
	}
	out_log << '+';
	for (map<int,double>::iterator it(coeff.begin()); it!=coeff.end(); ++it) { // try all
		if (!(val*it->second>0 && x(it->first)==split_prob->upper(it->first)) && !(val*it->second<0 && x(it->first)==split_prob->lower(it->first))) continue;
		double newval=val+coeff[it->first]*(split_prob->lower(it->first)+split_prob->upper(it->first)-2*x(it->first));
		if (fabs(val)>(eq ? fabs(newval) : newval)) { // closer to feasiblity -> accept
			x[it->first]=split_prob->lower(it->first)+split_prob->upper(it->first)-x(it->first);
			val=newval;
		}
		if (val<1E-4 && (!eq || val>-1E-4)) {
			out_log << '.';
			return true;
		}
	}
	out_log << '!';// << oldval << ' ' << val << '\t';
	return false;
}

void LagHeu2b::switch_binaries2_rek(double val, bool eq, double& bestval, dvector& pt, map<int,double>::iterator coeff_it, map<int,double>& coeff) {
	if (coeff_it==coeff.end()) {
		if ((eq && fabs(val)<fabs(bestval)) || (!eq && bestval>1E-4 && val<bestval)) {
			bestval=val;
			x=pt;
		}
		return;
	}

	map<int,double>::iterator next(coeff_it); next++;

	switch_binaries2_rek(val, eq, bestval, pt, next, coeff);
	if (bestval<1E-4 && (!eq || bestval>-1E-4)) return;

	pt[coeff_it->first]=split_prob->lower(coeff_it->first)+split_prob->upper(coeff_it->first)-pt(coeff_it->first);
	double val2=val+coeff_it->second*(split_prob->lower(coeff_it->first)+split_prob->upper(coeff_it->first)-2*x(coeff_it->first));

	switch_binaries2_rek(val2, eq, bestval, pt, next, coeff);

	pt[coeff_it->first]=split_prob->lower(coeff_it->first)+split_prob->upper(coeff_it->first)-pt(coeff_it->first);
}

void LagHeu2b::switch_binaries2(double val, bool eq, map<int,double>& coeff) {
	if (coeff.size()>8) { switch_binaries(val, eq, coeff); return; }

	dvector y(x);

	double best=val;
	switch_binaries2_rek(val, eq, best, y, coeff.begin(), coeff);

//	if (val>1E-4 || (eq && val<-1E-4)) { out_log << '!'; }
//	else out_log << '.';
//	out_log << val << "->" << best << '\t';
}

void LagHeu2b::move_continuous(double& val, bool eq, map<int,double>& cont_coeff) {
	for (map<int,double>::iterator it(cont_coeff.begin()); it!=cont_coeff.end(); ++it) {
		double move=-val/it->second;
		if (x(it->first)+move < node->lower(it->first)) move=node->lower(it->first)-x(it->first);
		if (x(it->first)+move > node->upper(it->first)) move=node->upper(it->first)-x(it->first);

		x[it->first]+=move;
		val+=it->second*move;

		if (node->lower(it->first)>x[it->first]+rtol || node->upper(it->first)<x[it->first]-rtol)
			out_err << "Warning: Moved outside box: " << x[it->first] << " not in " << node->lower(it->first) << ", " << node->upper(it->first) << endl;

		if (val<1E-4 && (!eq || val>-1E-4)) {
//			out_log << 'm';
			break;
		}
	}
}

void LagHeu2b::make_binaries_feasible() {
	for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->couple_con.begin()); it!=linear_relax->couple_con.end(); ++it) {
		map<int, double> coeff; // collect coefficients of binary variables
		map<int, double> cont_coeff; // collect coefficients of continuous variables
		double val=it->c;
		for (int k=0; k<it->b.size(); ++k)
			if (it->b[k]) {
				val+=2.* (*it->b[k]*x(linear_relax->obj->block[k]));
				for (int i=0; i<it->b[k]->size(); ++i) {
					int i0=linear_relax->obj->block[k][i];
					double co=2.*(*it->b[k])(i);
					if (!co) continue;
					if (node->lower(i0)==node->upper(i0)) continue;
					if (i0<split_prob->dim() && split_prob->discr[i0]) coeff.insert(pair<int, double>(i0,co));
					else cont_coeff.insert(pair<int, double>(i0,co));
				}
			}
		if ((val>1E-4 || (it->eq && val<-1E-4)) && !coeff.empty()) {
			if (!cont_coeff.empty()) move_continuous(val, it->eq, cont_coeff);
			if (val>1E-4 || (it->eq && val<-1E-4)) switch_binaries2(val, it->eq, coeff);
		}
	}
	for (int k=0; k<linear_relax->block_con.size(); ++k) {
		dvector xk(x, linear_relax->obj->block[k]);
		for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->block_con[k].begin()); it!=linear_relax->block_con[k].end(); ++it) {
			double val=it->c + 2.*(*it->b[0] * xk);
			if (val<1E-4 && ((!it->eq) || val>-1E-4)) continue; // satisfied
			map<int, double> coeff; // collect coefficients of binary variables
			map<int, double> cont_coeff; // collect coefficients of continuous variables
			for (int i=0; i<xk.dim(); ++i) {
				int i0=linear_relax->obj->block[k][i];
				double co=2.*(*it->b[0])(i);
				if (!co) continue;
				if (node->lower(i0)==node->upper(i0)) continue;
				if (i0<split_prob->dim() && split_prob->discr[i0]) coeff.insert(pair<int, double>(i0,co));
				else cont_coeff.insert(pair<int, double>(i0,co));
			}
			if (!cont_coeff.empty()) move_continuous(val, it->eq, cont_coeff);
			if (val>1E-4 || (it->eq && val<-1E-4)) switch_binaries2(val, it->eq, coeff);
		}
	}
//	if (linear_relax->point_feasible(NULL, x)) { out_log << 'f' << endl; }
//	else out_log << 'i' << endl;
}
