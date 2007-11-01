// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "linrelax.h"
#include "node.h"

// ----------------------------- LinearRelax -----------------------------------

LinearRelax::LinearRelax(Pointer<Param> param_)
: param(param_), max_cutsnr(500000), solver_node(NULL)
{ inactivetime_limit_global=param->get_i("Cut inactive time limit global", 10);
	inactivetime_limit_local=param->get_i("Cut inactive time limit local", 3);
	max_cutsnr=param->get_i("max cuts nr", 500000);
	levelcut_pos=couple_con.end();
}

LinearRelax::LinearRelax(const LinearRelax& linrelax, int k, Pointer<MinlpNode> node, Pointer<SepQcFunc> alt_obj)
: param(linrelax.param), max_cutsnr(linrelax.max_cutsnr), inactivetime_limit_global(linrelax.inactivetime_limit_global), inactivetime_limit_local(linrelax.inactivetime_limit_local), solver_node(NULL),
	lower(node->lower, linrelax.obj->block[k]), upper(node->upper, linrelax.obj->block[k]),
	block_con(1, linrelax.block_con[k]), cutpool(1, new CutPool(*linrelax.cutpool[k], node)), cutpoolcoupling(new CutPool(linrelax.inactivetime_limit_global, linrelax.inactivetime_limit_local))
{ if (alt_obj) obj=alt_obj;
	else {
		obj=new SepQcFunc(linrelax.obj->block[k].size());
		obj->b[0]=linrelax.obj->b[k];
		obj->c=linrelax.obj->c;
	}
	levelcut_pos=couple_con.end();
	solver=new LinearRelaxSolverMIP(*this);
}

int LinearRelax::core_size() const {
	int size=couple_con.size();
	for (int k=0; k<block_con.size(); k++) size+=block_con[k].size();
	return size;
}

void LinearRelax::add_cut(const Pointer<SimpleCut>& cut, int k, const Pointer<MinlpNode>& node) {
	assert(k<(int)cutpool.size());
	CutPool::CutInfo cutinfo((k>=0 ? cutpool[k] : cutpoolcoupling)->add_cut(cut, node, k)); // add the cut to the cutpool
	if (solver && ((!node) || solver_node==node)) solver->add_cut(cutinfo); // if it's the correct node, notify the solver
}

void LinearRelax::add_cut(const Pointer<LinearizationCut>& cut, int k, const Pointer<MinlpNode>& node) {
	assert(k<(int)cutpool.size());
	CutPool::CutInfo cutinfo((k>=0 ? cutpool[k] : cutpoolcoupling)->add_cut(cut, node, k)); // add the cut to the cutpool
	if (solver && ((!node) || solver_node==node)) solver->add_cut(cutinfo); // if it's the correct node, notify the solver
}

void LinearRelax::add_cut(Pointer<IntervalGradientCut> intgradcut, int k, Pointer<MinlpNode> node) {
	assert(k>=0 && k<cutpool.size());
	if (!intgradcut) return;
	CutPool::CutInfo cutinfo(cutpool[k]->add_cut(intgradcut, node, k)); // add the cut to the cutpool
	if (solver && ((!node) || solver_node==node)) solver->add_cut(cutinfo); // if it's the correct node, notify the solver
}

void LinearRelax::integrate(LinearRelax& linrelax, int k, Pointer<MinlpNode> node) {
	assert(k<cutpool.size());
	cutpool[k]->integrate(*linrelax.cutpool[0], node);
}

void LinearRelax::get_cuts(list<CutPool::CutInfo>& cutinfos, Pointer<MinlpNode> node) {
	for (int k=0; k<cutpool.size(); ++k) cutpool[k]->get_cuts(cutinfos, k, node);
	cutpoolcoupling->get_cuts(cutinfos, -1, node);
}

void LinearRelax::update_cuts(Pointer<MinlpNode> node, int k, IntervalGradientCutGenerator& generator, LinearizedConCutGenerator& linconcutgen) {
	cutpool[k]->update_cuts(node, k, node->lower(obj->block[k]), node->upper(obj->block[k]), generator, linconcutgen);
	// no IntervalGradientCuts in cutpoolcoupling
}

void LinearRelax::init(Pointer<MinlpProblem> convex_prob, const vector<int>& i_discr_, Pointer<dvector> feas_point, bool is_feas) {
	out_log << "Constructing linear relaxation" << endl;
	block_con.resize(convex_prob->block.size());
	cutpool.resize(convex_prob->block.size());
	for (int k=0; k<convex_prob->block.size(); k++) cutpool[k]=new CutPool(inactivetime_limit_global, inactivetime_limit_local);
	cutpoolcoupling=new CutPool(inactivetime_limit_global, inactivetime_limit_local);

	obj=convex_prob->obj;
	lower=convex_prob->lower;
	upper=convex_prob->upper;
	i_discr=i_discr_;

	LinearizedConCutGenerator linconcut(convex_prob);

	// sort linear constraints and linearize active convexified constraints, if feasible point of (C) is known
	for (int c=0; c<convex_prob->con.size(); c++) { // adding linear constraints
		bool linear=true; int block_con_nr=-1;
		for (int k=0; k<convex_prob->block.size(); k++) {
			linear&=!(convex_prob->con[c]->A[k] || convex_prob->con[c]->s[k]);
			if (convex_prob->con[c]->A[k] || convex_prob->con[c]->b[k] || convex_prob->con[c]->s[k]) {
				if (block_con_nr<0) block_con_nr=k; // first block, we found
				else block_con_nr=INF; // found a second block
			}
		}
		if (block_con_nr<0) { // constant constraint !
			if (convex_prob->con[c]->c<1E-4 && ((!convex_prob->con_eq[c]) || convex_prob->con[c]->c>-1E-4)) continue;
			out_out << "Warning: Infeasible constant constraint found!" << endl;
			block_con_nr=0;
		}
		if (linear) {
			if (block_con_nr<INF) block_con[block_con_nr].push_back(LinConstraint(convex_prob->con[c]->b[block_con_nr], convex_prob->con[c]->c, convex_prob->con_eq[c], convex_prob->con_names[c]));
			else couple_con.push_back(LinConstraint(convex_prob->con[c]->b, convex_prob->con[c]->c, convex_prob->con_eq[c], convex_prob->con_names[c]));
		} else if (feas_point) {
			double val=convex_prob->con[c]->eval(*feas_point);
			if (is_feas && val>1E-4) out_log << "Constraint " << convex_prob->con_names[c] << " violated: " << val << endl;
			if ((val>-1E-4 || !is_feas) && block_con_nr<INF) {
				Pointer<LinearizationCut> cut(new LinearizationCut(linconcut.get_cut((*feas_point)(convex_prob->block[block_con_nr]), c, block_con_nr, val)));
				if (cut->coeff) { 
//					out_log << "adding cut on con. " << convex_prob->con_names[c] << " block " << block_con_nr << " const.= " << cut->constant << "\t coeff.=" << *cut->coeff;
					add_cut(cut, block_con_nr);
				}
			}
		}
	}

//	solver=new LinearRelaxSolverGeneral(*this);
	solver=new LinearRelaxSolverMIP(*this);

	if (feas_point && is_feas) { // checking feasible point from (C) for feasiblity for (R)
		bool feas=point_feasible(NULL, *feas_point);
		if (!feas) out_out << "Warning: Solution of (Cext) not feasible for linear relaxation." << endl;
	}
	
//	out_log << nr_all_cuts() << " cuts. " << endl;
}

// computing level cut: using v_up = upper bound of initial linear relaxation (R)
double LinearRelax::get_upper_bound() {
	solver->construct();
	solver->set_objective(new SepQcFunc(*obj, true));

	int ret=solver->solve();
	double val=-solver->opt_val();
	out_out << "Computation of upper bound of linear relaxation. return: " << ret << "\t value: " << val << endl;

	solver->set_objective(obj);
	return (!ret) ? val : INFINITY;
}

void LinearRelax::add_level_cut(double level) {
	int block_nr=-1;
	for (int k=0; k<obj->block.size(); k++)
		if (obj->b[k])
			if (block_nr<0) block_nr=k;
			else block_nr=INF;
	if (block_nr<0) {
		block_nr=0;
		out_out << "Warning: Constant objective function. Not adding level cut." << endl;
		return;
	}
	if (block_nr<INF) {
		block_con[block_nr].push_back(LinConstraint(obj->b[block_nr], obj->c-level, false, strdup("level cut")));
		levelcut_pos=--block_con[block_nr].end();
	} else {
		couple_con.push_front(LinConstraint(obj->b, obj->c-level, false, strdup("level cut")));
		levelcut_pos=couple_con.begin();
	}

	if (solver) solver->init(); // reinit solver, because core linear relaxation changed
}

void LinearRelax::update_level_cut(double newlevel) {
	if (levelcut_pos!=couple_con.end()) levelcut_pos->c=obj->c-newlevel;
	if (solver) solver->update_levelcut(newlevel);
}

int LinearRelax::nr_local_cuts(Pointer<MinlpNode> node) const {
	int sum=cutpoolcoupling->nr_local_cuts(node);
	for (int k=0; k<cutpool.size(); k++) sum+=cutpool[k]->nr_local_cuts(node);
	return sum;
}

int LinearRelax::nr_global_cuts() const {
	int sum=cutpoolcoupling->nr_global_cuts();
	for (int k=0; k<cutpool.size(); k++) sum+=cutpool[k]->nr_global_cuts();
	return sum;
}

int LinearRelax::nr_all_cuts() const {
	int sum=cutpoolcoupling->nr_all_cuts();
	for (int k=0; k<cutpool.size(); k++) sum+=cutpool[k]->nr_all_cuts();
	return sum;
}

void LinearRelax::remove_node(Pointer<MinlpNode> node) {
	for (int k=0; k<cutpool.size(); k++) cutpool[k]->remove_node(node);
	cutpoolcoupling->remove_node(node);
	if (solver_node==node) solver->reset();
}

void LinearRelax::duplicate_nodeinfo(Pointer<MinlpNode> oldnode, Pointer<MinlpNode> newnode) {
	for (int k=0; k<cutpool.size(); k++) cutpool[k]->duplicate_nodeinfo(oldnode, newnode);
	cutpoolcoupling->duplicate_nodeinfo(oldnode, newnode);
}

void LinearRelax::clear_solver() {
	solver->reset();
	solver_node=NULL;
}

int LinearRelax::solve(dvector& sol_point, double& value, Pointer<MinlpNode> node, dvector* dual_point, double tol) {
	if ((!solver_node) || (solver_node!=node)) { solver->construct(node); solver_node=node; }
//	else solver->update(node);

	double oldtol=solver->tol;
	solver->tol=tol;
	
	int ret=solver->solve(node);
//	if (ret==0) { generate_cuts(node); ret=solver->solve(node); }

	solver->tol=oldtol;
	

	value=solver->opt_val();
	sol_point=solver->sol_point;
	if (dual_point) *dual_point=solver->duals;

	for (int k=0; k<cutpool.size(); k++) cutpool[k]->update_nodeinfo(solver->cutinfos, k, node);
	cutpoolcoupling->update_nodeinfo(solver->cutinfos, -1, node);
	solver->remove_cuts();

	return ret;
}

bool LinearRelax::feasible(Pointer<MinlpNode> node, double tol) {
	if ((!solver_node) || (solver_node!=node)) solver->construct(node);

	double oldtol=solver->tol;
	solver->tol=tol;
	bool feas=solver->feasible();
	solver->tol=oldtol;
	return feas;
}

bool LinearRelax::point_feasible(Pointer<MinlpNode> node, const dvector& x, double tol) const {
	vector<dvector> xb; xb.reserve(obj->block.size());
	for (int k=0; k<obj->block.size(); k++) xb.push_back(dvector(x, obj->block[k]));

	for (list<LinConstraint>::const_iterator it(couple_con.begin()); it!=couple_con.end(); it++) {
		double val=it->c;
		for (int k=0; k<it->b.size(); k++)
			if (it->b[k]) val+=2*(*it->b[k] * xb[k]);
		if (val>tol || (it->eq && val<-tol)) return false;
	}
	for (int k=0; k<block_con.size(); k++)
		for (list<LinConstraint>::const_iterator it(block_con[k].begin()); it!=block_con[k].end(); it++) {
			double val=it->c+2*(*it->b[0] * xb[k]);
			if (val>tol || (it->eq && val<-tol)) return false;
		}

	for (int k=0; k<cutpool.size(); k++) // checks cuts
		if (!cutpool[k]->feasible(node, xb[k], tol)) return false;
	if (!cutpoolcoupling->feasible(node, x, tol)) return false;

	return true;
}

bool LinearRelax::box_reduce(pair<double,double>& newbox, const pair<double,double>& oldbox, int k, int i, bool discrete, bool unknown_only) {
	int ret;
	newbox=oldbox;

//	if (oldup-oldlow<rtol) return newbox; // fixed variable

	Pointer<SepQcFunc> tempobj(new SepQcFunc(obj->block));
	tempobj->b[k]=new SparseVector<double>(obj->block[k].size(), i, .5); // for lower bound

	if ((!unknown_only) || oldbox.first<=-INFINITY) {
		solver->set_objective(tempobj);
		ret=solver->solve();
		if (!ret) newbox.first=project(solver->opt_val(), oldbox.first, oldbox.second);
		else if (ret==1) {
			solver->set_objective(obj);
			return false; // relaxation infeasible
		}
//		else out_log << "LinearRelax::box_reduce lower of " << obj->block[k][i] << " returned nonzero." << endl;
		if (discrete)
			if (newbox.first>oldbox.first+1E-4) {
//				out_log << oldbox.first << " -> " << newbox.first;
				newbox.first=MIN(upperint(newbox.first), oldbox.second);
//				out_log << " -> " << newbox.first << endl;
				boxreduce_reduced_integer.insert(pair<int,int>(k, i));
				if (newbox.first==oldbox.second) {
					solver->set_objective(obj);
					return true;
				}
			} else newbox.first=oldbox.first;  // avoid minimal (<1E-4) perturbation in bound of discrete variable
	}

	if ((!unknown_only) || oldbox.second>=INFINITY) {
		*tempobj->b[k]*=-1.;
		solver->set_objective(tempobj);
		ret=solver->solve();
		if (!ret) newbox.second=project(-solver->opt_val(), newbox.first, oldbox.second);
		else if (ret==1) {
			solver->set_objective(obj);
			return false;
		}
//		else out_log << "LinearRelax::box_reduce upper of " << obj->block[k][i] << " returned nonzero." << endl;
		if (discrete)
			if (newbox.second<oldbox.second-1E-4) {
//				out_log << oldbox.second << " -> " << newbox.second;
				newbox.second=MAX(lowerint(newbox.second), oldbox.first);
//				out_log << " -> " << newbox.second << endl;
				boxreduce_reduced_integer.insert(pair<int,int>(k, i));
				if (newbox.second==oldbox.first) {
					solver->set_objective(obj);
					return true;
				}
			}	else newbox.second=oldbox.second; // avoid minimal (<1E-4) perturbation in bound of discrete variable
	}

	solver->set_objective(obj);

	return true;
}

double LinearRelax::box_reduce(Pointer<MinlpNode> node, int k, const vector<bool>& discr, bool unknown_only, set<pair<int, IntervalReduction::which_bound_type> >* changed_var, double min_impr) {
	if ((!solver_node) || (solver_node!=node)) solver->construct(node);

	pair<double,double> newbox;
	double reduction=0.;
	int i0;

	dvector& low(node ? node->lower : lower);
	dvector& up(node ? node->upper : upper);

	for (int i=obj->block[k].size()-1; i>=0; --i) {
		i0=obj->block[k][i];
		if ((unknown_only && low(i0)>-INFINITY && up(i0)<INFINITY) || up(i0)-low(i0)<rtol) {
			reduction+=1.;
			continue;
		}

		if (!box_reduce(newbox, pair<double,double>(low(i0), up(i0)), k, i, discr[i0], unknown_only)) {
//			out_log << "Reducing box of variable " << i0 << " returned infeasible." << endl;
			return -INFINITY;
		}

//		out_log << "reduce box of " << i0 << " from " << interval<double>(low(i0), up(i0)) << " to " << interval<double>(newbox.first, newbox.second) << endl;

		if (low(i0)<=-INFINITY || up(i0)>=INFINITY)
			if (newbox.first<=-INFINITY || newbox.second>=INFINITY) reduction+=1.; else ;
		else reduction+=up(i0)-low(i0) ? (newbox.second-newbox.first)/(up(i0)-low(i0)) : 1.;

		if (changed_var) { // if we want specific information about the variables, which bound was changed
			if (newbox.first-low(i0)>min_impr*(up(i0)-low(i0))) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i0, IntervalReduction::LOWER));
			if (up(i0)-newbox.second>min_impr*(up(i0)-low(i0))) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i0, IntervalReduction::UPPER));
		}

		low[i0]=newbox.first;
		up[i0]=newbox.second;
	}

	return reduction/obj->block[k].size();
}

double LinearRelax::box_reduce(Pointer<MinlpNode> node, const vector<bool>& discr, bool unknown_only, set<pair<int, IntervalReduction::which_bound_type> >* changed_var, double min_impr) {
	double reduction=0.;
	double redblock;
	for (int k=0; k<obj->block.size(); ++k) {
		redblock=box_reduce(node, k, discr, unknown_only, changed_var, min_impr)*obj->block[k].size();
		if (redblock==-INFINITY) return -INFINITY;
		reduction+=redblock;
	}
	return reduction/obj->dim();
}

double LinearRelax::box_reduce(dvector& newlow, dvector& newup, const vector<bool>& discr, bool unknown_only, set<pair<int, IntervalReduction::which_bound_type> >* changed_var, double min_impr) {
	double reduction=box_reduce(NULL, discr, unknown_only, changed_var, min_impr);
	newlow=lower;
	newup=upper;
	return reduction;
}

int LinearRelax::generate_cuts(Pointer<MinlpNode> node) {
	list<Pointer<SimpleCut> > cuts;
	int nr=solver->generate_cuts(cuts);
	for (list<Pointer<SimpleCut> >::iterator it(cuts.begin()); it!=cuts.end(); ++it)
		add_cut(*it, -1, node);
		
	return nr;
}

ostream& operator<<(ostream& out, LinearRelax& linrelax) {
		out << "Objective: " << *linrelax.obj;
		out << "Coupling constraints: " << endl;
		for (list<LinearRelax::LinConstraint>::const_iterator it(linrelax.couple_con.begin()); it!=linrelax.couple_con.end(); ++it)
			out << *it;
		for (int k=0; k<linrelax.block_con.size(); ++k) {
			if (!linrelax.block_con[k].empty()) {
				out << "Constraints block " << k << ':' << endl;
				for (list<LinearRelax::LinConstraint>::const_iterator it(linrelax.block_con[k].begin()); it!=linrelax.block_con[k].end(); ++it)
					out << *it;								
			}			
		}

		list<CutPool::CutInfo> cutinfos;
		linrelax.get_cuts(cutinfos);
		out << "Cuts: " << endl;
		for (list<CutPool::CutInfo>::iterator it(cutinfos.begin()); it!=cutinfos.end(); ++it) {
			switch (it->type) {
				case CutPool::SIMPLE: out << it->it_simplecuts->get_cut(); break;
				case CutPool::LINEARIZATION: out << it->it_linearizationcuts->get_cut(); break;
				case CutPool::INTERVALGRADIENT: out << "IntervalgradientCut skipped" << endl; break;
			}
		}		
		return out;
}

ostream& operator<<(ostream& out, const LinearRelax::LinConstraint& lincon) {
	out << lincon.name << ':';
	for (int k=0; k<lincon.b.size(); ++k) {
		if (lincon.b[k])
			for (int i=0; i<lincon.b[k]->dim(); ++i)
				if ((*lincon.b[k])(i)) out << " +" << 2*(*lincon.b[k])(i) << "*x" << i;
	}
	out << " + " << lincon.c << (lincon.eq ? " =0" : " <=0") << endl;	
	return out;	
}


// ---------------------------------------- LinearRelaxSolver ----------------------------------------

LinearRelaxSolver::LinearRelaxSolver(LinearRelax& linrelax_)
: Solver(linrelax_.obj->dim()), linrelax(linrelax_)
{ }

LinearRelaxSolver::~LinearRelaxSolver() { }

int LinearRelaxSolver::solve() { return solve(NULL); }

// ------------------------------------- LinearRelaxSolverGeneral -------------------------------------

void LinearRelaxSolverGeneral::init() {
	dim_=linrelax.obj->dim();
	sol_point.resize(linrelax.obj->dim());

	prob=new MinlpProblem(linrelax.obj, linrelax.lower, linrelax.upper, false);
	for (list<LinearRelax::LinConstraint>::const_iterator it(linrelax.couple_con.begin()); it!=linrelax.couple_con.end(); it++) {
		prob->add_con(new SepQcFunc(prob->block), it->eq, it->name);
		prob->con.back()->b=it->b;
		prob->con.back()->c=it->c;
		if (it==linrelax.levelcut_pos) levelcut_pos=prob->con.size()-1;
	}
	for (int k=0; k<linrelax.block_con.size(); k++)
		for (list<LinearRelax::LinConstraint>::const_iterator it(linrelax.block_con[k].begin()); it!=linrelax.block_con[k].end(); it++) {
			prob->add_con(new SepQcFunc(prob->block), it->eq, it->name);
			prob->con.back()->b[k]=it->b[0];
			prob->con.back()->c=it->c;
			if (it==linrelax.levelcut_pos) levelcut_pos=prob->con.size()-1;
		}
	cutstart=prob->con.size();
}

void LinearRelaxSolverGeneral::reset() {
	assert(prob->con.size()>=cutstart);
	prob->con.resize(cutstart);
	prob->con_eq.resize(cutstart);
	prob->con_names.resize(cutstart);

	prob->obj=linrelax.obj;
	prob->lower=linrelax.lower;
	prob->upper=linrelax.upper;

	cutinfos.clear();
}

void LinearRelaxSolverGeneral::construct(Pointer<MinlpNode> node) {
	reset();
	linrelax.get_cuts(cutinfos, node);
	for (list<CutPool::CutInfo>::iterator it(cutinfos.begin()); it!=cutinfos.end(); it++) {
		assert(it->type==CutPool::LINEARIZATION || it->type==CutPool::SIMPLE);
		assert(it->block_nr>=0);
		const SimpleCut& cut(it->type==CutPool::LINEARIZATION ? it->it_linearizationcuts->get_cut() : it->it_simplecuts->get_cut());
		Pointer<SepQcFunc> con(new SepQcFunc(prob->block));
		con->b[it->block_nr]=cut.coeff;
		con->c=cut.constant;
		prob->add_con(con, false);
	}

	if (node) { // set bounds
		prob->lower=node->lower;
		prob->upper=node->upper;
	}
}

void LinearRelaxSolverGeneral::add_cut(const CutPool::CutInfo& cutinfo) {
	assert(cutinfo.type==CutPool::LINEARIZATION || cutinfo.type==CutPool::SIMPLE);
	assert(cutinfo.block_nr>=0);
	const SimpleCut& cut(cutinfo.type==CutPool::LINEARIZATION ? cutinfo.it_linearizationcuts->get_cut() : cutinfo.it_simplecuts->get_cut());
	cutinfos.push_back(cutinfo);
	Pointer<SepQcFunc> con(new SepQcFunc(prob->block));
	con->b[cutinfo.block_nr]=cut.coeff;
	con->c=cut.constant;
	prob->add_con(con, false);
}

void LinearRelaxSolverGeneral::set_objective(Pointer<SepQcFunc> obj) {
	prob->obj=obj;
}

void LinearRelaxSolverGeneral::update_levelcut(double newlevel) {
	if (levelcut_pos==-1) return;
	prob->con[levelcut_pos]->c=linrelax.obj->c-newlevel;
}

int LinearRelaxSolverGeneral::solve(Pointer<MinlpNode> node) {
	Pointer<LocOpt> solver(LocOpt::get_lp_solver(prob, linrelax.param, NULL, NULL, NULL));
	solver->tol=tol;
	int ret=solver->solve();
	opt_val_=solver->opt_val();
	sol_point=solver->sol_point;
	dvector allduals(solver->get_lag_multipliers());
	duals=allduals(0, cutstart-1);
	solver=NULL;//	if (ret>1) return 3;
	if (ret && prob->feasible(sol_point, tol, NULL)) return 1; // not feasible
	if (ret) return 3; // solver thinks infeasible, put point is feasible for us.

	if (linrelax.inactivetime_limit_global || linrelax.inactivetime_limit_local) {
		// update inactivity info
		list<CutPool::CutInfo>::iterator it(cutinfos.begin());
		for (int c=cutstart; c<allduals.size(); c++, it++) it->inactive=fabs(allduals(c))<rtol;
	}

	return 0; // solved
}

void LinearRelaxSolverGeneral::remove_cuts() {
	int removed=0;
	int c=cutstart;
	list<CutPool::CutInfo>::iterator it(cutinfos.begin());
	while (it!=cutinfos.end()) {
		if (it->removed) {
			prob->del_con(c); // very slow!!
			removed++;
			it=cutinfos.erase(it);
		} else {
			c++;
			it++;
		}
	}
//		if (removed) out_log << "Deleted " << removed << " cuts." << endl;
}

bool LinearRelaxSolverGeneral::feasible() {
	Pointer<SepQcFunc> obj=prob->obj;
	prob->obj=new SepQcFunc(prob->block);
	Pointer<LocOpt> solver(LocOpt::get_lp_solver(prob, linrelax.param, NULL, NULL, NULL));
	solver->tol=tol;
	int ret=solver->solve();
	sol_point=solver->sol_point;
	solver=NULL;
	prob->obj=obj;
	if (!ret) return true;
	return !prob->feasible(sol_point, tol, NULL);
}

// -------------------------------- LinearRelaxSolverMIP -----------------------------

LinearRelaxSolverMIP::LinearRelaxSolverMIP(LinearRelax& linrelax_)
: LinearRelaxSolver(linrelax_), mip(linrelax_.lower, linrelax_.upper, linrelax_.i_discr), levelcut_pos(-1)
{	tol=1E-4;
	init();
}


void LinearRelaxSolverMIP::init() {
	dim_=linrelax.obj->dim();
	sol_point.resize(linrelax.obj->dim());

	mip.reset();

	int coresize=linrelax.core_size();
	mip.resize_row(MAX(coresize,1)); // to avoid, that OSI-CPLEX doesn't create an LP if there are no rows

	int c=0;
	for (list<LinearRelax::LinConstraint>::iterator it(linrelax.couple_con.begin()); it!=linrelax.couple_con.end(); it++, c++) {
		SparseVector<double> b(it->b, linrelax.obj->block); b*=2;
		mip.setRow(c, b, it->eq ? -it->c : -INFINITY, -it->c, it->name);
		if (it==linrelax.levelcut_pos) levelcut_pos=c;
	}
	for (int k=0; k<linrelax.block_con.size(); k++)
		for (list<LinearRelax::LinConstraint>::iterator it(linrelax.block_con[k].begin()); it!=linrelax.block_con[k].end(); it++, c++) {
			Pointer<UserVector<double> > b(it->b[0]->getcopy()); *b*=2;
			mip.setRow(c, *b, linrelax.obj->block[k], it->eq ? -it->c : -INFINITY, -it->c, it->name);
			if (it==linrelax.levelcut_pos) levelcut_pos=c;
		}

	Pointer<UserVector<double> > b(new SparseVector<double>(linrelax.obj->b, linrelax.obj->block)); *b*=2;
	mip.setObj(b, linrelax.obj->c);

	mip.finish();

	duals.resize(mip.rows());
	sol_point.resize(mip.dim());

	solver=MIPSolver::get_solver(mip, linrelax.param);
	intervalgradcutvars_start=solver->nr_col();
}

void LinearRelaxSolverMIP::reset() {
	solver->reset();
	for (int i=0; i<mip.dim(); i++) solver->modify_col(i, linrelax.lower(i), linrelax.upper(i), mip.getColType(i));

	cutinfos.clear();
}

void LinearRelaxSolverMIP::add_intervalgradientcut(CutPool::CutInfo& cutinfo) {
	const IntervalGradientCut& cut(cutinfo.it_intgradcuts->get_cut());

	if (cut.coninfos.empty()) return;

	int size=cut.indices.dim();
	int bsize=cut.ref_x.dim();

	// add variables w and z
	int w_z_begin=solver->nr_col();
	dvector w_z_low(2*size); // zero
	solver->add_cols(cutinfo.colitems, w_z_low, cut.w_z_up);
//out_log << "w z up: " << cut.w_z_up;

	vector<pair<dvector, ivector> > rows(cut.coninfos.size()+2*size);
	dvector rowlb(rows.size());
	dvector rowub(rows.size());
	int row_index=0;
	
	// main constraints
	dvector b(bsize+2*size);
	ivector ind(b.dim());
	for (list<IntervalGradientCut::ConstraintInfo>::const_iterator it(cut.coninfos.begin()); it!=cut.coninfos.end(); ++it) {
		for (int i=0; i<bsize; ++i) {
			b[i]=(*it->b_x)(i);
			ind[i]=linrelax.obj->block[cutinfo.block_nr](i);
			if (i<size) {
				b[i+bsize]=(*it->b_w)(i);
				b[i+bsize+size]=(*it->b_z)(i);
				ind[i+bsize]=w_z_begin+i;
				ind[i+bsize+size]=w_z_begin+size+i;
			}
		}
//		out_log << "row: " << b << ind;
		rows[row_index].first=b;
		rows[row_index].second=ind;
		rowlb[row_index]=-INFINITY;
		rowub[row_index]=-it->c;
		++row_index;
//		cutinfo.rowitems.push_back(solver->add_row(b, ind, -INFINITY, -it->c));
	}

	// x + w - z = ref_x constraints
	dvector b2(3, 1.); b2[2]=-1.; // (x(i), w(i), z(i))
	ivector ind2(3);
	ind2[1]=w_z_begin; // w(i)
	ind2[2]=w_z_begin+size; // z(i)
	int i0;
	for (int i=0; i<size; ++i) {
		i0=cut.indices(i);
		ind2[0]=linrelax.obj->block[cutinfo.block_nr][i0];  // x(i)
//		cutinfo.rowitems.push_back(solver->add_row(b2, ind2, cut.ref_x(i0), cut.ref_x(i0)));
		rows[row_index].first=b2;
		rows[row_index].second=ind2;
		rowlb[row_index]=cut.ref_x(i0);
		rowub[row_index]=cut.ref_x(i0);
		++row_index;
		++ind2[1]; ++ind2[2];
	}

	// w + z <= max(...) constraints
	dvector b3(2, 1.);
	ivector ind3(2);
	ind3[0]=w_z_begin; // w(i)
	ind3[1]=w_z_begin+size; // z(i)
	for (int i=0; i<size; ++i) {
//		cutinfo.rowitems.push_back(solver->add_row(b3, ind3, -INFINITY, MAX(cut.w_z_up(i),cut.w_z_up(size+i))));
		rows[row_index].first=b3;
		rows[row_index].second=ind3;
		rowlb[row_index]=-INFINITY;
		rowub[row_index]=MAX(cut.w_z_up(i),cut.w_z_up(size+i));
		++row_index;
		++ind3[0]; ++ind3[1];
	}
	
	solver->add_rows(cutinfo.rowitems, rows, rowlb, rowub);
}

void LinearRelaxSolverMIP::construct(Pointer<MinlpNode> node) {
	reset();
	linrelax.get_cuts(cutinfos, node);
	for (list<CutPool::CutInfo>::iterator it(cutinfos.begin()); it!=cutinfos.end(); it++)
		switch (it->type) {
			case CutPool::SIMPLE:
			case CutPool::LINEARIZATION: {
				const SimpleCut& cut(it->type==CutPool::SIMPLE ? it->it_simplecuts->get_cut() : it->it_linearizationcuts->get_cut());
				Pointer<UserVector<double> > b(cut.coeff->getcopy()); *b*=2;
				if (it->block_nr>=0) it->rowitems.push_back(solver->add_row(*b, linrelax.obj->block[it->block_nr], -INFINITY, -cut.constant));
				else { b->resize(solver->nr_col());
					it->rowitems.push_back(solver->add_row(*b, -INFINITY, -cut.constant));
//					out_log << -cut.constant << " >= " << *b;
				}
			} break;
			case CutPool::INTERVALGRADIENT: {
				add_intervalgradientcut(*it);
			} break;
			default: out_err << "LinearRelaxSolverMIP: CutType unknown. Aborting." << endl; exit(-1);
		}

	if (node) { // set bounds
		for (int i=0; i<node->lower.dim(); ++i)
			solver->modify_col(i, node->lower(i), node->upper(i), mip.getColType(i));
	}
}

void LinearRelaxSolverMIP::add_cut(const CutPool::CutInfo& cutinfo) {
	cutinfos.push_back(cutinfo);
	switch (cutinfo.type) {
		case CutPool::SIMPLE:
		case CutPool::LINEARIZATION: {
			const SimpleCut& cut(cutinfo.type==CutPool::SIMPLE ? cutinfo.it_simplecuts->get_cut() : cutinfo.it_linearizationcuts->get_cut());
			Pointer<UserVector<double> > b(cut.coeff->getcopy()); *b*=2;
			if (cutinfo.block_nr>=0) cutinfos.back().rowitems.push_back(solver->add_row(*b, linrelax.obj->block[cutinfo.block_nr], -INFINITY, -cut.constant));
			else { b->resize(solver->nr_col());
				cutinfos.back().rowitems.push_back(solver->add_row(*b, -INFINITY, -cut.constant));
//				out_log << -cutinfo.it_cuts->get_cut().second << " >= " << *b;
			}
		} break;
		case CutPool::INTERVALGRADIENT: {
			add_intervalgradientcut(cutinfos.back());
		} break;
		default: out_err << "LinearRelaxSolverMIP: CutType unknown. Aborting." << endl; exit(-1);
	}
}

void LinearRelaxSolverMIP::set_objective(Pointer<SepQcFunc> obj) {
	SparseVector<double> b(obj->b, obj->block); b*=2;
	b.resize(solver->nr_col());
	solver->set_obj(b, obj->c);
}

void LinearRelaxSolverMIP::update_levelcut(double newlevel) {
	if (levelcut_pos==-1) return;
	mip.setRowBounds(levelcut_pos, -INFINITY, newlevel-linrelax.obj->c);
	solver->modify_row(levelcut_pos, -INFINITY, newlevel-linrelax.obj->c);
}

int LinearRelaxSolverMIP::solve(Pointer<MinlpNode> node) {
//	solver->set_tol(tol);
	MIPSolver::SolutionStatus ret=solver->solve();
	opt_val_=solver->get_optval();
//	out_log << "MIPSolver returned " << ret << "\t optval: " << opt_val_ << "\t size: " << solver->nr_col() << " cols, " << solver->nr_row() << " rows" << endl;
	solver->get_primal(sol_point);
	solver->get_dual(duals); // get the dual values for the first rows
	duals*=-1; // to make them conform with ivo's notation
	switch (ret) {
		case MIPSolver::SOLVED:
			break;
		case MIPSolver::ITERATIONLIMITEXCEEDED:
			return 3;
		case MIPSolver::UNBOUNDED:
			return 2;
		case MIPSolver::FEASIBLE:
			return 4;
		case MIPSolver::INFEASIBLE:
			return 1;
		case MIPSolver::ABORTED:
		case MIPSolver::UNKNOWN:
		default:
			return 5;
	}

	if (linrelax.inactivetime_limit_global || linrelax.inactivetime_limit_local) {
		// update inactivity info
		for (list<CutPool::CutInfo>::iterator it(cutinfos.begin()); it!=cutinfos.end(); ++it) {
			if (it->rowitems.empty()) continue; // should never happen
			switch(it->type) {
				case CutPool::SIMPLE:
				case CutPool::LINEARIZATION: it->inactive=fabs(solver->get_dual(*it->rowitems.front()))<rtol; break;
				case CutPool::INTERVALGRADIENT: { it->inactive=true;
					list<const MIPSolver::RowItem*>::iterator itrow(it->rowitems.begin());
					for (int c=0; it->inactive && c<it->it_intgradcuts->get_cut().coninfos_size; ++c, ++itrow)
						it->inactive&=fabs(solver->get_dual(**itrow))<rtol;
				} break;
				default: out_err << "LinearRelaxSolverMIP: CutType unknown. Aborting." << endl; exit(-1);
			}
		}
	}

	return 0; // solved
}

int LinearRelaxSolverMIP::generate_cuts(list<Pointer<SimpleCut> >& cuts) {
	int nr=solver->generate_cuts(cuts);
	if (intervalgradcutvars_start==solver->nr_col()) return nr; // nothing todo
	for (list<Pointer<SimpleCut> >::iterator it(cuts.begin()); it!=cuts.end(); ++it) {
		for (int i=intervalgradcutvars_start; i<solver->nr_col(); ++i) {
			double coeff=2*(*(*it)->coeff)(i);
			if (coeff>0) (*it)->constant+=coeff*solver->get_collow(i);
			else if (coeff<0) (*it)->constant+=coeff*solver->get_colup(i);
		}
		(*it)->coeff->resize(intervalgradcutvars_start);
	}
	return nr;
}

void LinearRelaxSolverMIP::remove_cuts() {
	int removed=0;
	list<const MIPSolver::RowItem*> solver_remove; // the list for the MIPSolver of the RowItem's to remove
	list<const MIPSolver::ColItem*> solver_remove2; // the list for the MIPSolver of the ColItem's to remove
	for (list<CutPool::CutInfo>::iterator it(cutinfos.begin()); it!=cutinfos.end(); )
		if (it->removed) {
			solver_remove.insert(solver_remove.end(), it->rowitems.begin(), it->rowitems.end());
			solver_remove2.insert(solver_remove2.end(), it->colitems.begin(), it->colitems.end());
			it=cutinfos.erase(it);
			removed++;
		} else ++it;
	solver->delete_rows(solver_remove);
	solver->delete_cols(solver_remove2);
//	if (removed) out_log << "Deleted " << removed << " cuts." << endl;
}

bool LinearRelaxSolverMIP::feasible() {
	return solver->feasible()!=MIPSolver::INFEASIBLE;
}
