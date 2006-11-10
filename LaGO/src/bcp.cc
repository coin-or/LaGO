// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "bcp.h"
#include "column.h"
#include "lagheu.h"

extern "C" {
#include "ranlib.h"
}

// ----------------------------- ExtremePoint ------------------------------------

void ExtremePoint::set_rmpdata(Pointer<LinearRelax> linear_relax, int k) {
	rmpcolumn=0.; // clear
	rmpcolumn.resize(linear_relax->core_size());

	int c=0;
	for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->couple_con.begin()); it!=linear_relax->couple_con.end(); it++, c++)
		if (it->b[k]) rmpcolumn.SetElement(c, 2.* *it->b[k] * *this);
	for (int k_=0; k_<k; k_++) c+=linear_relax->block_con[k_].size();
	for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->block_con[k].begin()); it!=linear_relax->block_con[k].end(); it++, c++)
		rmpcolumn.SetElement(c, 2.* *it->b[0] * *this);

	rmpobjcoeff=linear_relax->obj->b[k] ? 2.* *linear_relax->obj->b[k] * *this : 0.;
}

// --------------------------- BCP Algorithm -------------------------------------

MinlpBCP::MinlpBCP(Pointer<MinlpProblem> orig_prob_, Pointer<MinlpProblem> split_prob_, Pointer<LinearRelax> linear_relax_,
	bool is_gams_prob, double closeval_tol_, Pointer<dvector> diam_, Pointer<Param> param_,
	Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: RelaxationSolver(orig_prob_, is_gams_prob, closeval_tol_, diam_, param_, out_solver_p_, out_solver_log_p_),
	intgrad_cutgen(split_prob_), linconcutgen(split_prob_), ExtremePoints(split_prob_->block.size()),
	bound_impr_tol(.1), bound_failed(0), bound_computed(0), lagprob_solves(0), bound_time(0.), init_RMP_time(0), max_time(INFINITY),
	find_solcand_time(0.), nr_solcand_found(0),
	subdiv_time(0.), update_ExtremePoints_count(0), nr_subdiv_contvar(0), nr_subdiv_bisect(0),
	prob_is_convex(false), timer(new Timer())
{	set_split_prob(split_prob_);
	set_linear_relax(linear_relax_);
	init();
}

MinlpBCP::~MinlpBCP() { };

void MinlpBCP::set_convex_relax(Pointer<MinlpProblem> convex_prob_, Pointer<dvector> sol_C_, bool sol_C_is_solution_) {
	RelaxationSolver::set_convex_prob(convex_prob_);
	sol_C=sol_C_;
	sol_C_is_solution=sol_C_is_solution_;
	
	if (!reform) linconcutgen.set_problem(convex_prob);
}

void MinlpBCP::set_reform(Pointer<Reformulation> reform_, Pointer<dvector> sol_C_, bool sol_C_is_solution_) {
	set_reform(reform_);
	sol_C=sol_C_;
	sol_C_is_solution=sol_C_is_solution_;
}

void MinlpBCP::set_reform(Pointer<Reformulation> reform_) {
	RelaxationSolver::set_reform(reform_);
	
	linconcutgen.set_problem(reform->ext_convex_prob);
	linconcutgen.set_reform(reform);
	intgrad_cutgen.set_problem(reform->ext_prob);
	if (lagheu) lagheu->set_reform(reform);
	if (intervalreduction) intervalreduction->set_problem(reform->ext_prob);
}

void MinlpBCP::set_linear_relax(Pointer<LinearRelax> linrelax_) {
	RelaxationSolver::set_linear_relax(linrelax_);
	if (lagheu) lagheu->set_linear_relax(linrelax_);
}

void MinlpBCP::set_MINLPData(Pointer<MINLPData> minlpdata_) {
	RelaxationSolver::set_MINLPData(minlpdata_);
	
	linconcutgen.set_MINLPData(minlpdata);
}

void MinlpBCP::set_timer(Pointer<Timer> timer_) {
	timer=timer_;
}


void MinlpBCP::init() {
	iter_max=param->get_i("MinlpBCP max iter", 10000);
	max_time=param->get_d("MinlpBCP max time", 3600); // one hour
	gap_tol=param->get_d("MinlpBCP gap tol", 0.01);
	tol=1E-4;
	lag_cuts=param->get_i("Lagrangian cuts", 1);
	intgrad_cuts=param->get_i("IntervalGradient cuts", 0);
	mip_cuts=param->get_i("MIP cuts", 1);

	mem_limit=param->get_i("Memory limit", 0);

	if (!strcmp(param->get("BCP bound type", "LP"), "LP")) maj_bound_type=LP_bound;
	else if (!strcmp(param->get("BCP bound type"), "NLP")) maj_bound_type=NLP_bound;
	else if (!strcmp(param->get("BCP bound type"), "RMP")) maj_bound_type=RMP_bound;
	else if (!strcmp(param->get("BCP bound type"), "LP-RMP")) maj_bound_type=LP_RMP_bound;
	else if (!strcmp(param->get("BCP bound type"), "stop")) maj_bound_type=stop_bound;
	else maj_bound_type=LP_bound;

	pre_bb_max_iter=param->get_i("BCP preprocess max iter", 0);
	pre_bound_type=RMP_bound;

	lagsolve_type=BranchCut;

	if (!strcmp(param->get("BCP subdiv type", "Binary"), "Binary")) subdiv_type=BinSubdiv;
	else if (!strcmp(param->get("BCP subdiv type"), "Cost")) subdiv_type=CostSubdivLag;
	else if (!strcmp(param->get("BCP subdiv type"), "Bisection")) subdiv_type=BisectSubdiv;
	else if (!strcmp(param->get("BCP subdiv type"), "Violation")) subdiv_type=ViolSubdiv;
	else subdiv_type=CostSubdivLag;

	is_maxcut=param->get_i("maxcut", 0);

	upper_bound_effort_level=param->get_i("BCP upper bound effort", is_maxcut ? 2 : 0);

//	conv_rate_cntrl=param.get_i("control convergence rate", 1);
	conv_rate_cntrl_stopping_rho=param->get_d("stopping rho", 0.1);
	conv_rate_cntrl_minor_iter=param->get_i("minor iterations", 5);
	conv_rate_cntrl_last_major_val=-INFINITY;
	conv_rate_cntrl_last_val=-INFINITY;
	conv_rate_cntrl_improve_iter=-1;

//	if (split_prob->block.size()>1) bound_print=new ofstream("bounds.log");

	if (param->get_i("BCP IntervalReduction", 1)) {
		intervalreduction=new IntervalReduction();
//		intervalreduction=new IntervalReduction(orig_prob);
//		intervalreduction->do_print=true;
	}

	opt_val_=INFINITY;
}

void MinlpBCP::init_block_problems() {
	Pointer<MinlpProblem> prob(reform ? reform->ext_prob : split_prob);
	Pointer<MinlpProblem> conv(reform ? reform->ext_convex_prob : convex_prob);
	assert(prob->block.size()==conv->block.size());
	
	block_prob.resize(prob->block.size());
	block_convex_prob.resize(prob->block.size());

	for (int k=0; k<block_prob.size(); k++) {
		block_prob[k]=new MinlpProblem(*prob, k);
		block_convex_prob[k]=new MinlpProblem(*conv, k);
		for (int c=0; c<prob->con.size(); c++) {
			if (prob->con[c]->A[k] || prob->con[c]->s[k]) {
				block_prob[k]->add_con(new SepQcFunc(prob->con[c]->A[k], prob->con[c]->b[k], prob->con[c]->s[k], prob->con[c]->c), prob->con_eq[c], prob->con_names[c]);
				if (conv->con[c]->A[k] || conv->con[c]->s[k] || conv->con[c]->b[k])
					block_convex_prob[k]->add_con(new SepQcFunc(conv->con[c]->A[k], conv->con[c]->b[k], conv->con[c]->s[k], conv->con[c]->c), conv->con_eq[c], conv->con_names[c]);
				else {
					block_convex_prob[k]->add_con(new SepQcFunc(prob->block[k].size()));
					block_convex_prob[k]->con.back()->c=conv->con[c]->c;
				}
			}
		}

		for (int c=prob->con.size(); c<conv->con.size(); c++)
			if (conv->con[c]->A[k] || conv->con[c]->s[k])
				block_convex_prob[k]->add_con(new SepQcFunc(conv->con[c]->A[k], conv->con[c]->b[k], conv->con[c]->s[k], conv->con[c]->c), conv->con_eq[c], conv->con_names[c]);

		for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->block_con[k].begin()); it!=linear_relax->block_con[k].end(); it++) {
			Pointer<SepQcFunc> con(new SepQcFunc(NULL, it->b[0], NULL, it->c));
			block_prob[k]->add_con(con, it->eq, it->name);
			block_convex_prob[k]->add_con(con, it->eq, it->name);
		}
	}
}

void MinlpBCP::clean_sub_problems() {
	RelaxationSolver::clean_sub_problems();
	block_sub_convex_prob.clear();
	lag_problem.clear();
	colgen=NULL;
}

// -------------------------------------- RMP points -----------------------------

int MinlpBCP::init_ExtremePoints(Pointer<MinlpNode> node) {
	ExtremePoints.resize(split_prob->block.size());
	node->i_ExtremePoints.resize(split_prob->block.size());
	node->i_ExtremePoints_limit.resize(split_prob->block.size());
	for (int k=0; k<node->i_ExtremePoints_limit.size(); ++k)
		node->i_ExtremePoints_limit[k]=ExtremePoints[k].end();

	int ret;
	if (sol_C) {
		ret=primal_init_ExtremePoints(*sol_C, node);
		out_solver_log << "primal_init_ExtremePoints using feasible point from (C) returned: " << ret << endl;
		if (!ret) {
			project_ExtremePoints(*sol_C, node);
			return 0;
/*			if (!colgen) colgen=new ColumnGenerator(node, *this);
			if (colgen->RMP_feasible()) return 0;
			out_solver_log << "RMP not feasible, trying to generate more points." << endl;
*/		}
	}

	ret=set_LP_bound(node);
	if (ret) {
		out_solver << "Cannot initialize Extreme points" << endl;
		return ret;
	}

	ret=dual_init_ExtremePoints(node);
	out_solver_log << "dual_init_ExtremePoints using dual of (R) returned: " << ret << endl;
/*	if (!ret) {
		if (!colgen) colgen=new ColumnGenerator(node, *this);
		if (colgen->RMP_feasible()) return 0;
		out_solver_log << "RMP not feasible, trying to generate more points." << endl;
	}
*/
	ret=primal_init_ExtremePoints(node->ref_point, node);
	out_solver_log << "primal_init_ExtremePoints using feasible point from (R) returned: " << ret << endl;
	if (ret) {
		out_solver << "Cannot initialize Extreme points." << endl;
		return ret;
	}
	project_ExtremePoints(node->ref_point, node);
	return 0;
/*
	if (!colgen) colgen=new ColumnGenerator(node, *this);
	if (colgen->RMP_feasible()) return 0;
	out_solver << "RMP still not feasible." << endl;

	return 1;
*/}
//-----------------------------------------------------------------

int MinlpBCP::dual_init_ExtremePoints(Pointer<MinlpNode> node) {
	if (!lag_problem.size()) init_lag_problems(node);
	update_lag_problems(node->dual_point);
	for (int k=0; k<split_prob->block.size(); k++) {
		LagSolveStatus status;
		solve_lag_problem(status, k, node);
		out_solver_log << "dual_init_ExtremePoints block " << k << "\t lagsolve return: " << status.ret << "\t points: " << status.solset.size() << "\t blocksize: " << split_prob->block[k].size() << endl;
		if (!status.solset.size()) return (status.ret ? status.ret : 1);

		for (set<SolCandidate>::iterator it(status.solset.begin()); it!=status.solset.end(); it++)
			add_ExtremePoint(it->second, k, node);
	}

	return 0;
}
//-----------------------------------------------------------------

int MinlpBCP::primal_init_ExtremePoints(const dvector& x, Pointer<MinlpNode> node) {
	Pointer<MinlpProblem> prob=reform ? reform->ext_prob : split_prob; assert(prob);
	int count_allblocks=0;
	for (int k=0; k<prob->block.size(); k++) {
		set<int> K;
		dvector x_k(x, prob->block[k]);
		dvector low(prob->lower, prob->block[k]);
		dvector up(prob->upper, prob->block[k]);
		dvector rounded(x_k);
		dvector w(x_k.dim());
		for (int i=0; i<prob->block[k].size(); i++) { // move x, if necessary, a bit inside the box
			if (up(i)-low(i)<rtol) continue;
			if (fabs(x_k(i)-low(i))<rtol) x_k[i]+=random(.1, .2)*(up(i)-low(i)); // on lower bound
			else if (fabs(up(i)-x_k(i))<rtol) x_k[i]-=random(.1, .2)*(up(i)-low(i)); // on upper bound
			else if (fabs(x_k(i)-.5*(low(i)+up(i)))<rtol) x_k[i]+=random(-.1, .1)*(up(i)-low(i)); // on center
   		K.insert(i);
		}

		int count=0;
		while (K.size()) {
			double t=INFINITY;
			double newval;
			set<int>::iterator j(K.end());
			for (set<int>::iterator it(K.begin()); it!=K.end(); it++) {
				rounded[*it] = x_k(*it)<=0.5*(low(*it)+up(*it)) ? low(*it) : up(*it);
				newval=(1-t)*rounded(*it)+t*x_k(*it);
				if (newval<low(*it)) {
					t=(low(*it)-rounded(*it))/(x_k(*it)-rounded(*it));
					j=it;
				}	else if (newval>up(*it)) {
					t=(up(*it)-rounded(*it))/(x_k(*it)-rounded(*it));
					j=it;
				}
			}
			assert(j!=K.end());
			assert(t>-rtol);
			w=rounded;
			if ((!reform) || reform->set_t_block(w, k)) // try to make feasible by adjusting t's
				if (add_ExtremePoint(w, k, node).second) count++;
//				out_log << k << ": added:\t" << w;

			x_k*=t; x_k.AddMult(1-t, rounded); // x = w + t (x-w)
			rounded=x_k;
			K.erase(j);
			for (set<int>::iterator it(K.begin()); it!=K.end(); )
				if (fabs(x_k(*it)-low(*it))<rtol || fabs(x_k(*it)-up(*it))<rtol) {
					set<int>::iterator next(it); next++;
					K.erase(it);
					it=next;
				} else it++;
		}
		if (reform) reform->set_t_block(x_k, k); // try to make feasible by adjusting t's
		if (add_ExtremePoint(x_k, k, node).second) count++;
//			out_log << k << ": added:\t" << x_k;

		out_log << "primal_init_ExtremePoints block " << k << "\t points: " << count << "\t blocksize: " << split_prob->block[k].size() << endl;
		count_allblocks+=count;
	}

	return (!count_allblocks);
}

//-----------------------------------------------------------------

pair<list<ExtremePoint>::iterator, bool> MinlpBCP::add_ExtremePoint(const dvector &w, int k, Pointer<MinlpNode> node) {
	bool this_is_new=true;
	list<ExtremePoint>::iterator it_old(ExtremePoints[k].begin());
	for (; it_old!=ExtremePoints[k].end(); ++it_old) {
		this_is_new=false;
		for (int i=0; (!this_is_new) && i<w.dim(); i++)
			this_is_new=fabs((*it_old)(i)-w(i))/(fabs(w(i))+1)>0.01;
		if (!this_is_new) break;
	}
//out_log << "Point is new: " << this_is_new << endl;
	if (!this_is_new) {
		if (!node) return pair<list<ExtremePoint>::iterator, bool>(it_old, false);
		this_is_new=true;
		for (list<list<ExtremePoint>::iterator>::iterator it_node(node->i_ExtremePoints[k].begin()); this_is_new && it_node!=node->i_ExtremePoints[k].end(); ++it_node)
			this_is_new=(*it_node!=it_old);
		if (!this_is_new) // point already known for node
			return pair<list<ExtremePoint>::iterator, bool>(it_old, false);

		node->i_ExtremePoints[k].push_back(it_old); //??? conflict with pruning method
		if (node->i_ExtremePoints_limit[k]==ExtremePoints[k].end()) // first point
			node->i_ExtremePoints_limit[k]=it_old;
		else { // if point was the next one to check
			list<ExtremePoint>::iterator next(node->i_ExtremePoints_limit[k]); next++;
			if (it_old==next) ++node->i_ExtremePoints_limit[k];
			else node->i_ExtremePoints_limit[k]=--ExtremePoints[k].end();
		}
		return pair<list<ExtremePoint>::iterator, bool>(it_old, true);
	}

	ExtremePoints[k].push_back(ExtremePoint(w, linear_relax, k));
	if (node) {
		node->i_ExtremePoints_limit[k]=--ExtremePoints[k].end();
		node->i_ExtremePoints[k].push_back(node->i_ExtremePoints_limit[k]);
	}

	return pair<list<ExtremePoint>::iterator, bool>(--ExtremePoints[k].end(), true);
}

//-----------------------------------------------------------------

void MinlpBCP::prune_ExtremePoints(multimap<double, Pointer<MinlpNode> >& bb_tree2) {
	if (!bb_tree.size() && !bb_tree2.size()) return; // or prune all points?
	int count=0;
	int points=0;
	for (int k=0; k<ExtremePoints.size(); k++) {
		if (!ExtremePoints[k].size()) continue;
		points+=ExtremePoints[k].size();

		map<ExtremePoint*, list<ExtremePoint>::iterator> prunepool;
		for (list<ExtremePoint>::iterator it(ExtremePoints[k].begin()); it!=ExtremePoints[k].end(); ++it)
			prunepool.insert(pair<ExtremePoint*, list<ExtremePoint>::iterator>(&*it, it));

		multimap<double, Pointer<MinlpNode> >::iterator it_tree(bb_tree.begin());
		while (it_tree!=bb_tree2.end()) {
			if (it_tree==bb_tree.end()) {
				if (!bb_tree2.size()) break;
				it_tree=bb_tree2.begin();
			}
			for (list<list<ExtremePoint>::iterator>::iterator it_pt(it_tree->second->i_ExtremePoints[k].begin()); it_pt!=it_tree->second->i_ExtremePoints[k].end(); it_pt++)
				prunepool.erase(&**it_pt);
			it_tree++;
		}

		// check "last checked"-elements
		it_tree=bb_tree.begin();
		while (it_tree!=bb_tree2.end()) {
			if (it_tree==bb_tree.end()) {
				if (!bb_tree2.size()) break;
				it_tree=bb_tree2.begin();
			}
			while (it_tree->second->i_ExtremePoints_limit[k]!=ExtremePoints[k].end() && prunepool.count(&*it_tree->second->i_ExtremePoints_limit[k])) {
				if (it_tree->second->i_ExtremePoints_limit[k]==ExtremePoints[k].begin()) { // all RMP points removed
					it_tree->second->i_ExtremePoints_limit[k]=ExtremePoints[k].end();
					break;
				}
				--it_tree->second->i_ExtremePoints_limit[k];
			}
			it_tree++;
		}
		for (map<ExtremePoint*, list<ExtremePoint>::iterator>::iterator it(prunepool.begin()); it!=prunepool.end(); it++) {
			ExtremePoints[k].erase(it->second);
			count++;
		}
	}
	if (count) out_solver_log << "Pruned " << count << " unused Extreme points out of " << points << endl;
}

//-----------------------------------------------------------------

void MinlpBCP::project_ExtremePoints(dvector& x, Pointer<MinlpNode> node) {
	if (!lag_problem.size()) init_lag_problems(node);
	int count=0;
	int ret;
	for (int k=0; k<ExtremePoints.size(); ++k) {
		dvector mid(lag_problem[k]->dim());
		int size=0;
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it, ++size)
			mid+=**it;
		mid/=size;
		out_log << "project " << size << " extreme points block " << k << ": ";
		list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin());
		for (int i=0; i<size; ++i) {
			if (!lag_problem[k]->feasible(**it, tol, NULL)) { ++it; out_log << 'f'; continue; } // feasible point
			dvector a(mid);
			a-=**it;
			update_lag_problem(k, a);
			LagSolveStatus status;
			solve_lag_problem(status, k, node);
			out_log << status.ret;
			if (status.ret) {
				it++;
				continue;
			}
			pair<list<ExtremePoint>::iterator, bool> ret(add_ExtremePoint(status.solset.begin()->second, k, node));
			if (!ret.second) { // point not new
				if (*it==ret.first) { // because it didn't change
					out_log << 'n';
					++it;
				} else {
					out_log << 'e';
					it=node->i_ExtremePoints[k].erase(it);
					--size;
				}
				count++;
				continue;
			}
			it=node->i_ExtremePoints[k].erase(it); // new point was added to node
//			*it=ret.first; // replace old by new point
			node->i_ExtremePoints_limit[k]=ret.first;
			count++;
//			it++;
		}
		out_log << "\t -> " << size << endl;
	}
	out_log << "Projected " << count << " extreme points." << endl;
}

//-----------------------------------------------------------------------
// Bounding
//-----------------------------------------------------------------------

int MinlpBCP::set_low_bound(Pointer<MinlpNode> node) {
	int ret=0;
	Timer t;
	switch(bound_type) {
		case NLP_bound : ret=set_NLP_bound(node); break;
  	case RMP_bound : ret=set_RMP_bound(node); break;
		case LP_bound  : ret=set_LP_bound(node); break;
		case LP_RMP_bound : ret=set_LP_RMP_bound(node); break;
	}
	bound_time+=t.stop();

	if (!ret) bound_computed++;
	return ret==1 ? 1 : 0;
}

//-----------------------------------------------------------------

pair<bool, double> MinlpBCP::improve_bound(Pointer<MinlpNode> node) {
	double oldbound=node->low_bound;
	pair<bool, double> ret(true, 0.);
	int subret;

	if (node->update_subdiv_bound_called) return ret;

	if (intgrad_cuts) {
		MinlpProblem& prob(reform ? *reform->ext_prob : *split_prob);
		int added=0;
		for (int k=0; k<prob.block.size(); ++k) {
			Pointer<IntervalGradientCut> intgradcut(intgrad_cutgen.get_cuts(node->ref_point(prob.block[k]), k, node->lower(prob.block[k]), node->upper(prob.block[k])));
			if (intgradcut) {
				added+=intgradcut->coninfos_size;
				linear_relax->add_cut(intgradcut, k, node);
			}
		}
		if (added) out_solver_log << "Added " << added << " IntervalGradientCuts." << endl;
	} //else return ret;
//	if (!mip_cuts) return ret;


	Timer t;
	switch(bound_type) {
	  case NLP_bound : ret.first=!set_NLP_bound(node, true); break;
	  case RMP_bound : ret.first=!set_RMP_bound(node); break;
	  case LP_RMP_bound : ret.first=!set_LP_RMP_bound(node); break;
	  case LP_bound  : ret.first=!improve_LP_bound(node); break;
		case stop_bound : break;
		default: out_err << "Bound type not known. Aborting" << endl; exit(-1); // panic
	}

	bound_time+=t.stop();

	if (ret.first) { // bound computed
		bound_computed++;
		ret.second=(node->low_bound-oldbound)/(1+fabs(oldbound));
	}

	return ret;
}

//-----------------------------------------------------------------

int MinlpBCP::set_NLP_bound(Pointer<MinlpNode> node, bool improve) {
	bool loccuts=false;
//	for (int k=0; (!loccuts) && k<node->local_cuts.size(); k++) loccuts=!node->local_cuts[k].empty();
	for (int k=0; (!loccuts) && k<node->part_con.size(); k++) loccuts=!node->part_con[k].empty();

	if ((!sub_convex_prob) || loccuts) {
		sub_convex_prob=NULL;
		make_sub_prob(sub_convex_prob, reform ? reform->ext_convex_prob : convex_prob, node);
// reform sub_convex_prob by moving nonlinear constraints as summand to objective
//out_log << *sub_convex_prob;
if (is_maxcut) {
		sub_convex_prob->obj=new SepQcFunc(*sub_convex_prob->obj);
		int c=0;
		for (int k=0; k<sub_convex_prob->block.size(); k++)
			if (!sub_convex_prob->obj->b[k]) sub_convex_prob->obj->b[k]=new SparseVector<double>(sub_convex_prob->block[k].size());
			else sub_convex_prob->obj->b[k]=new SparseVector<double>(*sub_convex_prob->obj->b[k]);
		while (c<sub_convex_prob->con.size()) {
			bool nonlin=false;
			for (int k=0; (!nonlin) && k<sub_convex_prob->block.size(); k++)
				nonlin=sub_convex_prob->con[c]->A[k] || sub_convex_prob->con[c]->s[k];
			if (!nonlin) { c++; continue; }
			for (int k=0; k<sub_convex_prob->block.size(); k++) {
				if (sub_convex_prob->con[c]->A[k])
					sub_convex_prob->obj->A[k]=new SumMatrix(Pointer<const UserMatrix>(sub_convex_prob->con[c]->A[k]), Pointer<const UserMatrix>(sub_convex_prob->obj->A[k]));
				if (sub_convex_prob->con[c]->s[k])
					sub_convex_prob->obj->s[k]=new SumFunc(sub_convex_prob->con[c]->s[k], sub_convex_prob->obj->s[k]);
				if (sub_convex_prob->con[c]->b[k])
					*sub_convex_prob->obj->b[k]+=*sub_convex_prob->con[c]->b[k];
			}
			sub_convex_prob->obj->c+=sub_convex_prob->con[c]->c;
			sub_convex_prob->del_con(c);
		}
}
	} else { // just update convex relaxation
		sub_convex_prob->lower=node->lower;
		sub_convex_prob->upper=node->upper;
		for (int i=0; i<node->lower.dim(); ++i)
			if (node->lower(i)>sub_convex_prob->primal_point(i)) sub_convex_prob->primal_point[i]=node->lower(i);
			else if (node->upper(i)<sub_convex_prob->primal_point(i)) sub_convex_prob->primal_point[i]=node->upper(i);
	}

	Pointer<LocOpt> locopt=LocOpt::get_solver(sub_convex_prob, param, "ConvexSolve", NULL, NULL);
//	locopt->iter_max=20*(sub_convex_prob->con.size()+sub_convex_prob->dim());
	int ret=locopt->solve(node->ref_point);
	if (ret==1) {
		out_solver_log << "Solving (Cext[U]) return: 1" << endl;
		return 1;
	}

	out_solver_log << "Solving (Cext[U]) return: " << ret << "\t value: " << locopt->opt_val() << endl;
	if (!ret) {
		node->low_bound=MAX(node->low_bound, locopt->opt_val()-1E-2);
		node->ref_point=locopt->sol_point;
//		node->dual_point=locopt->get_lag_multipliers(); // not sure, if correct :( causes a bug!
/*set_LP_bound(node);*/		return 0;
	}
	if (ret==1) return 1;
/*if (ret==3) {
	out_log << *sub_convex_prob;
	exit(0);
}
*//*	if (!sub_convex_prob->feasible(locopt->sol_point, tol, NULL)) {
		out_solver_log << "yes" << endl;
		node->ref_point=locopt->sol_point;
		node->dual_point=locopt->get_lag_multipliers(); // not sure, if correct :(
		return 0;
	}

	out_solver_log << "no" << endl;
*/	bound_failed++;

	if (improve) return improve_LP_bound(node);
	return set_LP_bound(node);
}

int MinlpBCP::set_RMP_bound(Pointer<MinlpNode> node) {
	if (!colgen) {
		colgen=new ColumnGenerator(node, *this);
		init_lag_problems(node);
	} else
		colgen->clear_RMP();

	return colgen->generate_RMP();
}

int MinlpBCP::set_LP_RMP_bound(Pointer<MinlpNode> node) {
	int ret=set_LP_bound(node);
	if (node->low_bound>opt_val()) return 0;
	if (ret) return ret;

	if (!colgen) {
		colgen=new ColumnGenerator(node, *this);
		init_lag_problems(node);
	}

	ret=colgen->solve_RMP();
	if (ret) out_solver_log << "RMP solve return: " << ret << endl;
	if (ret==1) return 1;

	return 0;
}

int MinlpBCP::set_LP_bound(Pointer<MinlpNode> node) {
	double lower;
	int ret=linear_relax->solve(node->ref_point, lower, node, &node->dual_point);
	out_solver_log << "Solving (R[U]) return: " << ret << "\t value: " << lower << endl;
	if (ret) return ret;

	node->low_bound=MAX(node->low_bound, lower);
	return 0;
}

//-----------------------------------------------------------------

int MinlpBCP::improve_LP_bound(Pointer<MinlpNode> node) {
	// generate new linearization cuts, ref_point was updated at end of the subdivision which created this node
	if (linear_relax->cutlimit_reached()) return 0;

	Pointer<MinlpProblem> prob;
	Pointer<MinlpProblem> conv;
	if (reform) {
		prob=reform->ext_prob;
		conv=reform->ext_convex_prob;
	} else {
		prob=split_prob ? split_prob : orig_prob;
		conv=convex_prob;
	}

	if (conv) { // adding linearization cuts
		list<pair<LinearizationCut, pair<int, bool> > > cuts;
		
		Project::project(node->ref_point, node->ref_point, node->lower, node->upper);
		linconcutgen.get_cuts(cuts, node->ref_point, node->lower, node->upper, true);
		
		int local_cuts_nr=0, global_cuts_nr=0;
		for (list<pair<LinearizationCut, pair<int, bool> > >::iterator cutit(cuts.begin()); cutit!=cuts.end() && !linear_relax->cutlimit_reached(); ++cutit) {
			linear_relax->add_cut(Pointer<LinearizationCut>(new LinearizationCut(cutit->first)), cutit->second.first, cutit->second.second ? NULL : node);
			if (cutit->second.second) ++global_cuts_nr; else ++local_cuts_nr;
		}
			
		out_log << "Added " << global_cuts_nr << " global and " << local_cuts_nr << " local LinearizationCuts." << endl;
	}
	
	if (mip_cuts) {
		int ret=set_LP_bound(node);
		if (ret) return ret;
		linear_relax->generate_cuts(node);
	}
	return set_LP_bound(node);
}
//----------------------------------------------------------------------

bool MinlpBCP::boxreduce(Pointer<MinlpNode> node, int index, IntervalReduction::which_bound_type which_bound) {
	if (!intervalreduction) return true;

	set<pair<int, IntervalReduction::which_bound_type> > startset;
	startset.insert(pair<int, IntervalReduction::which_bound_type>(index, which_bound));

	set<int> changed_blocks; // the blocks we changed, i.e. we need to update the ExtremePoints there

	do {
		// intervalreduction
		intervalreduction->compute(node->lower, node->upper, node->lower, node->upper, startset);
		if (intervalreduction->empty_boxes) return false;
		for (set<pair<int,int> >::iterator it(intervalreduction->fixed_binaries.begin()); it!=intervalreduction->fixed_binaries.end(); ++it)
			node->bcp_fixed_var[it->first].insert(it->second);

		// updating cuts
		for (int k=0; k<intervalreduction->reduction_by_block.size(); ++k)
			if (intervalreduction->reduction_by_block[k]<1.) {
				changed_blocks.insert(k);
				if (intervalreduction->reduction_by_block[k]<.8) linear_relax->update_cuts(node, k, intgrad_cutgen, linconcutgen);
			}
					
		// reduction using linear relaxation
		linear_relax->clear_solver();
		linear_relax->boxreduce_fixed_binaries.clear();
		startset.clear();
		out_log << "LinearRelax boxreduction: ";
		for (int k=0; k<intervalreduction->reduction_by_block.size(); ++k) {
			if (intervalreduction->reduction_by_block[k]>0.8) continue;
			double red=linear_relax->box_reduce(node, k, (reform ? reform->ext_prob : split_prob)->discr, false, &startset);
			out_log << 'b' << k;
			if (red==-INFINITY) {
				out_log << ": Relaxation infeasible." << endl;
				return false;
			}
			if (red!=1.) out_log << ':' << ' ' << red;
			out_log << ' ';
		}
		if (!linear_relax->boxreduce_fixed_binaries.empty()) {
			out_log << "fixed binaries: ";// << linear_relax->boxreduce_fixed_binaries.size() << '\t';
			for (set<pair<int,int> >::iterator it(linear_relax->boxreduce_fixed_binaries.begin()); it!=linear_relax->boxreduce_fixed_binaries.end(); ++it) {
				out_log << split_prob->var_names[split_prob->block[it->first](it->second)] << '(' << node->lower(split_prob->block[it->first](it->second)) << ')' << ' ';
				node->bcp_fixed_var[it->first].insert(it->second);
			}
		}
		out_log << endl;
	} while (!startset.empty());

	if (node->i_ExtremePoints.size()) {
		for (set<int>::iterator it_b(changed_blocks.begin()); it_b!=changed_blocks.end(); ++it_b) {
			list<dvector> newpoints;
			for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[*it_b].begin()); it!=node->i_ExtremePoints[*it_b].end();)
				if (!node->inside_part_set(**it, *it_b, linear_relax->obj->block)) {
					newpoints.push_back(**it); // project onto node
					int i0;
					for (int i=0; i<(**it).dim(); ++i) {
						i0=linear_relax->obj->block[*it_b][i];
						if ((**it)(i)<node->lower(i0)) newpoints.back()[i]=node->lower(i0);
						else if ((**it)(i)>node->upper(i0)) newpoints.back()[i]=node->upper(i0);
					}
					if (reform) reform->set_t_block(newpoints.back(), *it_b); // make feasible
					it=node->i_ExtremePoints[*it_b].erase(it);
				} else ++it;
			for (list<dvector>::iterator it(newpoints.begin()); it!=newpoints.end(); ++it) {
				add_ExtremePoint(*it, *it_b, node);
			}
//			out_log << "changed " << changed_blocks.size() << " blocks \t" << newpoints.size() << " points out of node" << endl;
		}
	}

	return true;
}

int MinlpBCP::update_subdiv_bound(int k, int i, Pointer<MinlpNode> node) {
	out_solver_log << "Updating bound after subdivision." << endl;
	Timer t;
	int ret=0;

	switch(bound_type) {
		case LP_bound: {
			node->update_subdiv_bound_called=true;
			if (!linear_relax->cutlimit_reached()) {
				Pointer<MinlpProblem> prob;
				Pointer<MinlpProblem> conv;
				if (reform) {
					prob=reform->ext_prob;
					conv=reform->ext_convex_prob;
				} else {
					prob=split_prob ? split_prob : orig_prob;
					conv=convex_prob;
				}
			
				if (conv) { // adding linearization cuts
					list<pair<LinearizationCut, pair<int, bool> > > cuts;
					
					Project::project(node->ref_point, node->ref_point, node->lower, node->upper);
					linconcutgen.get_cuts(cuts, node->ref_point, node->lower, node->upper, true);
					
					int local_cuts_nr=0, global_cuts_nr=0;
					for (list<pair<LinearizationCut, pair<int, bool> > >::iterator cutit(cuts.begin()); cutit!=cuts.end() && !linear_relax->cutlimit_reached(); ++cutit) {
						linear_relax->add_cut(Pointer<LinearizationCut>(new LinearizationCut(cutit->first)), cutit->second.first, cutit->second.second ? NULL : node);
						if (cutit->second.second) ++global_cuts_nr; else ++local_cuts_nr;
					}
						
					out_log << "Added " << global_cuts_nr << " global and " << local_cuts_nr << " local LinearizationCuts." << endl;
				}

				if (mip_cuts) {
					ret=set_LP_bound(node);
					if (ret) break;
					linear_relax->generate_cuts(node);
				}
			}
			ret=set_LP_bound(node);
		} break;
	}

	bound_time+=t.stop();

	return ret;
}

//----------------------------------------------------------------
// Branching
//----------------------------------------------------------------

int MinlpBCP::subdivide(list<Pointer<MinlpNode> >&nodes, Pointer<MinlpNode> node) {
	Timer t;
	out_solver_log << "Subdividing..." << endl;
	int subdiv_var;

//  if not all binaries fixed: call bin_subdiv
	int ret=bin_subdiv(nodes, subdiv_var, node);
	if (ret==1) return subdiv_var; // linear relaxation infeasible

	if (nodes.empty() && (!prob_is_convex)) { // further subdivision only needed if no binary subdiv and problem is not convex
		switch (subdiv_type) {
			case CostSubdivLag:
			case CostSubdivNewton: cost_subdiv(nodes, subdiv_var, node); break;
			case BisectSubdiv: bisect_subdiv(nodes, subdiv_var, node); break;
			case ViolSubdiv: viol_subdiv(nodes, subdiv_var, node); break;
			case BinSubdiv: if (ret==0 && lower_bound==-INFINITY) lower_bound=node->low_bound; break; // fix this lower bound
		}
	}

	subdiv_time+=t.stop();
	return subdiv_var;
}

//----------------------------------------------------------------------

void MinlpBCP::rect_subdiv(list<Pointer<MinlpNode> >& nodes, Pointer<MinlpNode> node, int k_star, int i_star, double cut) {
	int i0=split_prob->block[k_star][i_star];
	out_solver_log << "Subdivide at variable ";
	if (split_prob->var_names[i0]) { out_solver_log << split_prob->var_names[i0]; }
	else out_solver_log << i0;
	out_solver_log << "\t binary: " << (split_prob->discr[i0] ? "yes" : "no") << "\t cut at: " << cut << endl;

	Pointer<MinlpNode> left=new MinlpNode(*node);
	Pointer<MinlpNode> right=new MinlpNode(*node);
	if (split_prob->discr[i0]) {
		if (!node->bcp_fixed_var.count(k_star)) {
			set<int> fix; fix.insert(i_star);
			left->bcp_fixed_var.insert(pair<int, set<int> >(k_star, fix));
			right->bcp_fixed_var.insert(pair<int, set<int> >(k_star, fix));
		} else {
			left->bcp_fixed_var.find(k_star)->second.insert(i_star);
			right->bcp_fixed_var.find(k_star)->second.insert(i_star);
		}
		left->ref_point[i0]=left->upper[i0]=left->lower[i0];
		right->ref_point[i0]=right->lower[i0]=right->upper[i0];
	} else {
		nr_subdiv_contvar++;
		left->upper[i0]=cut;
		right->lower[i0]=cut;
	}
//init_lag_problems(node);
	if (node->i_ExtremePoints.size()) {
		left->i_ExtremePoints.resize(node->i_ExtremePoints.size());
		right->i_ExtremePoints.resize(node->i_ExtremePoints.size());
		left->i_ExtremePoints_limit=node->i_ExtremePoints_limit;
		right->i_ExtremePoints_limit=node->i_ExtremePoints_limit;
		for (int k=0; k<node->i_ExtremePoints.size(); ++k) {
			if (k==k_star) continue;
			left->i_ExtremePoints[k]=node->i_ExtremePoints[k];
			right->i_ExtremePoints[k]=node->i_ExtremePoints[k];
		}
		list<dvector> addleft, addright;
		if (split_prob->discr[i0]) {
			for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k_star].begin()); it!=node->i_ExtremePoints[k_star].end(); ++it) {
				if ((**it)(i_star)<split_prob->lower(i0)+rtol) { // point belongs to left node
					left->i_ExtremePoints[k_star].push_back(*it);
					addright.push_back(**it);
					addright.back()[i_star]=split_prob->upper(i0); // project onto right node
					if (reform) reform->set_t_block(addright.back(), k_star); // make feasible
				} else {
					right->i_ExtremePoints[k_star].push_back(*it);
					addleft.push_back(**it);
					addleft.back()[i_star]=split_prob->lower(i0); // project onto right node
					if (reform) reform->set_t_block(addleft.back(), k_star); // make feasible
				}
			}
		} else {
			for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k_star].begin()); it!=node->i_ExtremePoints[k_star].end(); ++it) {
				if ((**it)(i_star)<cut-rtol) { // point belongs to left node, and not right one
					left->i_ExtremePoints[k_star].push_back(*it);
					addright.push_back(**it); addright.back()[i_star]=cut;
					if (reform) reform->set_t_block(addright.back(), k_star); // make feasible
				} else if ((**it)(i_star)>cut+rtol) { // point belongs to right node, and not to left one
					right->i_ExtremePoints[k_star].push_back(*it);
					addleft.push_back(**it); addleft.back()[i_star]=cut;
					if (reform) reform->set_t_block(addleft.back(), k_star); // make feasible
				} else { // point feasible for both nodes
					left->i_ExtremePoints[k_star].push_back(*it);
					right->i_ExtremePoints[k_star].push_back(*it);
				}
			}
		}
		for (list<dvector>::iterator it(addleft.begin()); it!=addleft.end(); ++it) add_ExtremePoint(*it, k_star, left);
		for (list<dvector>::iterator it(addright.begin()); it!=addright.end(); ++it) add_ExtremePoint(*it, k_star, right);
		left->i_ExtremePoints_limit[k_star]=right->i_ExtremePoints_limit[k_star]=--ExtremePoints[k_star].end();

/*
			for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k_star].begin()); it!=node->i_ExtremePoints[k_star].end(); ++it) {
				dvector proj(**it); // projected RMP point
				if (proj(i_star)<split_prob->lower(i0)+rtol) { // point belongs to left node
					left->i_ExtremePoints[k_star].push_back(*it);

					proj[i_star]=split_prob->upper(i0); // project onto right node
					if (reform) reform->set_t_block(proj, k_star); // make feasible
					pair<list<ExtremePoint>::iterator, bool> ret(add_ExtremePoint(proj, k_star)); // add to ExtremePoints pool
					if (ret.second) { // if new at pool, add to right node
						right->i_ExtremePoints[k_star].push_back(ret.first);
						right->i_ExtremePoints_limit[k_star]=ret.first;
					}
				} else {
					right->i_ExtremePoints[k_star].push_back(*it);

					proj[i_star]=split_prob->lower(i0); // project onto left node
					if (reform) reform->set_t_block(proj, k_star); // make feasible
					pair<list<ExtremePoint>::iterator, bool> ret(add_ExtremePoint(proj, k_star)); // add to ExtremePoints pool
					if (ret.second) {  // if new at pool, add to right node
						left->i_ExtremePoints[k_star].push_back(ret.first);
						left->i_ExtremePoints_limit[k_star]=ret.first;
					}
				}
			}
		} else {
			for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k_star].begin()); it!=node->i_ExtremePoints[k_star].end(); ++it) {
				dvector proj(**it); proj[i_star]=cut; // project onto cut
				if (reform) reform->set_t_block(proj, k_star); // make feasible
				if ((**it)(i_star)<cut-rtol) { // point belongs to left node, and not right one
					left->i_ExtremePoints[k_star].push_back(*it);
					left->i_ExtremePoints_limit[k_star]=*it;
					pair<list<ExtremePoint>::iterator, bool> ret(add_ExtremePoint(proj, k_star)); // add to ExtremePoints pool
					if (ret.second) {
						right->i_ExtremePoints[k_star].push_back(ret.first); // if new at pool, add to right node
						right->i_ExtremePoints_limit[k_star]=ret.first;
					}
				} else if ((**it)(i_star)>cut+rtol) { // point belongs to right node, and not to left one
					right->i_ExtremePoints[k_star].push_back(*it);
					pair<list<ExtremePoint>::iterator, bool> ret(add_ExtremePoint(proj, k_star)); // add to ExtremePoints pool
					if (ret.second) {  // if new at pool, add to left node
						left->i_ExtremePoints[k_star].push_back(ret.first);
						left->i_ExtremePoints_limit[k_star]=ret.first;
					}
				} else { // point feasible for both nodes
					left->i_ExtremePoints[k_star].push_back(*it);
					right->i_ExtremePoints[k_star].push_back(*it);
				}
			}
		}
*/	}
	clean_sub_problems();

	MinlpProblem& prob(reform ? *reform->ext_prob : *split_prob);

	out_solver_log << "Updating cuts and boxreduction for left node: ";
	if (split_prob->var_names[i0]) { out_solver_log << split_prob->var_names[i0]; }	else out_solver_log << i0;
	out_solver_log << " = [" << left->lower(i0) << ',' << left->upper(i0) << ']' << endl;
	linear_relax->duplicate_nodeinfo(node, left);
	linear_relax->update_cuts(left, k_star, intgrad_cutgen, linconcutgen); // update existing intervalgradient cuts
	if (left->ref_point(i0)<=cut) { // refpoint in left node -> generate new intervalgradient cuts
		if (intgrad_cuts) {
			Pointer<IntervalGradientCut> intgradcut(intgrad_cutgen.get_cuts(left->ref_point(prob.block[k_star]), k_star, left->lower(prob.block[k_star]), left->upper(prob.block[k_star])));
			linear_relax->add_cut(intgradcut, k_star, left);
			if (intgradcut && intgradcut->coninfos_size) out_solver_log << "Added " << intgradcut->coninfos_size << " IntervalGradientCuts." << endl;
		}
		if (update_subdiv_bound(k_star, i_star, left)) {
			linear_relax->remove_node(left);
			left=NULL;
		}
	} else left->ref_point[i0]=cut;
	if (left && !boxreduce(left, i0, IntervalReduction::UPPER)) {
		linear_relax->remove_node(left);
		left=NULL;
	}
	if (left) nodes.push_back(left);

	out_solver_log << "Updating cuts and boxreduction for right node: ";
	if (split_prob->var_names[i0]) { out_solver_log << split_prob->var_names[i0]; }	else out_solver_log << i0;
	out_solver_log << " = [" << right->lower(i0) << ',' << right->upper(i0) << ']' << endl;
	linear_relax->duplicate_nodeinfo(node, right);
	linear_relax->update_cuts(right, k_star, intgrad_cutgen, linconcutgen);
	if (right->ref_point(i0)>=cut) { // refpoint in right node
		if (intgrad_cuts) {
			Pointer<IntervalGradientCut> intgradcut(intgrad_cutgen.get_cuts(right->ref_point(prob.block[k_star]), k_star, right->lower(prob.block[k_star]), right->upper(prob.block[k_star])));
			linear_relax->add_cut(intgradcut, k_star, right);
			if (intgradcut && intgradcut->coninfos_size) out_solver_log << "Added " << intgradcut->coninfos_size << " IntervalGradientCuts." << endl;
		}
		if (update_subdiv_bound(k_star, i_star, right)) {
			linear_relax->remove_node(right);
			right=NULL;
		}
	} else right->ref_point[i0]=cut;
	if (right && !boxreduce(right, i0, IntervalReduction::LOWER)) {
		linear_relax->remove_node(right);
		right=NULL;
	}
	if (right) nodes.push_back(right);
}

//----------------------------------------------------------------------

int MinlpBCP::bin_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node) {
	int fixed=0;
	for (map<int, set<int> >::iterator it(node->bcp_fixed_var.begin()); it!=node->bcp_fixed_var.end(); it++)
		fixed+=it->second.size();

	if (fixed==orig_prob->i_discr.size()) {
		out_solver_log << "All binaries fixed, so no more binary subdivision here." << endl;
		return 0;
	}

	if (bound_type==LP_bound /*|| bound_type==RMP_bound*/) { // (R) might have been improved
		int ret=set_LP_bound(node);
		if (ret==1) {
			out_solver_log << "(R[U]) infeasible" << endl;
			return 1;
		}
	}

	int k_star=-1; int i_star=-1; double max_cost=-1E-4; // to cope with almost feasible box constraints
	for (int k=0; k<split_prob->block.size(); k++) {
		dvector b_lag(linear_relax->obj->block[k].size());
		if (linear_relax->obj->b[k]) b_lag=*linear_relax->obj->b[k];
		if (node->dual_point.size()) {
			int c=0;
			for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->couple_con.begin()); it!=linear_relax->couple_con.end(); it++, c++)
				if (it->b[k])
					b_lag.AddMult(node->dual_point(c), *it->b[k]);
		}
		b_lag*=2;

		for (int i=0; i<split_prob->block[k].size(); i++) {
			int i0=split_prob->block[k][i];
			if (!split_prob->discr[i0]) continue; // no binary
			if (node->upper(i0)-node->lower(i0)<rtol) continue; // already fixed
			double cost=MIN(fabs(node->ref_point(i0)-node->lower(i0)), fabs(node->upper(i0)-node->ref_point(i0))) / (1+(node->upper(i0)-node->lower(i0)));
 			if (b_lag(i)>0) cost+=10*b_lag(i);
			if (cost>max_cost) {
				max_cost=cost;
				k_star=k;
				i_star=i;
			}
		}
	}
	if (k_star<0 || i_star<0) {
		out_log << node->ref_point;
		out_log << node->lower << node->upper;
		out_log << node->dual_point;
		out_log << k_star << ' ' << i_star << '\t' << max_cost << endl;
		for (int k=0; k<split_prob->block.size(); ++k)
			for (int i=0; i<split_prob->block[k].size(); ++i) {
				int i0=split_prob->block[k][i];
				out_log << i0 << ':' << split_prob->var_names[i0] << '\t' << node->lower(i0) << '\t' << node->upper(i0) << '\t' << node->ref_point(i0)
				<< '\t' << MIN(fabs(node->ref_point(i0)-node->lower(i0)), fabs(node->upper(i0)-node->ref_point(i0))) / (1+(node->upper(i0)-node->lower(i0)))
				 << endl;
			}
	}
//	out_log << "subdiv at " << k_star << " " << i_star << endl;
	assert(k_star>=0 && i_star>=0);
	subdiv_var=linear_relax->obj->block[k_star][i_star];
	rect_subdiv(nodes,node,k_star,i_star,node->lower[subdiv_var]);

	if (nodes.empty()) return 1;
	return 0;
}

//----------------------------------------------------------------------

void MinlpBCP::bisect_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node) {
	int istar,kstar;
	double dmax=-rtol;
	double cut;

	for (int k=0; k<split_prob->block.size(); k++) {
		for(int i=0; i<split_prob->block[k].size(); i++) {
			int i0=split_prob->block[k][i];
			double dist=split_prob->upper(i0)-split_prob->lower(i0);
			if (dist<rtol) continue;
			double dist_local=node->upper(i0)-node->lower(i0);
			if (dist_local<rtol) continue; // fixed variable in node
			double d=dist_local/dist;
			if(d>dmax) {
				dmax=d;
				kstar=k; istar=i;
				cut=node->lower(i0)+dist_local/2;
			}
		}
	}

	if (dmax>1E-4) {
		nr_subdiv_bisect++;
		rect_subdiv(nodes, node, kstar, istar, cut);
		subdiv_var=linear_relax->obj->block[kstar][istar];
	} // TODO: should the lower bound be fixed otherwise?
}

//----------------------------------------------------------------------

void MinlpBCP::cost_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node) {
	if (!lag_problem.size()) init_lag_problems(node);
	update_lag_problems(node->dual_point);

	int k_star=-1, i_star;
	int k_free=-1, i_free;
	double maxcost=-INFINITY;
	double maxcost_free=-INFINITY;
	double cut, cut_free;
	int i0_free;

	for (int k=0; k<lag_problem.size(); k++) {
		//project $x_{J_k}$ onto $G_k$ using lagsolve obtaining x
		dvector xest(lag_problem[k]->dim());
//---------------------------------------------------------------------
		if (subdiv_type==CostSubdivLag) {
			LagSolveStatus status;
			solve_lag_problem(status, k, node);
			if (status.ret) {
				out_solver_log << "Solving " << k << "th lag problem return: " << status.ret << endl;
				continue;
			}
			xest=status.solset.begin()->second;
		} else {
//-----------------------------------------------------------------------
//perform a Newton based estimation:
//$g_i(x_i+t_i\nabla g_i(x))=0$, hence  $t_i=-g_i(x)/\norm{\nabla g_i(x)}^2$
//set $\tilde x_{J_k}=\sum_{i\in I_{k,viol}} x_i/|I_{k,viol}|$
			dvector x_k(node->ref_point, linear_relax->obj->block[k]);
			dvector g_k(lag_problem[k]->dim());
//			xest=x_k;
			int viol=0;
			for(int i=0; i<lag_problem[k]->con.size(); i++) {
				double val=lag_problem[k]->con[i]->eval(x_k);
				if (val<1E-4 && (val>-1E-4 || !lag_problem[k]->con_eq[i])) continue; // not violated
				lag_problem[k]->con[i]->grad(g_k, x_k);
				double gradnorm=g_k.sq_norm2();
				if (gradnorm<rtol) continue; // stationary point of constraint, bad luck
				xest.AddMult(-val/gradnorm, g_k);
				viol++;
			}
			if (!viol) continue; // point feasible for lag problem
			xest*=1.0/((double)viol);
			xest+=x_k;
			Project::project(xest, xest, lag_problem[k]->lower, lag_problem[k]->upper);
		}
//-----------------------------------------------------------------------
		double pseudo_cost_block=lag_problem[k]->obj->eval(xest-node->ref_point(linear_relax->obj->block[k]));
//		out_log << "Pseudo costs block " << k << ": " << pseudo_cost_block << "\t " << lag_problem[k]->dim() << endl;
		if (fabs(pseudo_cost_block)<rtol) continue;

		for (int i=0; i<split_prob->block[k].size(); i++) {
			int i0=split_prob->block[k][i];
			double cost=2*(*lag_problem[k]->obj->b[0])(i)*(xest(i)-node->ref_point[i0]);///fabs(pseudo_cost_block);
//			double cost=fabs(xest(i)-node->ref_point[i0]);
//			out_log << "Pseudo costs " << lag_problem[k]->var_names[i] << ": " << cost << endl;
			if (cost<rtol) continue; // no or negative pseudo costs on this variable

			double diam=split_prob->upper(i0)-split_prob->lower(i0);
			if (node->upper(i0)-node->lower(i0)<=.01*diam) continue; // almost fixed variable
			double delta_i=xest(i)-node->lower(i0);
			if (delta_i > node->upper(i0)-xest(i))
				delta_i=node->upper(i0)-xest(i);
			delta_i/=diam;
			cost*=delta_i;

			if(cost>maxcost) {
				k_star=k;
				i_star=i;
				cut=xest(i);
				maxcost=cost;
			}
			if(!node->fix_branch_var.count(i0)&&(cost>maxcost_free)) {
				k_free=k;
				i_free=i;
				i0_free=i0;
				cut_free=xest(i);
				maxcost_free=cost;
			}
		}
	}
	if (k_star==-1) {
		out_solver_log << "No Cost subdivision possible." << endl;
		bisect_subdiv(nodes, subdiv_var, node);
		return;
	}
	out_solver_log <<  "pseudo costs: " << maxcost << "\t free: " << maxcost_free << endl;
	if(maxcost_free<=0.1*maxcost){
		out_solver_log << "Release " << node->fix_branch_var.size() << " branching-fixed variables." << endl;
		node->fix_branch_var.clear();
		k_free=k_star;
		i_free=i_star;
		cut_free=cut;
	}
	node->fix_branch_var.insert(i0_free);

	if (lag_problem[k_star]->discr[i_star]) { // discrete variable
		rect_subdiv(nodes, node, k_free, i_free, cut_free);
		subdiv_var=i0_free;
		return;
	}
	int i0_star=split_prob->block[k_star][i_star];
	double diam=node->upper(i0_star)-node->lower(i0_star);
	if (cut_free-node->lower(i0_star)<.05*diam)
		rect_subdiv(nodes, node, k_free, i_free, .75*node->lower(i0_star)+.25*node->upper(i0_star));
	else if (node->upper(i0_star)-cut_free<.05*diam)
		rect_subdiv(nodes, node, k_free, i_free, .25*node->lower(i0_star)+.75*node->upper(i0_star));
	else rect_subdiv(nodes, node, k_free,i_free,cut_free);
	subdiv_var=i0_free;
}

void MinlpBCP::viol_subdiv(list<Pointer<MinlpNode> >& nodes, int& subdiv_var, Pointer<MinlpNode> node) {
	int k_star=-1, i_star=-1, k2=-1, i2=-1;
	double val, maxviol=rtol, cut;

	dvector x(reform ? reform->get_short_vector(node->ref_point) : node->ref_point);
	assert(x.dim()==split_prob->dim());
	dvector g(x.dim());

	for (int c=0; c<split_prob->con.size(); c++) {
		val=split_prob->con[c]->eval(x);
		if (val<tol && (val>-tol || !split_prob->con_eq[c])) continue;

		split_prob->con[c]->grad(g, x);
		double gradnorm=sqrt(g.sq_norm2());
		if (gradnorm<rtol) continue;
		if (fabs(val)<maxviol*gradnorm) continue;

		maxviol=fabs(val/gradnorm);
		k_star=-1; i_star=-1; k2=-1; i2=-1;

		double maxdelta=0, g_min=INFINITY;
		for (int k=0; k<split_prob->block.size(); k++) {
			for (int i=0; i<split_prob->block[k].size(); i++) {
				int i0=split_prob->block[k][i];
				if (fabs(g(i0))<rtol) continue;
				double dist=split_prob->upper(i0)-split_prob->lower(i0);
				if (dist<rtol) continue; // fixed variable in problem
				double dist_local=node->upper(i0)-node->lower(i0);
				if (dist_local<rtol) continue; // fixed variable in node
				double delta_i=MIN(x(i0)-node->lower(i0), node->upper(i0)-x(i0))/dist;
				double sigma_i=dist_local/dist;
				if (sigma_i<1E-4) continue; // already very small (compared to beginning)
				if ((delta_i>0.2) && (fabs(g(i0))/sigma_i<g_min)) { // delta_i > .2 -> sigma_i>0
					g_min=fabs(g(i0))/sigma_i;
					k_star=k;
					i_star=i;
					cut=x(i0);
				}
				if(delta_i*sigma_i>maxdelta) {
					maxdelta=delta_i*sigma_i;
					k2=k;
					i2=i;
				}
			}
		}
		if (k_star==-1 && k2>=0) {
			int i0=split_prob->block[k2][i2];
			double dist=split_prob->upper(i0)-split_prob->lower(i0);
			double dist_local=node->upper(i0)-node->lower(i0);
			if (dist_local/dist>.01) { // cut only if diameter still more than 1% over original boxsize
				k_star=k2;
				i_star=i2;
				cut=x(i0);
				if (cut>node->upper(i0)-.2*dist_local)
					cut=node->upper(i0)-.2*dist_local;
				else if (cut<node->upper(i0)+.2*dist_local)
					cut=node->lower(i0)+.2*dist_local;
			}
		}
	}

	if (k_star==-1) {
		out_solver_log << "No Violation subdivision possible. Falling back to Bisection subdivision." << endl;
		bisect_subdiv(nodes, subdiv_var, node);
		return;
	}

	subdiv_var=linear_relax->obj->block[k_star][i_star];
//	out_solver_log << reform->ext_prob->var_names[subdiv_var] << ": max violation: " << maxviol << "\t dist_local: " << node->upper(subdiv_var)-node->lower(subdiv_var) << endl;
	rect_subdiv(nodes, node, k_star, i_star, cut);
}


/*
void MinlpBCP::cost_subdiv2(list<Pointer<MinlpNode> >& nodes, Pointer<MinlpNode> node) {
	if (!lag_problem.size()) init_lag_problems(node);
	update_lag_problems(node->dual_point);

	int ret;
	int k_star=-1, i_star;
	double max_pseudo_cost=-INFINITY;
	double max_midness;
	double cut;

	//for $k=1,\dots,p$:
	for (int k=0; k<split_prob->block.size(); k++) {
		//project $x_{J_k}$ onto $G_k$ using lagsolve obtaining x
		set<SolCandidate> sol_cand;
		ret=solve_lag_problem(sol_cand, k, node);
		if (ret) {
			out_solver_log << "Solving " << k << "th lag problem return: " << ret << endl;
			continue;
		}

		double pseudo_cost=lag_problem[k]->obj->eval(xest-node->ref_point(linear_relax->prob->block[k]));
		if (pseudo_cost+rtol<max_pseudo_cost) continue;

		max_midness=-INFINITY;
		for (int i=0; i<split_prob->block[k].size(); i++) {
			double diam=lag_problem[k]->upper(i)-lag_problem[k]->lower(i);
			if (diam<rtol) continue; // fixed variable
			double delta_i=xest(i)-lag_problem[k]->lower(i);
			if (delta_i > lag_problem[k]->upper(i)-xest(i))
				delta_i=lag_problem[k]->upper(i)-xest(i);
			delta_i/=diam;
			if (delta_i < max_midness+rtol) continue;
			max_midness=delta_i;
			k_star=k;
			i_star=i;
			cut=xest(i);
		}
		if (k_star!=-1) max_pseudo_cost=pseudo_cost;
	}
	if (k_star==-1) {
		out_solver_log << "No Cost subdivision possible." << endl;
		return;
	}

	max_midness=MIN(cut-lag_problem[k_star]->lower[i_star], lag_problem[k_star]->upper[i_star]-cut);
	max_midness/=lag_problem[k_star]->upper[i_star]-lag_problem[k_star]->lower[i_star];
	out_solver_log <<  "pseudo costs: " << max_pseudo_cost << "\t midness: " << max_midness << endl;

	if (max_midness<tol) {
		bisect_subdiv(nodes, node);
		return;
	}

	rect_subdiv(nodes, node, k_star,i_star,cut);
}
*/


//-------------------------------------------------------------
// Solving Lagrange problems
//-------------------------------------------------------------

void MinlpBCP::init_lag_problems(Pointer<MinlpNode> node) {
	MinlpProblem& prob(reform ? *reform->ext_prob : *split_prob);
	lag_problem.clear();
	lag_problem.reserve(prob.block.size());
	block_sub_convex_prob.clear();
	block_sub_convex_prob.reserve(prob.block.size());
	
	if (!block_prob.size()) init_block_problems();

	for (int k=0; k<prob.block.size(); k++) {
		lag_problem.push_back(new MinlpProblem(*block_prob[k]));
		block_sub_convex_prob.push_back(new MinlpProblem(*block_convex_prob[k]));

		lag_problem[k]->add_obj(new SepQcFunc(NULL, new SparseVector<double>(lag_problem[k]->dim()), NULL));
// new MinlpPartLagFunc(opt.linear_prob, dual_point, k, false)));

		// partition cuts from this node
		if (node->part_con.size()) {
			for (list<Pointer<SepQcFunc> >::iterator it(node->part_con[k].begin()); it!=node->part_con[k].end(); it++) {
				lag_problem[k]->add_con(new SepQcFunc((*it)->A[k], (*it)->b[k], (*it)->s[k], (*it)->c), false, "part con");
				block_sub_convex_prob[k]->add_con(lag_problem[k]->con.back(), lag_problem[k]->con_eq.back(), lag_problem[k]->con_names.back());
			}
		}

		// local cuts from this node for (C_k)
// todo: add local cuts for this node again to block convex relax

		for (int i=0; i<prob.block[k].size(); ++i) {
			block_sub_convex_prob[k]->lower[i]=lag_problem[k]->lower[i]=node->lower(prob.block[k][i]);
			block_sub_convex_prob[k]->upper[i]=lag_problem[k]->upper[i]=node->upper(prob.block[k][i]);
		}
/*
		block_sub_convex_prob.push_back(new MinlpProblem(*lag_problem[k]));
		param->add("Relax bounds update", "0");
		param->add("Relax decomposition", "0");
		Pointer<ostream> temp_log(out_log_p), temp_out(out_out_p);
		out_log_p=NULL; out_out_p=NULL;
		Relax(param).relax(*block_sub_convex_prob.back());
		out_log_p=temp_log; out_out_p=temp_out;
		block_sub_convex_prob.back()->i_discr.clear();
		block_sub_convex_prob.back()->i_cont.clear();
		for (int i=0; i<block_sub_convex_prob.back()->dim(); i++) {
			block_sub_convex_prob.back()->i_cont.push_back(i);
			block_sub_convex_prob.back()->discr[i]=false;
		}
*/	}

	for (int k=0; k<block_sub_convex_prob.size(); k++)
		block_sub_convex_prob[k]->obj=lag_problem[k]->obj;
}

//----------------------------------------------------------------------

void MinlpBCP::update_lag_problems(const dvector& dual_point) {
	for (int k=0; k<lag_problem.size(); k++) {
		if (linear_relax->obj->b[k]) *lag_problem[k]->obj->b[0]=*linear_relax->obj->b[k];
		else *lag_problem[k]->obj->b[0]=0.;
		int c=0;
		for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->couple_con.begin()); it!=linear_relax->couple_con.end(); it++, c++)
			if (it->b[k]) lag_problem[k]->obj->b[0]->AddMult(dual_point(c), *it->b[k]);
	}
}

//----------------------------------------------------------------------

void MinlpBCP::update_lag_problem(int k, const dvector& a) {
	*lag_problem[k]->obj->b[0]=a;
	*lag_problem[k]->obj->b[0]*=.5;
}

//----------------------------------------------------------------------

void MinlpBCP::solve_lag_problem(LagSolveStatus& status, int k, Pointer<MinlpNode> node, Pointer<SepQcFunc> temp_cut) {
	assert(lag_problem.size()==linear_relax->obj->block.size() && lag_problem[k]);
Timer t;
	Pointer<LinearRelax> block_sub_linear_relax(new LinearRelax(*linear_relax, k, node, lag_problem[k]->obj));

	Pointer<dvector> diam(new dvector(*sol_cand_diam, linear_relax->obj->block[k]));

	int temp_cut_lag, temp_cut_lin, temp_cut_conv; // indices of temporary cut in block problems

	if (temp_cut) {
		temp_cut_lag=lag_problem[k]->con.size();
		lag_problem[k]->add_con(temp_cut, false, "temp cut");
//		temp_cut_lin=block_sub_linear_relax[k]->prob->con.size();
		temp_cut_conv=block_sub_convex_prob[k]->con.size();
//		block_sub_linear_relax[k]->prob->add_con(temp_cut, false, "temp cut");
		block_sub_linear_relax->add_cut(new SimpleCut(temp_cut->b[0], temp_cut->c), 0);
		block_sub_convex_prob[k]->add_con(temp_cut, false, "temp cut");
	}

//	out_log << "Solving " << k << "th lag problem, binaries: " << lag_problem[k]->i_discr.size() << ", dim: " << lag_problem[k]->dim() << endl;
//	if (lag_problem[k]->i_discr.size()<lag_problem[k]->dim())
//	out_log << *lag_problem[k];

	Pointer<RelaxationSolver> heu;
	switch(lagsolve_type) {
		case BranchCut : {
			Pointer<MinlpBCP> bcp=new MinlpBCP(lag_problem[k], lag_problem[k], block_sub_linear_relax, false, sol_cand_closeval_tol, diam, param, NULL/*out_solver_log_p*/, NULL);
			bcp->set_convex_prob(block_sub_convex_prob[k]);
			bcp->pre_bound_type=NLP_bound;
			bcp->maj_bound_type=NLP_bound;
			bcp->pre_bb_max_iter=0;
//			bcp->pre_bound_type=NLP_bound;
//			bcp->maj_bound_type=stop_bound;
			bcp->pre_bb_max_iter=bcp->iter_max;
			bcp->max_time=max_time-timer->stop();
			bcp->upper_bound_effort_level=0;
//			if (subdiv_type==CostSubdivLag) bcp->subdiv_type=CostSubdivNewton;
			bcp->subdiv_type=BinSubdiv;
			if (temp_cut) bcp->is_maxcut=false;
			bcp->intervalreduction=NULL; // since we do not have correct sparsity info in the lag problem
			heu=bcp;
		} break;

		default: {
			out_err << "MinlpNode::solve_lag_problem: MinlpBCP::lagsolve_type unknown. Aborting." << endl;
			exit(-1);
		}
	}
	if ((!temp_cut) && node->lagprob_solutions.size() && node->lagprob_solutions[k].dim()) status.ret=heu->solve(node->lagprob_solutions[k]);
	else status.ret=heu->solve();

//	out_log << ' ' << ((MinlpBCP*)(RelaxationSolver*)heu)->find_solcand_time/t.stop() << ' ';
/*if (status.ret) {
	heu->out_solver_log_p=out_log_p;
	heu->out_solver_p=out_log_p;
	heu->solve();
//out_log << *lag_problem[k];

//exit(0);
}*/
	status.iter=heu->iter();
	status.solset=heu->sol_cand;
/*	if (!status.ret) status.lowbound=heu->opt_val();
	else */status.lowbound=heu->lower_bound;
	status.value=heu->opt_val();
	if (status.lowbound==-INFINITY && status.value<INFINITY) status.lowbound=status.value; // not nice
	heu=NULL;

	if (temp_cut) {
		lag_problem[k]->del_con(temp_cut_lag);
		block_sub_convex_prob[k]->del_con(temp_cut_conv);
	} else if (status.solset.size()) {
		if (!node->lagprob_solutions.size()) node->lagprob_solutions.resize(lag_problem.size());
		node->lagprob_solutions[k]=status.solset.begin()->second;
		// keep global cuts from the subproblem as local cuts
//		out_log << "From " << linear_relax->nr_local_cuts(node);
		linear_relax->integrate(*block_sub_linear_relax, k, node);
//		out_log << " to " << linear_relax->nr_local_cuts(node) << endl;
	}

	lagprob_solves++;
}

int MinlpBCP::conv_rate_check(double val) {
  if (val>=conv_rate_cntrl_last_val-rtol) {  // serious step ?
	  conv_rate_cntrl_improve_iter++;
		conv_rate_cntrl_last_val=val;  // update last serious value
	}	else return 0;  // check only non-null-steps

	if (conv_rate_cntrl_improve_iter==0) { // first iteration
		conv_rate_cntrl_last_major_val=conv_rate_cntrl_first_major_val=val;
		return 0;
	}
	if (conv_rate_cntrl_improve_iter%conv_rate_cntrl_minor_iter == 0) { // major iteration
		if (conv_rate_cntrl_improve_iter==conv_rate_cntrl_minor_iter) { // first major iteration
			conv_rate_cntrl_max_rel_improvement=(val-conv_rate_cntrl_first_major_val)/(fabs(conv_rate_cntrl_first_major_val)+1);
			conv_rate_cntrl_last_major_val=val;
			return 0;
		}
		double conv_rate_cntrl_last_rel_improvement=(val-conv_rate_cntrl_last_major_val)/(fabs(conv_rate_cntrl_first_major_val)+1);  // compute relative improvment
		conv_rate_cntrl_last_major_val=val;
		if (conv_rate_cntrl_last_rel_improvement > conv_rate_cntrl_max_rel_improvement) {  // update max improvement
	   	conv_rate_cntrl_max_rel_improvement=conv_rate_cntrl_last_rel_improvement;
	    return 0;
	  }
		if (conv_rate_cntrl_last_rel_improvement <= conv_rate_cntrl_max_rel_improvement*conv_rate_cntrl_stopping_rho) {
			out_solver_log << "Convergence too slow: " << conv_rate_cntrl_last_rel_improvement << " <= " << conv_rate_cntrl_max_rel_improvement*conv_rate_cntrl_stopping_rho << endl;
			return 10;
		}
	}

  return 0;
}

void MinlpBCP::mem_check() {
	if (!mem_limit) return;

	unsigned int mem=get_mem();
	mem/=1024;
	out_solver_log << "Mem used: " << mem << "KB";

	mem/=1024;
	if (mem<mem_limit) {
		out_solver_log << endl;
		return;
	}
	out_solver_log << '\t';

	int toremove=bb_tree.size()/10;
	out_solver << "Memory limit reached. Removing " << toremove << " nodes.";

	multimap<double, Pointer<MinlpNode> >::iterator it;
	while (toremove--) {
		it=--bb_tree.end();
		linear_relax->remove_node(it->second);
		bb_tree.erase(it);
	}
/*
	multimap<double, Pointer<MinlpNode> >::reverse_iterator it(bb_tree.rbegin());
	for (int i=0; it!=bb_tree.rend() && i<100; ++i)
		linear_relax->remove_node(it->second);
*/
	out_solver << "\t Mem used now: " << get_mem()/1024 << "KB" << endl;
/*
	while (mem>.9*mem_limit && bb_tree.size()>10) {
		Pointer<MinlpNode> node=it->second;
		linear_relax->remove_node(node);
		node=NULL;
		bb_tree.erase(--bb_tree.end());
		it=bb_tree.rbegin();
		mem=get_mem(procfile, dummy)/1024;
		out_solver_log << mem << ' ';
		mem/=1024;
//		++it;
	}
	out_solver_log << endl;
*/
}

// ------------------------------------------- main loop ----------------------------------------------------

int MinlpBCP::solve(dvector& start) {
	Timer t;

	Pointer<MinlpProblem> prob(reform ? reform->ext_prob : split_prob);
	assert(start.dim()==prob->dim());

	dvector trialpoint(start.dim());
	Round::round(trialpoint, start, split_prob->i_discr, prob->lower, prob->upper);

	int ret=locopt_NLP(trialpoint).first;
	out_solver_log << "Local minimization for starting point: " << ret << endl;
	find_solcand_time+=t.stop();

	return solve();
}

int MinlpBCP::solve() {
	Pointer<MinlpNode> node1;
	int ret;
	iter_=0;
//	timer->start();
	double final_gap;

	map<int, set<int> > fixed;
	map<int, set<int> >::iterator it_block;
	for (int k=0; k<split_prob->block.size(); k++) {
		it_block=fixed.end();
		for (int i=0; i<split_prob->block[k].size(); i++) {
			int i0=split_prob->block[k][i];
			if (split_prob->discr[i0] && split_prob->upper(i0)-split_prob->lower(i0)<rtol) {
				if (it_block==fixed.end()) {
					set<int> newset;
					it_block=fixed.insert(pair<int, set<int> >(k, newset)).first;
				}
				it_block->second.insert(i);
				out_solver_log << "Recognized variable " << split_prob->var_names[i0] << " as already fixed to " << split_prob->lower[i0] << endl;
			}
		}
	}

	bound_type = pre_bb_max_iter ? pre_bound_type : maj_bound_type;
	ret=0;
	node1=new MinlpNode(reform ? reform->ext_prob->lower : split_prob->lower, reform ? reform->ext_prob->upper : split_prob->upper);
	node1->bcp_fixed_var=fixed;
	out_solver_log << "Computing lower bound of initial node. Current bound: " << node1->low_bound << endl;
	if (sol_C && sol_C_is_solution) node1->low_bound=reform ? reform->ext_convex_prob->obj->eval(*sol_C) : convex_prob->obj->eval(*sol_C);
	if (bound_type==RMP_bound || bound_type==LP_RMP_bound)
		if (init_ExtremePoints(node1))
			out_solver << "Initialization of Extreme points failed." << endl;
	if (bound_type==NLP_bound) {
		if (!(sol_C && sol_C_is_solution)) ret=set_low_bound(node1);
	} else
		ret=set_low_bound(node1) && (!(sol_C && sol_C_is_solution));
//	linear_relax->generate_cuts(NULL);
	if (ret) {
		out_solver << "Could not compute lower bound of root node." << endl;
		return 1;
	}
	out_solver_log << endl;

	bb_tree.insert(pair<double, Pointer<MinlpNode> >(node1->key(), node1));
	node1=NULL;
	final_gap=start_bb();

	double bcp_time=timer->stop();
	out_solver << "MinlpBCP iterations: " << iter() << endl;
	out_solver << "MinlpBCP last improvement iteration: " << last_impr_iter << endl;
	out_solver << "Final number of cuts: " << linear_relax->nr_all_cuts();
	if (linear_relax->cutlimit_reached()) out_solver << "(l)";
//	out_solver << "\t Global: " << linear_relax->nr_global_cuts();
	out_solver << endl;
	out_solver << "MinlpBCP final bb-tree size: " << bb_tree.size() << endl;
	if (bound_computed) out_solver << "Boundcomputation failures: " << (100*bound_failed)/bound_computed << "\\%" << endl;
	if (!lag_problem.empty()) out_solver << "Lagrangian subproblems solves: " << lagprob_solves << endl;
	if (bcp_time>1E-4) out_solver << "Bounds time: " << ((int)(1000*bound_time/bcp_time))/10. << "\\%" << endl;
	if (bound_time>1E-4 && colgen) out_solver << "Init RMP time: " << ((int)(1000*init_RMP_time/bound_time))/10. << "\\%" << endl;
	if (find_solcand_time<rtol) find_solcand_time=0.;
	if (bcp_time>1E-4) out_solver << "find_sol_candidates time: " << ((int)(1000*find_solcand_time/bcp_time))/10. << "\\%" << endl;
	if (subdiv_time<rtol) subdiv_time=0.;
	if (bcp_time>1E-4) out_solver << "Subdivision time: " << ((int)(1000*subdiv_time/bcp_time))/10. << "\\%" << endl;
	if (nr_subdiv_contvar) out_solver << "Subdivisions on cont. var by bisection: " << (100*nr_subdiv_bisect)/nr_subdiv_contvar << "\\%" << endl;
	if (point_conversion_warnings) out_solver << "Point conversion warnings: " << point_conversion_warnings << endl;
	if (last_impr_iter>=0) {
		if (fabs(final_gap)<.01) {
			out_solver << "Final Gap: $<1$\\%" << endl;
		} else {
			out_solver.precision(final_gap<.1 ? 2 : (final_gap<1. ? 3 : 4));
			out_solver << "Final Gap: " << 100.*final_gap << "\\%" << endl;
			out_solver.precision(6);
		}
	}
	if (colgen) {
		int ExtremePoints_nr=0; for (int k=0; k<ExtremePoints.size(); k++) ExtremePoints_nr+=ExtremePoints[k].size();
		out_solver << "Final number of extreme points: " << ExtremePoints_nr << endl;
	}
	
//	while (!bb_tree.empty()) {
//		multimap<double, Pointer<MinlpNode> >::iterator it(bb_tree.begin());
//		out_log << "Pruning node with " << linear_relax->nr_local_cuts(it->second) << " local cuts -> ";
//		linear_relax->remove_node(it->second);
//		bb_tree.erase(it);
//		out_log << linear_relax->nr_all_cuts() << "\t global: " << linear_relax->nr_global_cuts() << endl;		
//	}

	return (!sol_cand.size());
}
//----------------------------------------------------------------------

/* Start branch and bound loop */
double MinlpBCP::start_bb() {
	if (iter_>=iter_max) return 1.;

	Pointer<MinlpNode> node1;
	double gap=1.;
	int ret=0;

	int pre_bb_iter=0;
	if (pre_bb_max_iter) {
		out_solver << "Starting preprocessing" << endl;
		bound_type=pre_bound_type;
	}

	do {
		iter_++;
		out_solver << "BCP Iter. " << iter_ << " / " << iter_max << "\t Time: " << timer->stop();
		if (out_solver_p && max_time<INFINITY) { out_solver << " / "; Timer::print(*out_solver_p, max_time); }
		out_solver << "\t tree: " << bb_tree.size() << "\t cuts: " << linear_relax->nr_all_cuts();
		if (out_solver_p && linear_relax->cutlimit_reached()) out_solver << "(l)";
		out_solver << "\t lower: " << bb_tree.begin()->first << "\t upper: " << opt_val_ << "\t " << nr_solcand_found << "\t gap: " << gap << endl;
//		max_bb_tree_size=MAX(max_bb_tree_size, bb_tree.size());

		mem_check();

		pair<bool, double> bound_impr;
		do {
			// (selection) Take a partition element $U$ from $L$.
			/* take node from the tree */
			node1=bb_tree.begin()->second;
			bb_tree.erase(bb_tree.begin());
			clean_sub_problems();

			if (out_solver_log_p && node1->bcp_fixed_var.size()) {
				out_solver_log << "Fixation: ";
				for (map<int, set<int> >::iterator it_block(node1->bcp_fixed_var.begin()); it_block!=node1->bcp_fixed_var.end(); it_block++)
					for (set<int>::iterator it_var(it_block->second.begin()); it_var!=it_block->second.end(); it_var++)
						out_solver_log << split_prob->var_names[split_prob->block[it_block->first][*it_var]] << "(" << node1->lower(linear_relax->obj->block[it_block->first][*it_var]) << ") ";
				out_solver_log << endl;
			}
			out_solver_log << "Box diameter (2-norm): " << sqrt((node1->upper-node1->lower).sq_norm2());
/*			if (out_solver_log_p && node1->box_cuts.size()) {
				out_solver_log << "Box cuts: ";
				for (map<int, map<int, pair<double, double> > >::iterator it_block(node1->box_cuts.begin()); it_block!=node1->box_cuts.end(); it_block++)
					for (map<int, pair<double, double> >::iterator it_var(it_block->second.begin()); it_var!=it_block->second.end(); it_var++)
						out_solver_log << split_prob->var_names[split_prob->block[it_block->first][it_var->first]] << "(" << it_var->second.first << ", " << it_var->second.second << ") ";
				out_solver_log << endl;
			}
*/			out_solver_log << "\t Local cuts: " << linear_relax->nr_local_cuts(node1) << endl;

			bound_impr=improve_bound(node1);
		} while (((!bound_impr.first) || (node1->low_bound>=opt_val()-rtol)) && bb_tree.size());
		if ((!bound_impr.first) || (node1->low_bound>=opt_val()-rtol)) { // all nodes infeasible or their lower bound exceeds the upper bound
			if (lower_bound==-INFINITY) lower_bound=opt_val();
			out_solver_log << endl;
			break;
		}

		if (bound_impr.second) out_solver_log << "Bound improved by " << bound_impr.second << endl;

		//(solution search)\\
		//If $U\cap X_{\rm cand}=\emptyset$,
		set<SolCandidate>::const_iterator sol_candidate_it(node1->outside_part_set(sol_cand));
		if (sol_candidate_it==sol_cand.end() || bound_type==stop_bound) {
			ret=find_sol_candidates(node1);
			out_solver_log << "find_sol_candidates return: " << ret << endl;
		}

		if (bound_type!=stop_bound && opt_val()>node1->low_bound+rtol) { // some room left for improvement
			if (bound_impr.second>bound_impr_tol) {
				bb_tree.insert(pair<double, Pointer<MinlpNode> >(node1->key(orig_prob->i_discr.size()), node1));
			} else {
				list<Pointer<MinlpNode> > nodes;
				int subdiv_var=subdivide(nodes, node1);
				linear_relax->remove_node(node1);
				node1=NULL;

				for (list<Pointer<MinlpNode> >::iterator it(nodes.begin()); it!=nodes.end(); it++) {
					if (!(*it)->update_subdiv_bound_called) {
						out_solver_log << "Compute lower bound of new node." << endl;
						clean_sub_problems();
						ret=set_low_bound(*it);
					} else ret=0;
					if (!ret) bb_tree.insert(pair<double, Pointer<MinlpNode> >((*it)->key(orig_prob->i_discr.size()), *it));
					else {
						linear_relax->remove_node(*it);
						out_solver_log << "Node not added to bb-tree." << endl;
					}
				}
				// no subdivided node added, and no other node left
				if ((!bb_tree.size()) && (lower_bound==-INFINITY)) lower_bound=opt_val();
			}
		}

		/** Prune Branch- and Bound-Tree */
		multimap<double, Pointer<MinlpNode> >::iterator prunestart(bb_tree.lower_bound(opt_val_));
		if (prunestart!=bb_tree.end()) {
			int nrpruned=0;
			for (multimap<double, Pointer<MinlpNode> >::iterator it(prunestart); it!=bb_tree.end(); ++it, ++nrpruned)
				linear_relax->remove_node(it->second);
			out_solver_log << "Pruning " << nrpruned << " nodes." << endl;
			bb_tree.erase(prunestart, bb_tree.end());
		}

		if (bb_tree.size()) gap=(opt_val_-bb_tree.begin()->first)/(1+fabs(bb_tree.begin()->first));
		else {
			gap=0.;
			if (lower_bound==-INFINITY) lower_bound=node1 ? node1->low_bound : opt_val();
		}
		prune_ExtremePoints(bb_tree);

		if (bound_print && opt_val_<INFINITY) {
			*bound_print << (double)timer->stop() << "\t ";
			if (lower_bound==-INFINITY) *bound_print << bb_tree.begin()->first;
			else *bound_print << lower_bound;
			*bound_print << "\t " << opt_val_ << "\t ";
			if (lower_bound==-INFINITY) *bound_print << gap << endl;
			else *bound_print << (opt_val_-lower_bound)/(1+fabs(lower_bound)) << endl;
		}

		out_solver_log << endl;
		if (pre_bb_max_iter) {
			if (iter()>=pre_bb_max_iter || timer->stop()>.8*max_time || conv_rate_check(bb_tree.begin()->first)) {
				out_solver << "Finished preprocessing. Starting with next bound type." << endl;
				out_solver << "MinlpBCP preprocessing iterations: " << iter() << endl;
				bound_type=maj_bound_type;
				pre_bb_max_iter=0;
				if (maj_bound_type==stop_bound && lower_bound==-INFINITY && bb_tree.size())
					lower_bound=bb_tree.begin()->first;
			}
		}

		if (node1) linear_relax->remove_node(node1);
	} while (bb_tree.size() && (gap>gap_tol) && iter_<iter_max && (timer->stop()<max_time));
	if (pre_bb_max_iter) out_solver << "MinlpBCP preprocessing iterations: " << iter() << endl;
	if (lower_bound==-INFINITY && bb_tree.size()) lower_bound=bb_tree.begin()->first;
	out_solver_log << "Final lower bound: " << lower_bound << endl;

	if (last_impr_iter>=0) return (opt_val_-lower_bound)/(1+fabs(lower_bound));
	return 1.;
}
//----------------------------------------------------------------------

bool MinlpBCP::add_sol_candidate(const dvector& x) {
	bool already_known=RelaxationSolver::add_sol_candidate(x);

	if (already_known) return true;
	
	if (!linear_relax || linear_relax->cutlimit_reached()) return false;

	Pointer<MinlpProblem> prob;
	Pointer<MinlpProblem> conv;
	if (reform && (x.dim()==reform->ext_prob->dim())) { // if we use the extended version
		prob=reform->ext_prob;
		conv=reform->ext_convex_prob;
	} else {
		prob=split_prob ? split_prob : orig_prob;
		conv=convex_prob;
	}
	dvector& lower(current_node ? current_node->lower : prob->lower);
	dvector& upper(current_node ? current_node->upper : prob->upper);

	if (conv) { // adding linearization cuts
		list<pair<LinearizationCut, pair<int, bool> > > cuts;
		
		linconcutgen.get_cuts(cuts, x, lower, upper);
		
		int local_cuts_nr=0, global_cuts_nr=0;
		for (list<pair<LinearizationCut, pair<int, bool> > >::iterator cutit(cuts.begin()); cutit!=cuts.end() && !linear_relax->cutlimit_reached(); ++cutit) {
			linear_relax->add_cut(Pointer<LinearizationCut>(new LinearizationCut(cutit->first)), cutit->second.first, cutit->second.second ? NULL : current_node);
			if (cutit->second.second) ++global_cuts_nr; else ++local_cuts_nr;
		}
			
		out_log << "Added " << global_cuts_nr << " global and " << local_cuts_nr << " local LinearizationCuts." << endl;
	}

	if (intgrad_cuts) { // adding IntervalGradientCuts
		int cuts_nr=0;
		for (int k=0; k<prob->block.size() && !linear_relax->cutlimit_reached(); k++) {
			Pointer<IntervalGradientCut> intgradcut(intgrad_cutgen.get_cuts(x(prob->block[k]), k, lower(prob->block[k]), upper(prob->block[k])));
			if (intgradcut) {
				cuts_nr+=intgradcut->coninfos_size;
				linear_relax->add_cut(intgradcut, k, current_node);
			}
		}
		if (cuts_nr) out_solver_log << "Added " << cuts_nr << " IntervalGradientCuts." << endl;
	}

	return already_known;
}

int MinlpBCP::find_sol_candidates(Pointer<MinlpNode> node) {
	Timer t;

	current_node=node;

	// to get uniform rounding for copied variables (if x and y satisfy x==y, they should do so after rounding as well)
	for (int i=0; i<split_prob->i_discr.size(); i++) {
		int i0=split_prob->i_discr[i];
		double mid=.5*(split_prob->lower(i0)+split_prob->upper(i0));
		if (fabs(node->ref_point(i0)-mid)<rtol) node->ref_point[i0]=mid;
	}
	dvector trialpoint(node->ref_point.dim());
	Round::round(trialpoint, node->ref_point, split_prob->i_discr, node->lower, node->upper);

	int ret;
	// LocOpt from rounded solution of relaxation
	if (is_maxcut && (split_prob->block.size()==1)) {
		assert(split_prob->con.size());
		trialpoint[split_prob->dim()-1]=0.;
		trialpoint[split_prob->dim()-1]=split_prob->con[0]->eval(trialpoint);
		if (!split_prob->feasible(trialpoint, tol, NULL)) {
			add_sol_candidate(trialpoint);
			ret=0;
		} else ret=1;
		find_solcand_time+=t.stop();
		return ret;
	} else {
		ret=locopt_NLP(trialpoint).first;
		out_solver_log << "locopt NLP return: " << ret << endl;
		if (!ret) nr_solcand_found++;
	}

	// Preswitching
	if (upper_bound_effort_level>=1 && (!is_maxcut)) {
		int ret2=preswitching(node);
		if (ret && (!ret2)) ret=0;
	}

	// Lag Heu
	if (node->i_ExtremePoints.size()) {
		char* name=NULL;
		if ((!lagheu) && (name=param->get("LagHeu"))) {
			if (!strcmp(name, "first")) lagheu=new LagHeu1(orig_prob, linear_relax, is_gams_prob, sol_cand_closeval_tol, sol_cand_diam, param, out_solver_p, out_solver_log_p);
			else if (!strcmp(name, "second")) lagheu=new LagHeu2(orig_prob, linear_relax, is_gams_prob, sol_cand_closeval_tol, sol_cand_diam, param, out_solver_p, out_solver_log_p);
			else if (!strcmp(name, "second b")) lagheu=new LagHeu2b(orig_prob, linear_relax, is_gams_prob, sol_cand_closeval_tol, sol_cand_diam, param, out_solver_p, out_solver_log_p);
			else lagheu=new LagHeu_SimAnnealing(orig_prob, linear_relax, is_gams_prob, sol_cand_closeval_tol, sol_cand_diam, param, out_solver_p, out_solver_log_p);
			if (split_prob) lagheu->set_split_prob(split_prob);
			if (reform) lagheu->set_reform(reform);
		}
		if (lagheu) {
			int ret2=lagheu->solve(node);
			if (!ret2) {
				if (lagheu->opt_val()<opt_val_) {
					opt_val_=lagheu->opt_val();
					sol_point=lagheu->sol_point;
					last_impr_iter=iter();
				}
				sol_cand.insert(lagheu->sol_cand.begin(), lagheu->sol_cand.end()); // this will also try to add already added points!
				ret=0;
			}
		}
	}

	find_solcand_time+=t.stop();
	return ret;
}

int MinlpBCP::preswitching(Pointer<MinlpNode> node) {
	int midbinaries=0; // binaries, which are undecided yet
 	for (int i=0; i<split_prob->i_discr.size(); i++) {
  	int i0=split_prob->i_discr[i];
		if (split_prob->upper[i0]-split_prob->lower[i0]<rtol) continue;
		double mid=MIN(node->ref_point[i0]-split_prob->lower[i0], split_prob->upper[i0]-node->ref_point[i0])/(split_prob->upper[i0]-split_prob->lower[i0]);
		if (mid<.02) continue;
		midbinaries++;
	}
	dvector trialpoint(node->ref_point.dim());
	int ret=1;
	int samplesize=midbinaries;//(int)pow(log(midbinaries+1.), 2);
	out_solver_log << "Apply random switching " << samplesize << " times for " << midbinaries << " undecided binaries." << endl;
	while (samplesize--) {
		Round::round(trialpoint, node->ref_point, split_prob->i_discr, node->lower, node->upper);
		bool changed=false;
  	for (int i=0; i<split_prob->i_discr.size(); i++) {
	  	int i0=split_prob->i_discr[i];
			if (split_prob->upper[i0]-split_prob->lower[i0]<rtol) continue;
			double mid=MIN(node->ref_point[i0]-split_prob->lower[i0], split_prob->upper[i0]-node->ref_point[i0])/(split_prob->upper[i0]-split_prob->lower[i0]);
			if (mid<.02) continue;
			if (random(0., 1.)<mid) {
				trialpoint[i0]=split_prob->lower[i0]+split_prob->upper[i0]-trialpoint[i0]; // switch variable
				changed=true;
			}
		}
		if (changed) {
			int ret2;
			if (is_maxcut && split_prob->block.size()==1) {
				trialpoint[split_prob->dim()-1]=0.;
				trialpoint[split_prob->dim()-1]=split_prob->con[0]->eval(trialpoint);
				if (!split_prob->feasible(trialpoint, tol, NULL)) {
					add_sol_candidate(trialpoint);
					ret2=0;
				} else ret2=1;
			} else {
//			out_solver_log << "trying " << trialpoint;
				ret2=locopt_NLP(trialpoint).first;
				out_solver_log << ret2;
				if (!ret2) nr_solcand_found++;
			}

			if (ret && (!ret2)) ret=0;
			if (ret && (ret2!=1)) ret=ret2;
		}
	}
	out_solver_log << endl;

	return ret;
}
