// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "boxfind.h"

pair<int,int> BoundsFinder::compute_bounds(MinlpProblem& conv_prob, vector<bool>& discr) {
	pair<int,int> ret;
	switch (method) {
		case 0: ret=compute_bounds_guess(conv_prob);
			break;
		case 1: {
			dvector lower(conv_prob.lower), upper(conv_prob.upper);
			ret=compute_bounds_expensive(conv_prob, lower, upper, discr);
			conv_prob.lower=lower;
			conv_prob.upper=upper;
		} break;
		case 2: {
			dvector lower(conv_prob.lower), upper(conv_prob.upper);
			ret=compute_bounds_expensive2(conv_prob, lower, upper, discr);
			conv_prob.lower=lower;
			conv_prob.upper=upper;
		} break;
		default: {
			out_err << "BoundsFinder: method not known." << endl;
			exit(-1);
		}
	}
	if (ret.second) out_out << "Warning: " << ret.second << " variable bounds couldn't be determined and were guessed." << endl;
	return ret;
}

double BoundsFinder::compute_bound(MinlpProblem& conv_prob, int& ret, int index, int block, bool low) {
	conv_prob.obj->b[block]->SetElement(index, low ? .5 : -.5);

	Pointer<LocOpt> locopt=LocOpt::get_solver(Pointer<MinlpProblem>(&conv_prob, false), param, "BoundsFinder", NULL, NULL);
	ret=locopt->solve(); //locopt->solve(conv_prob.primal_point);
	if (ret) out_log << "BoundsFinder " << (low ? "low" : " up") << ": "
		<< (conv_prob.var_names[conv_prob.block[block][index]] ? (char*)conv_prob.var_names[conv_prob.block[block][index]] : "noname")
		<< " LocOpt: " << ret << " val=" << (low ? 1 : -1)*locopt->opt_val() << endl;
//	if (ret && conv_prob.feasible(locopt->sol_point, 1E-4, NULL)) ret=1; // solution point not feasible

	*conv_prob.obj->b[block]=0;
	return (low ? 1 : -1) * locopt->opt_val();
}

pair<int,int> BoundsFinder::compute_bounds_guess(MinlpProblem& conv_prob) {
	int guessed=0;
	bool missing_bounds=false;
	double min_lower=-1.;
	double max_upper=1.;
	for (int i=0; i<conv_prob.dim(); i++) {
		if (conv_prob.lower[i]>-INFINITY) {
			if (conv_prob.lower[i]<min_lower) min_lower=conv_prob.lower[i];
		} else missing_bounds=true;
		if (conv_prob.upper[i]<INFINITY) {
			if (conv_prob.upper[i]>max_upper) max_upper=conv_prob.upper[i];
		} else missing_bounds=true;
	}
	if (missing_bounds)
		for (int i=0; i<conv_prob.dim(); i++) {
			if (conv_prob.lower[i]<=-INFINITY || conv_prob.upper[i]>=INFINITY) ++guessed;
			if (conv_prob.lower[i]<=-INFINITY && low) conv_prob.lower[i]=10*min_lower;
			if (conv_prob.upper[i]>=INFINITY && up) conv_prob.upper[i]=10*max_upper;
		}

  return pair<int,int>(0, guessed);
}

pair<int,int> BoundsFinder::compute_bounds_expensive(MinlpProblem& conv_prob, dvector& new_lower, dvector& new_upper, vector<bool>& discr) {
	bool skip_binaries=param->get_i("Box Reduction skip binaries", 0);

	int ret=0;
	int guessed=0;
	bool missing_bounds=false;
	bool all_discr=true;
	double min_lower=-1.;
	double max_upper=1.;
	for (int i=0; i<conv_prob.dim(); i++) {
		if (conv_prob.lower[i]>-INFINITY) {
			if (conv_prob.lower[i]<min_lower) min_lower=conv_prob.lower[i];
		} else missing_bounds=true;
		if (conv_prob.upper[i]<INFINITY) {
			if (conv_prob.upper[i]>max_upper) max_upper=conv_prob.upper[i];
		} else missing_bounds=true;
		if (discr[i]) all_discr=false;
	}
	if (skip_binaries && all_discr) return pair<int,int>(0,0);

	Pointer<SepQcFunc> orig_obj(conv_prob.obj);
	conv_prob.obj=new SepQcFunc(conv_prob.block);
	for (int k=0; k<conv_prob.block.size(); k++)
		conv_prob.obj->b[k]=new SparseVector<double>(conv_prob.block[k].size());
	(*conv_prob.obj->b[0])[0]=1.;  // to get correct number of nonzeros while construction of SnOpt

	int locopt_ret;

	*conv_prob.obj->b[0]=0.;

	for (int k=0; k<conv_prob.block.size(); k++)
		for (int i=0; i<conv_prob.block[k].size(); i++) {
			bool guessed_lower=false;
			int index=conv_prob.block[k][i];
			if (skip_binaries && discr[index]) continue;
			if (conv_prob.upper[index]-conv_prob.lower[index]<rtol) continue; // fixed variable

			if (low && (known || conv_prob.lower[index]<=-INFINITY)) {
				double newlow=compute_bound(conv_prob, locopt_ret, i, k, true);
				if (locopt_ret==0) { // good return from SnOpt
					if (conv_prob.lower[index]<=-INFINITY) new_lower[index]=conv_prob.lower[index]=newlow;
					else if (newlow > conv_prob.lower[index]+rtol) {  // improved bound
						if (newlow>conv_prob.upper[index]) newlow=conv_prob.upper[index];
						if (discr[index]) newlow=newlow>conv_prob.lower[index]+1E-4 ? conv_prob.upper[index] : conv_prob.lower[index];
						new_lower[index]=conv_prob.lower[index]=newlow;
					}
				} else {
					if (conv_prob.lower[index]<=-INFINITY) {
						++guessed;
						guessed_lower=true;
						if (newlow>-1E+10) new_lower[index]=MIN(10*min_lower, newlow);
						else new_lower[index]=10*min_lower;
					}
					++ret;
				}
			}

			if (up && (known || conv_prob.upper[index]>=INFINITY)) {
				double newup=compute_bound(conv_prob, locopt_ret, i, k, false);
				if (locopt_ret==0) { // good return from SnOpt
					if (conv_prob.upper[index]>=INFINITY) new_upper[index]=conv_prob.upper[index]=newup;
					else if (newup < conv_prob.upper[index]-rtol) { // improved bound
						if (newup<conv_prob.lower[index]) newup=conv_prob.lower[index];
						if (discr[index]) newup=(newup<conv_prob.upper[index]-1E-4) ? conv_prob.lower[index] : conv_prob.upper[index];
						new_upper[index]=conv_prob.upper[index]=newup;
					}
				} else {
					if (conv_prob.upper[index]>=INFINITY) {
						if (!guessed_lower) ++guessed; // do not count twice
						if (newup<1E+10) new_upper[index]=MAX(newup, 10*max_upper);
						else newup=10*max_upper;
					}
					++ret;
				}
			}
		}

	conv_prob.obj = orig_obj;

	return pair<int,int>(ret, guessed);
}


pair<int,int> BoundsFinder::compute_bounds_expensive2(MinlpProblem& conv_prob, dvector& lower, dvector& upper, vector<bool>& discr) {
	int i=0;
	int old_errors; // problems in prior run
	pair<int,int> new_errors(2*conv_prob.dim(), 0);  // problems in last run
	do {
		old_errors=new_errors.first;
		new_errors=compute_bounds_expensive(conv_prob, lower, upper, discr);
		out_log << "BoundsFinder: problems in run number " << i++ << ": " << new_errors.first << endl;
	} while (new_errors.first && new_errors.first<old_errors); // while we improve the number of bounds, where we are sure

	return new_errors;
}

void BoundsFinder::compute_linbounds(MinlpProblem& prob) {
	for (int c=0; c<prob.con.size(); c++) {
		bool linear=true;
		map<int, double> coeff;
		for (int k=0; linear && k<prob.block.size(); k++) {
			if (prob.con[c]->A[k] || prob.con[c]->s[k]) {
				linear=false;
				break;
			}
			if (!prob.con[c]->b[k]) continue;
			UserVector<double>& b(*prob.con[c]->b[k]);
			for (int i=0; i<b.dim(); i++) {
				double bi=b(i);
				if (bi) coeff.insert(pair<int, double>(prob.block[k][i], 2*bi));
			}
		}
		if (!linear) continue;

		for (map<int,double>::iterator it(coeff.begin()); it!=coeff.end(); it++) {
			if ((!known) && prob.lower[it->first]>-INFINITY && prob.upper[it->first]<INFINITY) continue;
			bool min_known=true, max_known=true;
			double min=0., max=0.;
			for (map<int,double>::iterator it2(coeff.begin()); it2!=coeff.end(); it2++) {
				if (it==it2) continue;
				if (it2->second>0) {
					if (prob.lower[it2->first]>-INFINITY) min+=it2->second*prob.lower[it2->first];
					else min_known=false;
					if (prob.upper[it2->first]<INFINITY) max+=it2->second*prob.upper[it2->first];
					else max_known=false;
				} else {
					if (prob.upper[it2->first]<INFINITY) min+=it2->second*prob.upper[it2->first];
					else min_known=false;
					if (prob.lower[it2->first]>-INFINITY) max+=it2->second*prob.lower[it2->first];
					else max_known=false;
				}
			}
			if ((!min_known) && (!max_known)) continue;
			min+=prob.con[c]->c; max+=prob.con[c]->c;
			min/=-it->second; max/=-it->second;
			if (it->second>0) {
				if (min_known && prob.upper[it->first]>min+rtol && (known || prob.upper[it->first]>=INFINITY)) {
					if (prob.discr[it->first]) min=prob.lower[it->first];
//					out_log << "Constraint " << prob.con_names[c] << ": Reduced  upper bound of " << prob.var_names[it->first] << " from " << prob.upper[it->first] << " to " << min << endl;
					prob.upper[it->first]=min;
				}
				if (max_known && prob.con_eq[c] && prob.lower[it->first]<max-rtol && (known || prob.lower[it->first]<=-INFINITY)) {
					if (prob.discr[it->first]) max=prob.upper[it->first];
//					out_log << "Constraint " << prob.con_names[c] << ": Heighten lower bound of " << prob.var_names[it->first] << " from " << prob.lower[it->first] << " to " << max << endl;
					prob.lower[it->first]=max;
				}
			} else {
				if (min_known && prob.lower[it->first]<min-rtol && (known || prob.lower[it->first]<=-INFINITY)) {
					if (prob.discr[it->first]) min=prob.upper[it->first];
//					out_log << "Constraint " << prob.con_names[c] << ": Heighten lower bound of " << prob.var_names[it->first] << " from " << prob.lower[it->first] << " to " << min << endl;
					prob.lower[it->first]=min;
				}
				if (max_known && prob.con_eq[c] && prob.upper[it->first]>max+rtol && (known || prob.upper[it->first]>=INFINITY)) {
					if (prob.discr[it->first]) max=prob.lower[it->first];
//					out_log << "Constraint " << prob.con_names[c] << ": Reduced  upper bound of " << prob.var_names[it->first] << " from " << prob.upper[it->first] << " to " << max << endl;
					prob.upper[it->first]=max;
				}
			}
			if (prob.lower[it->first]>prob.upper[it->first]) {
				out_log << "Warning: Upper bound < Lower bound, variable " << prob.var_names[it->first] << " " << prob.lower[it->first] << " " << prob.upper[it->first] << endl;
				prob.lower[it->first]=prob.upper[it->first]=.5*(prob.lower[it->first]+prob.upper[it->first]);
			}
		}
	}
}

// -------------------------------- BoundsFinderLinear -------------------------

BoundsFinderLinear::BoundsFinderLinear(Pointer<MinlpProblem> prob_, Pointer<Param> param_)
: param(param_), prob(prob_), low(true), up(true), known(true)
{ MipProblem mip(prob->lower, prob->upper, vector<int>(0), prob->var_names);
	mip.resize_row(prob->con.size());
	int lincons=0;
	for (int c=prob->con.size()-1; c>=0; --c) {
		if (prob->con[c]->sparsity_available() && prob->con[c]->get_sparsity().nonlinear->size()) continue;
		bool linear=true;
		for (int k=prob->block.size()-1; linear && k>=0; --k) linear=(!prob->con[c]->A[k]) && (!prob->con[c]->s[k]);
		if (!linear) continue;
		SparseVector<double> row(prob->con[c]->b, prob->block);
		row*=2;
		mip.setRow(lincons++, row, prob->con_eq[c] ? -prob->con[c]->c : -INFINITY, -prob->con[c]->c, prob->con_names[c]);
	}
	if (!lincons) return;
	mip.resize_row(lincons);

	mip.setObj(new SparseVector<double>(mip.dim()));
	mip.finish();

	solver=MIPSolver::get_solver(mip, param);
//	solver->set_maxiter(10*mip.rows());
}

MIPSolver::SolutionStatus BoundsFinderLinear::compute_bound(double& bound, int index, bool low) {
	solver->set_obj(SparseVector<double>(solver->nr_col(), index, low ? 1. : -1.));

	MIPSolver::SolutionStatus ret(solver->solve());
	if (ret==MIPSolver::SOLVED || ret==MIPSolver::FEASIBLE) bound=(low ? 1. : -1.)*solver->get_optval();

	return ret;
}

int BoundsFinderLinear::compute(dvector& newlow, dvector& newup, set<pair<int, IntervalReduction::which_bound_type> >* changed_var, double min_impr) {
	newlow=prob->lower;
	newup=prob->upper;
	if (!solver) return 0; // cannot reduce any bounds if there are no constraints

	int nr_problems=0;
	double newbound;

	for (int i=0; i<prob->dim(); ++i) {
		if (prob->upper(i)-prob->lower(i)<rtol) continue; // fixed variable

		if (low && (known || prob->lower(i)<=-INFINITY)) {
			MIPSolver::SolutionStatus ret=compute_bound(newbound, i, true);
			switch(ret) {
				case MIPSolver::SOLVED: newlow[i]=MIN(newup(i), MAX(prob->lower(i), newbound)); // should be newbound always
					break;
				case MIPSolver::UNBOUNDED: break; // ignore
				default: ++nr_problems;
					out_err << "BoundsFinderLinear lower of " << prob->var_names[i] << " returned " << ret << endl;
					break;
			}
		}
		if (up && (known || prob->upper(i)>=INFINITY)) {
			MIPSolver::SolutionStatus ret=compute_bound(newbound, i, false);
			switch(ret) {
				case MIPSolver::SOLVED: newup[i]=MAX(newlow(i), MIN(prob->upper(i), newbound)); // should be newbound always
					break;
				case MIPSolver::UNBOUNDED: break; // ignore
				default: ++nr_problems;
					out_err << "BoundsFinderLinear upper of " << prob->var_names[i] << " returned " << ret << endl;
					break;
			}
		}

//		if (newlow(i)>prob->lower(i)+1E-4 || newup(i)<prob->upper(i)-1E-4) out_log << "BoundsFinderLinear reduce box of " << prob->var_names[i]
//			<< " from [" << prob->lower(i) << ", " << prob->upper(i) << "] to ["  << newlow(i) << ", " << newup(i) << ']' << endl;

		if (prob->discr[i]) { // maybe fixing binary
			if (newlow(i)>prob->lower(i)+1E-4) newlow[i]=newup[i]=prob->upper(i); //newup[i];
			else {
				newlow[i]=prob->lower(i);
				if (newup(i)<prob->upper(i)-1E-4) newup[i]=newlow[i]=prob->lower(i); //newlow[i];
				else newup[i]=prob->upper(i);
			}
		}

		if (changed_var) { // if we want specific information about the variables, which bound was changed
			if (prob->lower(i)<=-INFINITY)
				if (newlow(i)>-INFINITY) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i, IntervalReduction::LOWER));
				else ;
			else
				if (prob->upper(i)>=INFINITY)
					if (newlow(i)-prob->lower(i)>rtol) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i, IntervalReduction::LOWER));
					else ;
				else
					if (newlow(i)-prob->lower(i)>min_impr*(prob->upper(i)-prob->lower(i))) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i, IntervalReduction::LOWER));
					else ;
			if (prob->upper(i)>=INFINITY)
				if (newup(i)<INFINITY) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i, IntervalReduction::UPPER));
				else ;
			else
				if (prob->lower(i)<=-INFINITY)
					if (prob->upper(i)-newup(i)>rtol) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i, IntervalReduction::UPPER));
					else ;
				else
					if (prob->upper(i)-newup(i)>min_impr*(prob->upper(i)-prob->lower(i))) changed_var->insert(pair<int, IntervalReduction::which_bound_type>(i, IntervalReduction::UPPER));
					else ;
		}
	}

	return nr_problems;
}

// --------------------------------- DualBounds --------------------------------

DualBounds::DualBounds(Pointer<MinlpProblem> S_, Pointer<MinlpProblem> C_, Pointer<Param> param_, int type_)
: S(S_), C(C_), param(param_), type(type_)
{	// constructing R
	R=new MinlpProblem(*C);
	dvector mid(S->lower); mid+=S->upper; mid*=.5;
	R->taylor_approx(mid, 1);

	// constructing partial lagrangian problems (without objective)
	Pointer<MinlpProblem> prob(type==1 ? C : S);
	lin_con2.resize(prob->con.size(), true);
	for (int k=0; k<S->block.size(); k++) {
		lag_prob.push_back(new MinlpProblem(new SepQcFunc(S->block[k].size()), S->lower(S->block[k]), S->upper(S->block[k]), false));
		if (type==2) { // define binaries
			lag_prob.back()->i_cont.clear();
			for (int i=0; i<S->block[k].size(); i++)
				if (S->discr[S->block[k][i]]) {
					lag_prob.back()->discr[i]=true;
					lag_prob.back()->i_discr.push_back(i);
				} else lag_prob.back()->i_cont.push_back(i);
		}
		lin_con.push_back(list<int>(0));

		for (int c=0; c<prob->con.size(); c++)
			if (prob->con[c]->A[k] || prob->con[c]->s[k]) { // adding nonlinear constraints of this block
				lag_prob.back()->add_con(new SepQcFunc(prob->con[c]->A[k], prob->con[c]->b[k], prob->con[c]->s[k], prob->con[c]->c),
					prob->con_eq[c], prob->con_names[c]);
				lin_con2[c]=false;
			}
			else if (prob->con[c]->b[k]) {
				lin_con.back().push_back(c);
				bool is_block_con=true;
				for (int k_=0; is_block_con && k_<S->block.size(); k_++) {
					if (k_==k) continue;
					is_block_con=!prob->con[c]->b[k_];
				}
				if (is_block_con) lag_prob.back()->add_con(new SepQcFunc(NULL, prob->con[c]->b[k], NULL, prob->con[c]->c),
					prob->con_eq[c], prob->con_names[c]);
			}
	}
}

double DualBounds::dual_bound(UserVector<double>& a) {
	for (int k=0; k<R->block.size(); k++) R->obj->b[k]=a.getcopy(R->block[k]);
	Pointer<LocOpt> R_solver(LocOpt::get_lp_solver(R, param, NULL, NULL, NULL));
	int ret=R_solver->solve();
	if (ret) out_log << endl << "Solving (R[a]) returned " << ret << endl;
	if (ret) return -INFINITY;
	dvector mu(R_solver->get_lag_multipliers());
	if (mu==0) out_log << "(=0) ";
//	out_log << "mu: " << mu;

	double dual_bound=0;
	Pointer<MinlpProblem> prob(type==1 ? C : S);
	for (int k=0; k<lag_prob.size(); k++) {
		// construct objective L_k
		lag_prob[k]->obj->b[0]=new SparseVector<double>(a, S->block[k]);
//		*lag_prob[k]->obj->b[0]*=.5;
		for (list<int>::iterator it(lin_con[k].begin()); it!=lin_con[k].end(); it++)
			if (fabs(mu(*it))>rtol) lag_prob[k]->obj->b[0]->AddMult(mu(*it), *prob->con[*it]->b[k]);

		if (*lag_prob[k]->obj->b[0]==0) continue;

			Pointer<LocOpt> locopt=LocOpt::get_solver(lag_prob[k], param, "DualBounds", NULL, NULL);
			dvector start(R_solver->sol_point, prob->block[k]);
			int ret=locopt->solve(start);
			if (ret) out_log << endl << "Solving partial lag prob, block " << k << " return: " << ret << "\t value: " << locopt->opt_val() << endl;
			// SnOpt failed -> returning
			if (ret!=0 && ret!=9) return -INFINITY;
			if (ret && lag_prob[k]->feasible(locopt->sol_point, 1E-4, NULL)) return -INFINITY;
			dual_bound+=locopt->opt_val(); // ret==0 or (ret==9 and feasibility checked)

	}
	for (int c=0; c<S->con.size(); c++)
		if (lin_con2[c]) dual_bound+=mu(c)*prob->con[c]->c;

//	out_log << "Dual bound: " << dual_bound << " for objective: ";
//	if (out_log_p) R->obj->print(*out_log_p, R->var_names);

	return dual_bound;
}

int DualBounds::update_box(int n) {
	int problems=0;
	int skip_binaries=param->get_i("Box Reduction skip binaries", 0);
	if (skip_binaries && S->i_discr.size()==S->dim()) return 0; // binary-problem
	SparseVector<double> ei(S->dim());
	for (int i=0; i<n; i++) {
		if (skip_binaries && S->discr[i]) continue; // skipping binary variable
		out_log << "Updating box for " << S->var_names[i] << ": lower ";
		double old_low=S->lower[i], old_up=S->upper[i];
		ei.SetElement(i, .5); // lower bound
		double bound=dual_bound(ei);
		if (bound>-INFINITY && bound>C->lower[i]+rtol)
			if (S->discr[i]) S->lower[i]=C->lower[i]=C->upper[i]; // fixing binary
			else S->lower[i]=C->lower[i]=MIN(bound, C->upper[i]);
		if (bound==-INFINITY) problems++;

		out_log << "and upper ";
		ei*=-1; // upper bound
		bound=-dual_bound(ei);
		if (bound<INFINITY && bound<C->upper[i]-rtol)
			if (S->discr[i]) S->upper[i]=C->upper[i]=C->lower[i];
			else C->upper[i]=S->upper[i]=MAX(bound, C->lower[i]);
		ei.DelElement(i);
		if (bound==INFINITY) problems++;

		out_log << "reduced from [" << old_low << ", " << old_up << "] to [" << S->lower[i] << ", " << S->upper[i] << "]" << endl;
	}

	return problems;
}

double DualBounds::obj_bound() {
	dvector b(S->obj->b, S->block);
	double bnd=dual_bound(b);
	if (bnd>-INFINITY) return bnd+S->obj->c;
	return -INFINITY;
}


// ---------------------------- IntervalReduction --------------------------------

ostream& operator<<(ostream& out, const IntervalReduction::EdgeData& ed) {
	switch (ed.which_bound) {
		case IntervalReduction::LOWER: out << "low"; break;
		case IntervalReduction::UPPER: out << "up"; break;
		case IntervalReduction::WHATEVER: out << "both"; break;
	}
	out << ',' << ed.con_nr;
	return out;
};

const IntervalReduction::EdgeData& IntervalReduction::EdgeData::operator+=(const EdgeData& ed) const {
	out_err << "Adding EdgeData to const EdgeData not allowed. Aborting." << endl;
	exit(-1);
	return *this;
}

void IntervalReduction::set_problem(Pointer<MinlpProblem> prob_) {
	prob=prob_;
	reduction_by_block.resize(prob->block.size());

	dependency_graph.clear();
	vector<set<DependencyGraph::NodeType>::iterator> nodes(prob->dim()+1);
	int i0;
	for (int k=prob->block.size()-1; k>=0; --k)
		for (int i=prob->block[k].size()-1; i>=0; --i) {
			int i0=prob->block[k](i);
			nodes[i0]=dependency_graph.add_node(i0, NodeData(k,i));
		}
	set<DependencyGraph::NodeType>::iterator dummy_node(dependency_graph.add_node(prob->dim()));

	for (int c=0; c<prob->con.size(); ++c) {
		for (SepQcFunc::VariableIterator it(*prob->con[c]); it; ++it) { // iterate over all variables
			bool single=true;
			for (SepQcFunc::VariableIterator it_lin(*prob->con[c], true, false, false); it_lin; ++it_lin) { // iterate over linear variables
				if (it()==it_lin()) continue;

				single=false;
				which_bound_type wb;
				if (prob->con_eq[c]) wb=WHATEVER;
				else if (it.type()&VariableIterator_Type::NONLINEAR) wb=WHATEVER; // it is nonlinear -> cannot predict
				else wb=it.coeff_lin()*it.coeff_lin()>0 ? LOWER : UPPER; // it is linear

				dependency_graph.add_edge(nodes[it()], nodes[it_lin()], EdgeData(c, wb, it_lin.coeff_lin()));
			}
			if (single && it.type()==VariableIterator_Type::LINEAR) // constraint of form a*x + c <=/= 0 -> add dummy edge from dummy node to this variable
				dependency_graph.add_edge(dummy_node, nodes[it()], EdgeData(c, WHATEVER, it.coeff_lin()));
		}
	}

	out_log << "IntervalReduction dependency graph nr. edges: " << dependency_graph.m();
	if (dummy_node->get_adj().size()) out_log << ' ' << '(' << dummy_node->get_adj().size() << ')';
	out_log << endl;

//	out_log << "Dependency graph: " << dependency_graph;
}

void IntervalReduction::run(dvector& newlow, dvector& newup, const dvector& oldlow, const dvector& oldup, set<pair<const DependencyGraph::NodeType*, which_bound_type> >& nodeset) {
	int funceval=0;
	const int maxfunceval=10000;
	fixed_binaries.clear();
	empty_boxes=false;
	bool box_changed=false;
#ifdef FILIB_AVAILABLE
	IntervalVector box(oldlow, oldup);
	for (int i=0; i<box.dim(); ++i) {
		if (box(i).inf()==-INFINITY) box[i]=interval<double>(filib::fp_traits<double>::ninfinity(), box(i).sup());
		if (box(i).sup()==INFINITY) box[i]=interval<double>(box(i).inf(), filib::fp_traits<double>::infinity());
	}

	interval<double> oldbounds, newbounds, val;
	while (!nodeset.empty() && funceval<maxfunceval) {
		const DependencyGraph::NodeType& node(*nodeset.begin()->first);
		which_bound_type wb=nodeset.begin()->second;
		nodeset.erase(nodeset.begin());
		if (!nodeset.empty() && nodeset.begin()->first->idx()==node.idx()) { // same node two time in nodeset, so both bounds have changed
			nodeset.erase(nodeset.begin());
			wb=WHATEVER;
		}

//		out_log << "picked node " << node.idx() << " with bound " << wb << endl;
		const DependencyGraph::NodeType::adj_type& adj(node.get_adj());
		for (DependencyGraph::NodeType::adj_type::const_iterator it(adj.begin()); it!=adj.end() && funceval<maxfunceval; ++it) { // check each adjacent edge
			if (it->second->data.which_bound==LOWER && wb==UPPER) continue;
			if (it->second->data.which_bound==UPPER && wb==LOWER) continue;

			const DependencyGraph::EdgeType& edge(*it->second);
//			out_log << "taking edge to " << it->first << ": " << edge.data << endl;

			unsigned int index=edge.get_node2().idx();
			int con_nr=edge.data.con_nr;
			double coeff=edge.data.coeff;

			oldbounds=box(index);
			if (oldbounds.diam()<rtol) continue;

			box[index]=interval<double>::ZERO();
			val=-prob->con[con_nr]->eval(box)/coeff;
			++funceval;
			
			if (prob->con_eq[con_nr]) {
				newbounds=oldbounds.intersect(val);
				if (newbounds.isEmpty()) { // maybe rounding error, adding some tolerance
					val+=interval<double>(-rtol, rtol);
					newbounds=oldbounds.intersect(val);
					if (!newbounds.isEmpty()) newbounds=.5*(newbounds.inf()+newbounds.sup());
				}
			}
			else newbounds=coeff<0
				?	interval<double>(MAX(oldbounds.inf(), val.inf()), oldbounds.sup())
				:	interval<double>(oldbounds.inf(), MIN(oldbounds.sup(), val.sup()));

			if (prob->discr[index]) {
				if (newbounds.inf()>oldbounds.inf()+rtol) {
					newbounds=interval<double>(oldbounds.sup());
					const NodeData& node2(edge.get_node2().data);
					assert(node2.block_nr>=0 && node2.block_index>=0);
					fixed_binaries.insert(pair<int,int>(node2.block_nr, node2.block_index));
				} else if (newbounds.sup()<oldbounds.sup()-rtol) {
					newbounds=interval<double>(oldbounds.inf());
					const NodeData& node2(edge.get_node2().data);
					assert(node2.block_nr>=0 && node2.block_index>=0);
					fixed_binaries.insert(pair<int,int>(node2.block_nr, node2.block_index));
				} else newbounds=oldbounds; // no rtol's in bounds of binary variables please
			}

			box[index]=newbounds;

			if (oldbounds==newbounds) continue; // no change of box
			box_changed=true;

			if (do_print) out_log << prob->con_names[con_nr] << ": Reduce box of " << prob->var_names[index] << " from " << oldbounds << " to " << newbounds << "\t (" << val << ')';
			if (newbounds.isEmpty()) {
				empty_boxes=true;
				if (do_print) out_log << endl;
				break;
			}

			bool low_changed, up_changed;
			if (filib::fp_traits<double>::IsInf(oldbounds.inf())) low_changed=!filib::fp_traits<double>::IsInf(newbounds.inf());
			else if (filib::fp_traits<double>::IsInf(oldbounds.sup())) low_changed=newbounds.inf()-oldbounds.inf()>rtol;
			else low_changed=newbounds.inf()-oldbounds.inf()>min_impr*(oldbounds.diam()+1);
			if (filib::fp_traits<double>::IsInf(oldbounds.sup())) up_changed=!filib::fp_traits<double>::IsInf(newbounds.sup());
			else if (filib::fp_traits<double>::IsInf(oldbounds.inf())) up_changed=oldbounds.sup()-newbounds.sup()>rtol;
			else up_changed=oldbounds.sup()-newbounds.sup()>min_impr*(oldbounds.diam()+1);
			if (do_print) out_log << "\t low changed: " << low_changed << "\t up changed: " << up_changed << endl;

			if (low_changed || up_changed) nodeset.insert(
				pair<const DependencyGraph::NodeType*, which_bound_type>(&edge.get_node2(),
					low_changed && up_changed ? WHATEVER : (low_changed ? LOWER : UPPER) ) );
		}
	}

	out_log << "IntervalReduction: " << funceval << " function evaluations";
	if (!box_changed) {
		out_log << "\t no reduction" << endl;
		newlow=oldlow;
		newup=oldup;
		reduction_by_block=1.;
		reduction=1.;
		return;
	}

	if (!fixed_binaries.empty()) out_log << "\t fixed binaries: " << fixed_binaries.size();
	if (empty_boxes) {
		out_log << "\t empty boxes!" << endl;
		reduction=-INFINITY;
		return;
	}

	int i0;
	double olddiam_var, olddiam=0., newdiam=0.;
	reduction=0.;
	for (int k=prob->block.size()-1; k>=0; --k) {
		reduction_by_block[k]=0.;
		for (int i=prob->block[k].size()-1; i>=0; --i) {
			i0=prob->block[k][i];
			olddiam_var=oldup(i0)-oldlow(i0);
			if (box(i0).inf()==filib::fp_traits<double>::ninfinity()) newlow[i0]=-INFINITY;
			else newlow[i0]=box(i0).inf();
			if (box(i0).sup()==filib::fp_traits<double>::infinity()) newup[i0]=INFINITY;
			else newup[i0]=box(i0).sup();
			reduction_by_block[k]+=olddiam_var ? box(i0).diam()/olddiam_var : 1.;
		}
		reduction+=reduction_by_block[k];
		reduction_by_block[k]/=prob->block[k].size();
	}
	reduction/=prob->dim();

	out_log << "\t reduction: " << reduction;
	out_log << endl;
#else
	newlow=oldlow;
	newup=oldup;
	reduction=1.;
	reduction_by_block=1.;
#endif
}

void IntervalReduction::compute(dvector& newlow, dvector& newup, const dvector& oldlow, const dvector& oldup) {
	set<pair<const DependencyGraph::NodeType*, which_bound_type> > nodeset;
	for (set<DependencyGraph::NodeType>::iterator it(dependency_graph.nodes.begin()); it!=dependency_graph.nodes.end(); ++it)
		nodeset.insert(pair<const DependencyGraph::NodeType*, which_bound_type>(&*it, WHATEVER));

	run(newlow, newup, oldlow, oldup, nodeset);
}

void IntervalReduction::compute(dvector& newlow, dvector& newup, const dvector& oldlow, const dvector& oldup, set<pair<int, which_bound_type> >& startset) {
	if (startset.empty()) {
		out_log << "IntervalReduction: empty startset" << endl;
		return;
	}
	set<pair<const DependencyGraph::NodeType*, which_bound_type> > nodeset;
	for (set<pair<int, which_bound_type> >::iterator it(startset.begin()); it!=startset.end(); ++it)
		nodeset.insert(pair<const DependencyGraph::NodeType*, which_bound_type>(&dependency_graph.get_node(it->first), it->second));

	run(newlow, newup, oldlow, oldup, nodeset);
}

void IntervalReduction::print_small_boxes(dvector& low, dvector& up) {
	multimap<double, int> box;
	for (int i=0; i<low.dim(); ++i)
		box.insert(pair<double,int>((up(i)-low(i))/(1+fabs(low(i))+fabs(up(i))), i));

	out_log << "small boxes: " << endl;
	for (multimap<double,int>::iterator it(box.begin()); it!=box.end() && it->first<1e-4; ++it)
		out_log << prob->var_names[it->second] << ": " << it->first << "\t (" << low(it->second) << ',' << up(it->second) << ')' << endl;
}
