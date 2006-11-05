// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "ipopt2.h"
#ifdef IPOPT_AVAILABLE

IpOptProblem::IpOptProblem(const Pointer<MinlpProblem> prob_, IpOpt& ipopt_)
: prob(prob_), ipopt(ipopt_), nnz_jac_g(0), nnz_h_lag(0)
{
	for (int c=prob->con.size()-1; c>=0; --c)
		for (int k=prob->block.size()-1; k>=0; --k)
			if (prob->con[c]->b[k] || prob->con[c]->A[k] || prob->con[c]->s[k])
				if (prob->con[c]->sparsity_available(k)) nnz_jac_g+=prob->con[c]->get_sparsity(k).nr_var();
				else nnz_jac_g+=prob->block[k].size();

	prob->get_sparsity(sparsity);
	for (int k=prob->block.size()-1; k>=0; --k) {
		if (sparsity[k] && sparsity[k]->sparsity_pattern) {
			nnz_h_lag+=sparsity[k]->sparsity_pattern->size();
			nnz_h_lag+=sparsity[k]->nr_var(false);
		}
		else nnz_h_lag+=prob->block[k].size()*prob->block[k].size();
	}
}


bool IpOptProblem::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g_, Index& nnz_h_lag_, IndexStyleEnum &index_style) {
	n=prob->dim();
	m=prob->con.size();
	nnz_jac_g_=nnz_jac_g;
	nnz_h_lag_=nnz_h_lag;
	index_style=C_STYLE;

	return true;
}


bool IpOptProblem::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
	for (int i=n-1; i>=0; --i)
		if (prob->discr[i]) {
			x_l[i]=x_u[i]=ipopt.sol_point(i);
		}
		else {
			x_l[i]=prob->lower(i);
			x_u[i]=prob->upper(i);
		}

	for (int j=m-1; j>=0; --j) {
		g_l[j]=prob->con_eq[j] ? 0. : -INFINITY;
		g_u[j]=0.;
	}

	return true;
}


bool IpOptProblem::get_starting_point(Index n, bool init_x, Number* x,
			       bool init_z, Number* z_L, Number* z_U,
			       Index m, bool init_lambda, Number* lambda) {
	assert(init_z == false);
	assert(init_lambda == false);
	if (!init_x) return true; // no primal starting point needed

	// we initialize x in bounds, in the upper right quadrant
	for (int i=n-1; i>=0; --i) x[i] = ipopt.sol_point[i];

	return true;
}

bool IpOptProblem::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
	obj_value=prob->obj->eval(dvector(x,n));
	if (isnan(obj_value)) return false;
  return true;
}

bool IpOptProblem::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
	dvector xx(x, n);

	dvector g(n);
	prob->obj->grad(g, xx);
	int i=n; grad_f+=n;
	while (--i>=0) *(--grad_f)=g(i);

  return true;
}

bool IpOptProblem::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
	dvector xx(x, n);

	for (int c=m-1; c>=0; --c) {
		g[c]=prob->con[c]->eval(xx);
//		out_log << c << ' ' << prob->con_names[c] << ": " << g[c] << endl;
	}

  return true;
}

bool IpOptProblem::eval_jac_g(Index n, const Number* x, bool new_x,
		       Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values) {
	int count=0;
  if (values == NULL) { // return the structure of the jacobian of the constraints
		for (int k=0; k<prob->block.size(); ++k)
			for (int c=0; c<prob->con.size(); ++c)
				if (prob->con[c]->b[k] || prob->con[c]->A[k] || prob->con[c]->s[k]) {
					if (prob->con[c]->sparsity_available(k))
						for (VariableIterator it(prob->con[c]->get_sparsity(k)); it; ++it, ++count) {
							*(jCol++)=prob->block[k](it());
							*(iRow++)=c;
						}
					else
						for (int i=0; i<prob->block[k].size(); ++i, ++count) {
							*(jCol++)=prob->block[k](i);
							*(iRow++)=c;
						}
				}
		assert(count==nele_jac);

	} else {
		dvector* g=NULL;
		for (int k=0; k<prob->block.size(); ++k) {
			dvector xx(prob->block[k].size());
			for (int i=xx.dim()-1; i>=0; --i) xx[i]=x[prob->block[k][i]];
			for (int c=0; c<m; ++c) {
				SepQcFunc& func(*prob->con[c]);
				if (func.A[k] || func.s[k]) {
					if (!g) g=new dvector(xx.dim());
					func.grad(*g, xx, k);
					if (func.sparsity_available(k)) {
						for (VariableIterator it(func.get_sparsity(k)); it; ++it/*, ++count*/)
							*(values++)=(*g)(it());
					} else
						for (int i=0; i<func.block[k].size(); ++i/*, ++count*/)
							*(values++)=(*g)(i);
				} else if(func.b[k]) {
					if (func.sparsity_available(k)) {
						for (VariableIterator it(func.get_sparsity(k)); it; ++it/*, ++count*/) {
							*(values++)=it.coeff_lin();
						}
					} else {
						UserVector<double>& b(*func.b[k]);
						for (int i=0; i<func.block[k].size(); ++i/*, ++count*/)
							*(values++)=b(i);
					}
				}
			}
			if (g) { delete g; g=NULL; }
		}
//		assert(count==nele_jac);
	}

  return true;
}

template <bool add>
void IpOptProblem::set_hessian(const SepQcFunc& func, int blocknr, Number* values, double factor, dvector& xx, dvector& y, dvector& z) {
	const map<pair<int,int>, SparsityInfo::NonlinearConnection>& spall(*sparsity[blocknr]->sparsity_pattern);
	assert(func.sparsity_available(blocknr));
	const map<pair<int,int>, SparsityInfo::NonlinearConnection>& spfunc(*func.get_sparsity(blocknr).sparsity_pattern);

	map<pair<int,int>, SparsityInfo::NonlinearConnection>::const_iterator it_spall(spall.begin());
	map<pair<int,int>, SparsityInfo::NonlinearConnection>::const_iterator it_spfunc(spfunc.begin());
	VariableIterator it_nonlinvar(*sparsity[blocknr], false);
	Number* diagval=values+spall.size();
	while (it_spfunc!=spfunc.end()) {
		int col=it_spfunc->first.first;
		while (it_spall->first.first<col) { // skip columns in sp_all but not in spfunc
			int spallcol=it_spall->first.first;
			while (it_spall->first.first==spallcol) {
				if (!add) *values=0.;
				++values;
				++it_spall; //assert(it_spall!=spall.end());
			}
		}
		while (it_nonlinvar()<col) { // handle nonlinear columns not in spfunc (diagonal elements)
			if (func.get_sparsity(blocknr).nonlinear->count(it_nonlinvar())) {
				z[it_nonlinvar()]=1.;
				func.HessMult(y, xx, z, blocknr);
				if (add) *diagval+=factor*y[it_nonlinvar()]; else *diagval=factor*y[it_nonlinvar()]; // diagonal values
			} else {
				if (!add) *diagval=0.;
			}
			++diagval;
			++it_nonlinvar;
		}
		z[col]=1.;
		func.HessMult(y, xx, z, blocknr);
		if (add) *diagval+=factor*y[col]; else *diagval=factor*y[col]; // diagonal values
		while (it_spfunc!=spfunc.end() && (it_spfunc->first.first==col)) {
			while (it_spall->first.second!=it_spfunc->first.second) { // skip nonzeros not in this function
				++it_spall; //assert(it_spall!=spall.end()); 
				if (!add) *values=0.; ++values;
			}
			if (add) *values+=factor*y[it_spfunc->first.second]; else *values=factor*y[it_spfunc->first.second];
			++it_spall; ++it_spfunc; ++values;
		}
		z[col]=0.;
		++diagval;
		++it_nonlinvar;
	}
	while (it_nonlinvar) { // handle nonlinear columns not in spfunc (diagonal elements)
		if (func.get_sparsity(blocknr).nonlinear->count(it_nonlinvar())) {
			z[it_nonlinvar()]=1.;
			func.HessMult(y, xx, z, blocknr);
			if (add) *diagval+=factor*y[it_nonlinvar()]; else *diagval=factor*y[it_nonlinvar()]; // diagonal values
		} else {
			if (!add) *diagval=0.;
		}
		++diagval;
		++it_nonlinvar;
	}
	if (!add) {
		while (it_spall!=spall.end()) { *(values++)=0.; ++it_spall; }
	}
}

template <bool add>
void IpOptProblem::set_hessianquad(const SepQcFunc& func, int blocknr, Number* values, double factor) {
	const map<pair<int,int>, SparsityInfo::NonlinearConnection>& spall(*sparsity[blocknr]->sparsity_pattern);
	assert(func.sparsity_available(blocknr));
	const map<pair<int,int>, SparsityInfo::NonlinearConnection>& spfunc(*func.get_sparsity(blocknr).sparsity_pattern);
	Number* diagval=values+spall.size();

	map<pair<int,int>, SparsityInfo::NonlinearConnection>::const_iterator it_spall(spall.begin());
	map<pair<int,int>, SparsityInfo::NonlinearConnection>::const_iterator it_spfunc(spfunc.begin());
	while (it_spfunc!=spfunc.end()) {
		while (it_spall->first.second!=it_spfunc->first.second) { // skip nonzeros not in this function
			++it_spall; //assert(it_spall!=spall.end());
			if (!add) *values=0.; ++values;
		}
		// ipopt copies these elements to the upper right part; so set only to half values
		if (add) *values+=factor*it_spfunc->second.coeff; 
		else *values=factor*it_spfunc->second.coeff;
		++it_spfunc; ++it_spall; ++values;
	}
	if (!add) while (it_spall!=spall.end()) { *(values++)=0.; ++it_spall; }

	VariableIterator it_varall(*sparsity[blocknr], false);
	VariableIterator it_varfunc(func.get_sparsity(blocknr), false, false, true);
	while (it_varfunc) {
		while (it_varall()<it_varfunc()) {
			++it_varall; //assert(it_varall);
			if (!add) *diagval=0.; ++diagval;
		}
		if (add) *diagval+=2*factor*it_varfunc.coeff_quad();
		else *diagval=2*factor*it_varfunc.coeff_quad();
		++it_varall; ++it_varfunc; ++diagval;
	}
	if (!add) while (it_varall) { *(diagval++)=0.; ++it_varall; }
}


bool IpOptProblem::eval_h(Index n, const Number* x, bool new_x,
		   Number obj_factor, Index m, const Number* lambda,
		   bool new_lambda, Index nele_hess, Index* iRow,
		   Index* jCol, Number* values) {
	int count=0, oldcount;
  if (values == NULL) {
		for (int k=prob->block.size()-1; k>=0; --k) {
			if (sparsity[k] && sparsity[k]->sparsity_pattern) {
				for (map<pair<int,int>, SparsityInfo::NonlinearConnection>::iterator it(sparsity[k]->sparsity_pattern->begin()); it!=sparsity[k]->sparsity_pattern->end(); ++it, ++count) {
					jCol[count]=prob->block[k][it->first.first];
					iRow[count]=prob->block[k][it->first.second];
//					out_log << prob->var_names[prob->block[k][it->first.first]] << '\t' << prob->var_names[prob->block[k][it->first.second]] << endl;
				}
				for (VariableIterator it(*sparsity[k], false); it; ++it, ++count) {
					jCol[count]=prob->block[k][it()];
					iRow[count]=prob->block[k][it()];
//					out_log << prob->var_names[prob->block[k][it()]] << ' ';
				}
//				out_log << endl;
			}	else {
				for (int i=0; i<prob->block[k].size(); ++i)
					for (int j=0; j<prob->block[k].size(); ++j, ++count) {
						jCol[count]=prob->block[k][i];
						iRow[count]=prob->block[k][j];
					}
			}
		}
		assert(count==nele_hess);
//		out_log << "jCol: " << ivector(jCol, count) << "iRow: " << ivector(iRow, count);
	} else {
		dvector* xx=NULL; dvector* y; dvector* z;
		for (int k=prob->block.size()-1; k>=0; --k) {
			if (sparsity[k] && sparsity[k]->sparsity_pattern) {
				if (obj_factor && (prob->obj->A[k] || prob->obj->s[k])) {
					if (prob->obj->s[k]) {
						if (!xx) {
							xx=new dvector(prob->block[k].size());
							for (int i=xx->dim()-1; i>=0; --i) (*xx)[i]=x[prob->block[k][i]];
							z=new dvector(xx->dim());
							y=new dvector(xx->dim());
						}
						set_hessian<false>(*prob->obj, k, values+count, obj_factor, *xx, *y, *z);
					} else set_hessianquad<false>(*prob->obj, k, values+count, obj_factor);
				} else {
					for (int i=count+sparsity[k]->sparsity_pattern->size()+sparsity[k]->nr_var(false)-1; i>=count; --i) values[i]=0.;
				}

				for (int c=0; c<m; ++c) {
					if (!lambda[c]) continue;
					if (prob->con[c]->s[k]) {
						if (!xx) { 
							xx=new dvector(prob->block[k].size());
							for (int i=xx->dim()-1; i>=0; --i) (*xx)[i]=x[prob->block[k][i]];
							z=new dvector(xx->dim());
							y=new dvector(xx->dim());
						}
						set_hessian<true>(*prob->con[c], k, values+count, lambda[c], *xx, *y, *z);
					} else if (prob->con[c]->A[k]) set_hessianquad<true>(*prob->con[c], k, values+count, lambda[c]);
				}
				count+=sparsity[k]->sparsity_pattern->size()+sparsity[k]->nr_var(false); // this block is done

			} else {
				oldcount=count;
				if (obj_factor && (prob->obj->A[k] || prob->obj->s[k])) {
					if (!xx) { 
						xx=new dvector(prob->block[k].size());
						for (int i=xx->dim()-1; i>=0; --i) {
							(*xx)[i]=x[prob->block[k][i]];
						}
						z=new dvector(xx->dim());
						y=new dvector(xx->dim());
					}
					for (int i=0; i<xx->dim(); ++i) {
						(*z)[i]=1.;
						prob->obj->HessMult(*y, *xx, *z, k);
						for (int j=0; j<xx->dim(); ++j, ++count) {
							if (j!=i) values[count]=.5*obj_factor*(*y)(j); // looks light IPOPT do not like hessians with lower left and upper right part set :-(
							else values[count]=obj_factor*(*y)(j);
						}
						(*z)[i]=0.;
					}
				} else {
					for (int i=0; i<prob->block[k].size(); ++i)
						for (int j=0; j<prob->block[k].size(); ++j, ++count)
							values[count]=0.;
				}
				count=oldcount;

				for (int c=0; c<m; ++c) {
					oldcount=count;
					if (lambda[c] && (prob->con[c]->A[k] || prob->con[c]->s[k])) {
						if (!xx) { 
							xx=new dvector(prob->block[k].size());
							for (int i=xx->dim()-1; i>=0; --i) (*xx)[i]=x[prob->block[k][i]];
							z=new dvector(xx->dim());
							y=new dvector(xx->dim());
						}
						for (int i=0; i<xx->dim(); ++i) {
							(*z)[i]=1.;
							prob->con[c]->HessMult(*y, *xx, *z, k);
							for (int j=0; j<xx->dim(); ++j, ++count) {
								if (i==j) values[count]+=lambda[c]*(*y)(j);
								else values[count]+=.5*lambda[c]*(*y)(j);
							}
							(*z)[i]=0.;
						}
						count=oldcount;
					}
				}
				count+=prob->block[k].size()*prob->block[k].size();
			}
			if (xx) { delete xx; delete y; delete z; xx=NULL; }
		}
		assert(count==nele_hess);
	}

  return true;
}


void IpOptProblem::finalize_solution(SolverReturn status,
	Index n, const Number* x, const Number* z_L, const Number* z_U,
	Index m, const Number* g, const Number* lambda, Number obj_value)
{
	ipopt.sol_point.set(x, n);
	ipopt.lambda.set(lambda,m);

	switch (status) {
		case DIVERGING_ITERATES:
		case INVALID_NUMBER_DETECTED:
			ipopt.opt_val_=INFINITY;
			break;
		default:
			ipopt.opt_val_=prob->obj->eval(ipopt.sol_point);
	}

}

// ---------------------------------- IpOpt -------------------------------------

IpOpt::IpOpt(const Pointer<MinlpProblem> prob_, Pointer<Param> param_, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: LocOpt(prob_->dim(), out_solver_p_, out_solver_log_p_), param(param_), ipopt(new IpoptApplication()),
  lambda(prob_->con.size())
{	tol=1E-4;

	SmartPtr<OptionsList> options=ipopt->Options();

	ipopt->Options()->SetNumericValue("nlp_lower_bound_inf", -INFINITY);
	ipopt->Options()->SetNumericValue("nlp_upper_bound_inf", INFINITY);
	ipopt->Options()->SetIntegerValue("print_level", 2);

	ipoptproblem=new IpOptProblem(prob_, *this);
	sol_point=prob_->primal_point;
	reinit();
	
	ipopt->Initialize();
}

void IpOpt::reinit() {
}

int IpOpt::solve(dvector& start) {
	bool need_reinit=false;

	MinlpProblem& prob(*((IpOptProblem*)GetRawPtr(ipoptproblem))->prob);

	for (int i=0; (!need_reinit) && i<prob.i_discr.size(); ++i)
		need_reinit=(fabs(start(prob.i_discr[i])-sol_point(prob.i_discr[i]))>0);
  sol_point=start;
  if (need_reinit) reinit();
	return solve();
}

int IpOpt::solve() {
	ipopt->Options()->SetNumericValue("tol", tol);
	if (iter_max<INF) ipopt->Options()->SetIntegerValue("max_iter", iter_max);

	ApplicationReturnStatus status;
	try {
		timer.start();
		status = ipopt->OptimizeTNLP(ipoptproblem);
	} catch (IpoptException e) {
		status = Unrecoverable_Exception;
		out_err << "Caught IpoptException: " << e.Message() << endl;
	}

	time_=timer.stop();

	out_solver << "IPOPT returned " << status;
	if (status==Solve_Succeeded) {
		iter_=ipopt->Statistics()->IterationCount();
// 		opt_val_=ipopt->Statistics()->FinalObjective();
		out_solver << "\t iter: " << iter() << "\t optval " << opt_val() << endl;
	} else out_solver << endl;

//	out_log << "sol: " << sol_point;
	MinlpProblem& prob(*((IpOptProblem*)GetRawPtr(ipoptproblem))->prob);
	if (status==Solve_Succeeded) prob.feasible(sol_point, 1E-4, out_log_p);

	switch (status) {
		case Solve_Succeeded: return 0;
		case Infeasible_Problem_Detected: return 1;
		case Restoration_Failed:
		case Not_Enough_Degrees_Of_Freedom:
		case Insufficient_Memory:
		case Search_Direction_Becomes_Too_Small:
		case Invalid_Problem_Definition:
		case Invalid_Option:
		case Unrecoverable_Exception:
		case NonIpopt_Exception_Thrown:
		case Internal_Error: return 2;
		case Maximum_Iterations_Exceeded: return 3;
		case Solved_To_Acceptable_Level: return 4;
	}
	return 1;
}

#endif

