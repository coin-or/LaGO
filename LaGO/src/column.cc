// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "column.h"
#include "node.h"
#include "bcp.h"

// --------------------------- column generator -------------------------------------

ColumnGenerator::ColumnGenerator(Pointer<MinlpNode> node_, MinlpBCP& bcp_)
: node(node_), bcp(bcp_), param(bcp_.param), rc_tol(bcp_.gap_tol)
{ max_major_iter=param->get_i("ColumnGenerator max major iter", 5);
	max_minor_iter=param->get_i("ColumnGenerator max minor iter", 0);
	if (!max_minor_iter) {
		double bnr=node->i_ExtremePoints.size();
		double bsize=(double)bcp.linear_relax->obj->dim()/bnr;
	 	max_minor_iter=MAX(5, (int)(1.39*bsize-.05*bsize*bsize+0.15*bnr));
//		out_log << "formular gives: " << (1.39*bsize-.05*bsize*bsize+0.15*bnr) << endl;
	}
	// iteration limit in init_RMP:
	max_initRMP_iter=param->get_i("ColumnGenerator max init RMP iter", 0);
	if (!max_initRMP_iter) {
		for (int k=0; k<bcp.linear_relax->obj->block.size(); k++) max_initRMP_iter=MAX(max_initRMP_iter, bcp.linear_relax->obj->block[k].size());
		max_initRMP_iter*=2; max_initRMP_iter+=2;
//		max_initRMP_iter+=MAX(2, max_initRMP_iter/4);
	}
}

void ColumnGenerator::update_ExtremePoints() {
	// update extreme points
	for (int k=0; k<bcp.ExtremePoints.size(); ++k) {
		if (bcp.ExtremePoints[k].empty()) continue;
		if (node->i_ExtremePoints_limit[k]==bcp.ExtremePoints[k].end()) { // first point
			node->i_ExtremePoints_limit[k]=bcp.ExtremePoints[k].begin();
			if (node->inside_part_set(*node->i_ExtremePoints_limit[k], k, bcp.linear_relax->obj->block)) {
				node->i_ExtremePoints[k].push_back(node->i_ExtremePoints_limit[k]);
				if (RMP) RMP->add_column(*node->i_ExtremePoints_limit[k], k);
				++bcp.update_ExtremePoints_count;
			}
		}
		list<ExtremePoint>::iterator last(node->i_ExtremePoints_limit[k]++);
		while (node->i_ExtremePoints_limit[k]!=bcp.ExtremePoints[k].end()) {
			if (node->inside_part_set(*node->i_ExtremePoints_limit[k], k, bcp.linear_relax->obj->block)) {
				node->i_ExtremePoints[k].push_back(node->i_ExtremePoints_limit[k]);
				if (RMP) RMP->add_column(*node->i_ExtremePoints_limit[k], k);
				++bcp.update_ExtremePoints_count;
			}
			last=node->i_ExtremePoints_limit[k]++;
		}
		node->i_ExtremePoints_limit[k]=last;
	}
	if (RMP) node->yz_RMP.resize(RMP->dim());
}

void ColumnGenerator::mult_W(UserVector<double>& x, const dvector &z, int k) {
	x=0;
	int i=0;
	for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end() && i<z.dim(); ++it, ++i)
		x.AddMult(z[i], **it);
}

void ColumnGenerator::add_cut(Pointer<SepQcFunc> cut, int k) {
	double scale=sqrt(cut->b[0]->sq_norm2());
	if (bcp.linear_relax->obj->b[k]) scale/=sqrt(bcp.linear_relax->obj->b[k]->sq_norm2());
	if (scale>rtol) scale=1/scale; else scale=1.;
	if (scale<1E-3 || scale>1E+3) out_log << "Scaling cut by " << scale << endl;
	*cut->b[0]*=scale; cut->c*=scale;

	if (bcp.block_sub_convex_prob.size()) bcp.block_sub_convex_prob[k]->add_con(cut, false, "ColumnGenerator cut");
	bcp.linear_relax->add_cut(new SimpleCut(cut->b[0], cut->c), k, node);
}

void ColumnGenerator::solve_lag_problem(MinlpBCP::LagSolveStatus& status, int k, Pointer<SepQcFunc> temp_cut) {
	bcp.solve_lag_problem(status, k, node, temp_cut);
	if (status.ret || (!status.solset.size())) {
		out_log << "Solved lag_prob " << k << " ret: " << status.ret << "\t points: " << status.solset.size() << "\t lowbound: " << status.lowbound << "\t iterations: " << status.iter << endl;
		if (!status.solset.size()) status.ret=1;
		return;
	}
//	val=w.begin()->first;
//	sol_w=w.begin()->second;

	pair<list<ExtremePoint>::iterator, bool> ret(bcp.add_ExtremePoint(status.solset.begin()->second, k, node));

	// add Lagrangian cut
	if (bcp.lag_cuts && !(*bcp.lag_problem[k]->obj->b[0]==0.) && (!temp_cut)) {
		Pointer<SepQcFunc> cut(new SepQcFunc(*bcp.lag_problem[k]->obj, true)); // -lag func
		cut->c+=status.lowbound;
		add_cut(cut, k);
	}

	// add first RMP point
	if (ret.second) {
  	status.new_points++;
		if (RMP) RMP->add_column(*ret.first, k);
	}

	for (set<SolCandidate>::iterator it_w(++status.solset.begin()); it_w!=status.solset.end(); it_w++) {
		ret=bcp.add_ExtremePoint(it_w->second, k, node);
		if (!ret.second) continue; // point was not new (for this node)

  	status.new_points++;
		if (RMP) RMP->add_column(*ret.first, k);
	}

	return;
}

int ColumnGenerator::get_search_dir2(dvector& a, double& c, int k) {
	Pointer<MinlpProblem> prob(new MinlpProblem());

	dvector xk(node->ref_point, bcp.linear_relax->obj->block[k]);

	int size=node->i_ExtremePoints[k].size();
	int i;
	for (i=0; i<size; i++)
		prob->add_var(i, 0, false, 0, 1);

	double maxel=0;
	prob->add_obj(new SepQcFunc(prob->block));
	Pointer<DenseMatrix> WW(new DenseMatrix(prob->dim()));
	prob->obj->A[0]=WW;
	prob->obj->b[0]=new dvector(prob->dim());
	i=0;
	for (list<list<ExtremePoint>::iterator>::iterator it_i(node->i_ExtremePoints[k].begin()); it_i!=node->i_ExtremePoints[k].end(); i++, it_i++) {
		(*prob->obj->b[0])[i]=-xk* **it_i;
		int j=0;
		for (list<list<ExtremePoint>::iterator>::iterator it_j(node->i_ExtremePoints[k].begin()); j<i; j++, it_j++) {
			(*WW)(j,i)=(*WW)(i,j)=**it_i * **it_j;
			maxel=MAX(maxel, fabs((*WW)(i,j)));
		}
		(*WW)(i,i)=(**it_i).sq_norm2();
		maxel=MAX(maxel, fabs((*WW)(i,i)));
	}
	prob->obj->c=xk.sq_norm2();

	maxel/=10.;
	if (maxel>1) { // scale by max element of matrix
		for (int i=0; i<WW->dim(); i++)
			for (int j=0; j<WW->dim(); j++)
				(*WW)(i,j)/=maxel;
		*prob->obj->b[0]/=maxel;
		prob->obj->c/=maxel;
	}

	prob->add_con(new SepQcFunc(prob->block));
	prob->con[0]->b[0]=new dvector(prob->dim(), .5);
	prob->con[0]->c=-1;

	Pointer<LocOpt> locopt(LocOpt::get_solver(prob, param, "get search dir", NULL, NULL));
	int ret=locopt->solve();
	if (ret && prob->feasible(locopt->sol_point, 1E-4, NULL)) {
		out_log << "get_search_dir2 failed, infeasible solution, return " << ret << ", iterations: " << locopt->iter() << "\t dim: " << prob->dim() << endl;
		return 1;
	}
	if (ret) out_log << "get_search_dir2 return: " << ret << ", feasible solution, iterations: " << locopt->iter() << "\t dim: " << prob->dim() << endl;
	double val=locopt->opt_val(); if (maxel>10) val*=maxel;
//	out_log << "get_search_dir2 return: " << ret << "\t value: " << val << endl;
	mult_W(a, locopt->sol_point, k);
	if (val<1E-6) {
		int i=0;
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it, ++i)
			node->yz_RMP[RMP->colindex(**it)]=locopt->sol_point(i);
		for (i=0; i<a.dim(); i++)
			node->yz_RMP[bcp.linear_relax->obj->block[k][i]]=a(i)-xk(i);
		return 2;
	}

	double t=.75;
	c=t*a.sq_norm2()+(1-2*t)*(a*xk)-(1-t)*prob->obj->c;
	a-=xk;
//	out_log << "a xk: " << a*xk << endl << "c: " << c << endl;
	return 0;
}

int ColumnGenerator::init_RMP() {
	Pointer<MinlpProblem> prob(bcp.reform ? bcp.reform->ext_prob : bcp.split_prob);
	dvector y(prob->dim());
	int iter=0;
	bool stop;
	MIPSolver::SolutionStatus ret;

	do {
		node->yz_RMP.resize(RMP->dim());
		stop=true;
		out_log << "init_RMP iteration " << iter << " / " << max_initRMP_iter << endl;
		RMP->update_s(node->ref_point, true);
		//Compute dual and primal solutions  $\mu^j$ and $z^j$ of RMP
		ret=RMP->solve(node->yz_RMP);

		if (ret!=MIPSolver::SOLVED) { // update s by computing a new refpoint
			out_log << "Solving RMP return: " << ret << endl;

			if (bcp.set_LP_bound(node)) return 3;
			if (node->low_bound>bcp.opt_val()) return 2;

			RMP->update_s(node->ref_point, true);
			ret=RMP->solve(node->yz_RMP);
			out_log << "Solving RMP again return: " << ret << "\t RMP dim: " << RMP->dim() << endl;

			if (ret!=MIPSolver::SOLVED) return 1;
		}

		RMP->get_yz(node->yz_RMP);
		RMP->get_x(node->ref_point, node->yz_RMP);

		bool y_changed=false;
		double y_2norm=0.;
		list<pair<int, pair<dvector, double> > > nonzero_blocks;
		for (int k=0; k<prob->block.size(); k++) {
			double norm=0;
			for (int i=0; i<prob->block[k].size(); i++) {
				int i0=prob->block[k][i];
				y_changed|=fabs(node->yz_RMP(i0)-y(i0))>1E-3*(prob->upper(i0)-prob->lower(i0))+rtol;
				y[i0]=node->yz_RMP(i0);
				norm+=y(i0)*y(i0);
				node->ref_point[i0]+=y(i0);
			}
			if (sqrt(norm)>1E-3) {
				pair<int, pair<dvector, double> > el(k, pair<dvector, double>(dvector(prob->block[k].size()), 0.));
				int ret=get_search_dir2(el.second.first, el.second.second, el.first);
				if (ret==2) { // could find a better convex comb. for \hat x_J_k.
					norm=0;
					for (int i=0; i<prob->block[k].size(); ++i) {
						int i0=prob->block[k][i];
						y[i0]=node->yz_RMP(i0);
						norm+=y(i0)*y(i0);
					}
				} else if (ret==1) {
					out_log << "Couldn't compute new search direction, return: " << ret << endl;
				} else
					nonzero_blocks.push_back(el);
			}
			y_2norm+=norm;
		}
		if (y_changed) stop=false;
		y_2norm=sqrt(y_2norm);

		RMP->get_x(node->ref_point, node->yz_RMP);
		node->ref_point+=y;

		if (nonzero_blocks.empty()) { // we are happy
			if (1 || y_2norm>rtol) out_log << "Solving RMP: y is almost zero. value: " << RMP->get_optval() << "\t ||y||: " << y_2norm << endl;
			RMP->get_dual(node->dual_point);

			RMP->restrict_y(y);
			ret=RMP->solve();

			if (ret==MIPSolver::SOLVED) {
				RMP->get_yz(node->yz_RMP);
				RMP->get_x(node->ref_point, node->yz_RMP);
				RMP->get_dual(node->dual_point);
				int pruned=RMP->prune();
//				out_log << "Pruned " << pruned << " Extreme points from node." << endl;
				node->yz_RMP.resize(RMP->dim());
				RMP->get_yz(node->yz_RMP);
			} else
				out_log << "Solving RMP with restricted y: " << ret << "\t value: " << RMP->get_optval() << endl;

			return 0;
		}

		if (++iter>=max_initRMP_iter) break;

		// trying to improve the RMP
		out_log << "Solving RMP: y is not zero. value: " << RMP->get_optval() << "\t ||y||: " << y_2norm << "\t problematic blocks: " << nonzero_blocks.size() << endl;
		int new_points=0;
		bool recompute_LP_bound=false;
		for (list<pair<int, pair<dvector, double> > >::iterator it(nonzero_blocks.begin()); it!=nonzero_blocks.end(); it++) {
			int k=it->first;
//			out_log << "y_k norm: " << sqrt(it->second.first.sq_norm2()) << endl;

			Pointer<SepQcFunc> cut(new SepQcFunc(bcp.lag_problem[k]->block));
			cut->c=-it->second.second;
			cut->b[0]=new SparseVector<double>(it->second.first); *cut->b[0]*=.5;

			bcp.update_lag_problem(k, it->second.first);
			MinlpBCP::LagSolveStatus status; status.new_points=0;
			solve_lag_problem(status, k, cut);
			new_points+=status.new_points;
//			out_log << "lag prob " << k << ":\t return: " << status.ret << "\t new points: " << status.new_points << endl;
			if (status.ret && status.lowbound==-INFINITY) { // add cut
				*cut->b[0]*=-1; cut->c*=-1.;
				add_cut(cut, k);
//				out_log << "hat x_k in new cut: " << cut->eval(node->ref_point(bcp.linear_relax->prob->block[k])) << endl;
				recompute_LP_bound=true;
			} else if(status.lowbound>it->second.first*node->ref_point(prob->block[k])) {
				*cut->b[0]*=-1; cut->c=status.lowbound;
				add_cut(cut, k);
				recompute_LP_bound=true;
			}
		}
		if (recompute_LP_bound) {
			if (bcp.set_LP_bound(node)) return 3;
//			out_log << node->ref_point;
			stop=false;
		} else if (!new_points) {
			out_log << "No new RMP-points, but all lag problems had lower bounds!" << endl;
		}
		if (new_points) {
//			out_log << "Added " << new_points << " new points." << endl;
			update_ExtremePoints();
			stop=false;
		}
	} while (!stop);
	if (iter==max_initRMP_iter) out_log << "init_RMP: Iteration limit exceeded" << endl;

	RMP->restrict_y(y);

	ret=RMP->solve();
	out_log << "Solving RMP with restricted y: " << ret << "\t value: " << RMP->get_optval() << endl;
	if (ret==MIPSolver::SOLVED) {
		RMP->get_yz(node->yz_RMP);
		RMP->get_x(node->ref_point, node->yz_RMP);
		RMP->get_dual(node->dual_point);
		int pruned=RMP->prune();
//		out_log << "Pruned " << pruned << " Extreme points from node." << endl;
		node->yz_RMP.resize(RMP->dim());
		RMP->get_yz(node->yz_RMP);
		return 1;
	}

	return 1;
}

int ColumnGenerator::generate_RMP() {
	if (!bcp.lag_problem.size()) bcp.init_lag_problems(node);
	if (!RMP) {
		update_ExtremePoints();
		RMP=new RMPManager(bcp.linear_relax, node, param);
	}
	int iter=0, minor_iter;
	int ret;
	double old_RMP_val;
//	bool bound_computation_succeeded=false;
	bool RMP_solve_succeeded=false;
	double lower;
	double max_rc,rc,val_inner,val_outer;
	dvector x_RMP(node->ref_point.dim());

	if (bcp.set_LP_bound(node)) return 1;
	if (node->low_bound>bcp.opt_val()-rtol) return 0;

	do {
		Timer t;
		int RMP_ret=init_RMP();
		bcp.init_RMP_time+=t.stop();
		if (RMP_ret==2) return 0; // (R[U]) is feasible, but threshold reashed
		if (RMP_ret==3) return 1; // (R[U]) is infeasible
		if (RMP_ret==1) { // (R[U]) feasible, but couldn't make (RMP) feasible
			out_log << "Couldn't make (RMP) feasible. Proceeding with solution of (R[U])." << endl;
//			return 0;
		}
		x_RMP=node->ref_point;
		if (RMP_ret==0) {
			RMP_solve_succeeded=true;

			double new_RMP_val=RMP->get_optval();
			if (fabs(old_RMP_val-new_RMP_val)/fabs(new_RMP_val+1)<.001) {
				out_log << "No improvement of the RMP." << endl;
				break;
			}
			old_RMP_val=new_RMP_val;
		}
		old_RMP_val=RMP->get_optval();

		bcp.update_lag_problems(node->dual_point);
		max_rc=0;
// compute LP-reduced costs and sort lagproblems into subprobl_list

		// Indicates which Lagrangian sub-problems are active
		multimap<double, int> subprobl_list;
		for(int k=0; k<node->i_ExtremePoints.size(); ++k) {
			dvector x_k(node->ref_point, bcp.linear_relax->obj->block[k]);
			val_inner=bcp.lag_problem[k]->obj->eval(x_k);

			LinearRelax block_sub_linear_relax(*bcp.linear_relax, k, node, bcp.lag_problem[k]->obj);
			dvector dummy(x_k.dim());
			ret=block_sub_linear_relax.solve(dummy, val_outer, NULL);
			if (ret) {
				out_log << "LP subproblem " << k << " not feasible." << endl;
				continue;
			}

			rc=2*(val_inner-val_outer)/(fabs(val_inner)+fabs(val_outer)+1);
			max_rc=MAX(max_rc, fabs(rc));
//			out_log << k << "\t outer: " << val_outer << "\t inner: " << val_inner << "\t rc: " << rc << endl;
			if (rc>rc_tol) subprobl_list.insert(pair<double, int>(subprobl_key(rc,k),k));
		}
		minor_iter = max_minor_iter;//*bcp.linear_relax->prob->con.size() * bcp.linear_relax->prob->dim();
		out_log << "max approx rc: " << max_rc << endl;
		max_rc=0;
		while (minor_iter && subprobl_list.size()) {
			int k=subprobl_list.begin()->second;
			double old_rc=-subprobl_list.begin()->first;
			subprobl_list.erase(subprobl_list.begin());
			//Solve the $k$-th Lagrangian subproblem
			MinlpBCP::LagSolveStatus status; status.new_points=0;
			solve_lag_problem(status, k);
			if (status.ret) continue;
			double val_inner=bcp.lag_problem[k]->obj->eval(node->ref_point(bcp.linear_relax->obj->block[k]));
//			max_rc=MAX(max_rc, fabs(rc));
			rc=2*(val_inner-status.lowbound)/(fabs(val_inner)+fabs(status.lowbound)+1);
			max_rc=MAX(max_rc, rc);
			out_log << k << "\t outer: " << status.lowbound << "\t inner: " << val_inner << "\t old rc: " << old_rc << "\t rc: " << rc
				<< "\t iter: " << status.iter << "\t gap: " << (status.value-status.lowbound)/(1+fabs(status.lowbound));
			if (old_rc-rc>.1) { minor_iter--; out_log << "\t " << minor_iter; }
			out_log << endl;
		}
		if (subprobl_list.size()) max_rc=MAX(max_rc, -subprobl_list.begin()->first);
		ret=bcp.set_LP_bound(node);
		if (ret) return ret;
		if (node->low_bound>bcp.opt_val()) {
			if (!RMP_solve_succeeded) bcp.bound_failed++;
			return 0;
		}
	} while(iter++<max_major_iter && max_rc>rc_tol); // or no significant improvement of  (R)
	node->yz_RMP.resize(RMP->dim());
	if (iter>max_major_iter) out_log << "Iteration limit exceeded." << endl;
	if (!RMP_solve_succeeded) bcp.bound_failed++;
	node->ref_point=x_RMP;
	return 0;
}

double  ColumnGenerator::subprobl_key(double rc, int k) {
// {\rm gap\_impr}[k]=|{\rm rc\_LP}[k]|/(n_k^2\cdot m_k)
	assert(bcp.lag_problem[k]!=NULL);
	return -rc;
//	return -rc/((double) bcp.lag_problem[k]->dim()*bcp.lag_problem[k]->dim()*bcp.lag_problem[k]->con.size());
}

int ColumnGenerator::solve_RMP() {
	RMP=NULL;
	update_ExtremePoints();
	RMP=new RMPManager(bcp.linear_relax, node, param);

	// solve RMP
	Timer t;
	int ret=init_RMP();
	bcp.init_RMP_time+=t.stop();
	return ret;
}

