// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "rmp.h"
#include "osi.h"

RMPManager::RMPManager(Pointer<LinearRelax> linear_relax_, Pointer<MinlpNode> node_, Pointer<Param> param_)
: linear_relax(linear_relax_), node(node_), param(param_),
	coresize(linear_relax_->core_size()), blocknr(node_->i_ExtremePoints.size()), linrelaxdim(linear_relax_->obj->dim())
{ y_penalty_delta=param->get_d("RMP delta factor", 1E+6);
	init();
}

Pointer<SparseVector<double> > RMPManager::compute_col(const ExtremePoint& w, int k) {
	Pointer<SparseVector<double> > col=new SparseVector<double>(w.rmpcolumn); // from (R)
	col->resize(RMP->nr_row());
	col->SetElement(coresize+k, 1.); // convex combination constraint
	for (int j=0; j<w.dim(); ++j) // box constraint for y
		col->SetElement(coresize+blocknr+linear_relax->obj->block[k](j), w(j));

	return col;
}

void RMPManager::init() {
	// y-variables
	// constraint from core (R), convex-combination constraints, box on Wz+y
	MipProblem basermp(linrelaxdim, coresize+blocknr+linrelaxdim);

	// add core of (R)
	int c=0;
	for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->couple_con.begin()); it!=linear_relax->couple_con.end(); ++it, ++c)
		for (int k=0; k<it->b.size(); ++k)
			if (it->b[k]) {
				Pointer<UserVector<double> > row(it->b[k]->getcopy()); *row*=2.;
				basermp.setRow(c, *row, linear_relax->obj->block[k], it->eq ? -it->c : -INFINITY, -it->c);
			}
	for (int k=0; k<blocknr; ++k)
		for (list<LinearRelax::LinConstraint>::iterator it(linear_relax->block_con[k].begin()); it!=linear_relax->block_con[k].end(); ++it, ++c) {
			Pointer<UserVector<double> > row(it->b[0]->getcopy()); *row*=2.;
			basermp.setRow(c, *row, linear_relax->obj->block[k], it->eq ? -it->c : -INFINITY, -it->c);
		}

	assert(c==coresize);
	// add constant parts of convex-combination constraints
	for (int k=0; k<blocknr; ++k, ++c) basermp.setRowBounds(c, 1., 1.);

	// add y - boxconstraints ...
	for (int i=0; i<linrelaxdim; ++i, ++c) {
		SparseVector<double> row(basermp.dim(), i, 1.);
		basermp.setRow(c, row, node->lower(i), node->upper(i)); // \lb x(i) <= y_i <= \ub x(i)
	}

	// empty objective
	basermp.setObj(new SparseVector<double>(linrelaxdim), linear_relax->obj->c);

	basermp.finish();

	RMP=MIPSolver::get_solver(basermp, param);

	// compute columns
	vector<Pointer<UserVector<double> > > columns;
	vector<double> objcoeff;
	for (int k=0; k<blocknr; ++k)
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it) {
			columns.push_back((Pointer<UserVector<double> >)compute_col(**it, k));
			objcoeff.push_back((**it).rmpobjcoeff);
		}
	dvector low(columns.size()), up(columns.size(), 1.); // lower and upper bounds for z

	// add columns
	list<const MIPSolver::ColItem*> colitems;
	RMP->add_cols(colitems, columns, objcoeff, low, up);

	// store ColItems
	list<const MIPSolver::ColItem*>::iterator it_col(colitems.begin());
	for (int k=0; k<blocknr; ++k)
		for (list<list<ExtremePoint>::iterator>::iterator it_pt(node->i_ExtremePoints[k].begin()); it_pt!=node->i_ExtremePoints[k].end(); ++it_pt, ++it_col)
			(**it_pt).rmpcolitem=*it_col;
}

int RMPManager::colindex(const ExtremePoint& w) {
	assert(w.rmpcolitem);
	return RMP->get_colindex(*w.rmpcolitem);
}

void RMPManager::add_column(ExtremePoint& w, int k) {
	w.rmpcolitem=RMP->add_col(*compute_col(w, k), w.rmpobjcoeff, 0., 1., MipProblem::CONTINUOUS);
}

void RMPManager::update_column(const ExtremePoint& w, int k) {
	RMP->modify_col(*w.rmpcolitem, *compute_col(w, k), w.rmpobjcoeff, 0., 1., MipProblem::CONTINUOUS);
}

void RMPManager::update_s(const dvector& x, bool set_start) {
	ivector s(x.dim(), 1);
	for (int k=0; k<blocknr; ++k) {
		dvector mid(linear_relax->obj->block[k].size());
		int size=0;
		for (list<list<ExtremePoint>::iterator >::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it, ++size)
			mid+=**it;
		mid/=size;
		for (int i=0; i<mid.dim(); i++) {
			int i0=linear_relax->obj->block[k][i];
			if (x(i0)-mid(i)<-rtol) s[i0]=-1;
			else if (x(i0)-mid(i)<rtol) s[i0]=0;
			if (set_start) node->yz_RMP[i0]=x(i0)-mid(i);
		}
		if (set_start)
			for (list<list<ExtremePoint>::iterator >::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it)
				node->yz_RMP[colindex(**it)]=1./(double)size;
	}
	// penalty: delta * MAX(||c||,1)
	y_penalty_scale=1.;
//	for (int k=0; k<node->i_RMP_points.size(); k++)
//		if (RMP->obj->b[k]) y_penalty_scale+=RMP->obj->b[k]->sq_norm2();
//	if (y_penalty_scale>1) y_penalty_scale=sqrt(y_penalty_scale);
//	else y_penalty_scale=1;
	for (int i=0; i<linrelaxdim; ++i) {
		RMP->modify_obj(i, s(i)*y_penalty_delta*y_penalty_scale);
		RMP->modify_col(i, s(i)<0 ? -INFINITY : 0., s(i)>0 ? INFINITY : 0, MipProblem::CONTINUOUS);
	}
}

void RMPManager::restrict_y(const dvector& y) {
	for (int i=0; i<y.dim(); ++i) {
		if (y(i)>0) RMP->modify_col(i, 0, y(i), MipProblem::CONTINUOUS);
		else RMP->modify_col(i, y(i), 0, MipProblem::CONTINUOUS);
		RMP->modify_obj(i, 0.);
	}
}

void RMPManager::get_x(UserVector<double>& x, const dvector &yz) {
	x=0;
	for (int k=0; k<node->i_ExtremePoints.size(); ++k)
		for (list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin()); it!=node->i_ExtremePoints[k].end(); ++it) {
			int RMP_index=colindex(**it);
			if (!yz(RMP_index)) continue;
			for (int j=(**it).dim()-1; j>=0; --j)
				x[linear_relax->obj->block[k][j]]+=yz(RMP_index)*(**it)(j);
		}
}

int RMPManager::prune() {
	list<const MIPSolver::ColItem*> colitems;
	int count=0;
	for (int k=0; k<node->i_ExtremePoints.size(); ++k) {
//		out_log << k << '\t' << node->i_ExtremePoints[k].size() << ": ";
		double sum=0;
		list<list<ExtremePoint>::iterator>::iterator it(node->i_ExtremePoints[k].begin());
		while (it!=node->i_ExtremePoints[k].end()) {
//			out_log << RMP->get_primal(*(**it).rmpcolitem) << ' ';
			sum+=RMP->get_primal(*(**it).rmpcolitem);
			if (fabs(RMP->get_primal(*(**it).rmpcolitem))<rtol) {
				colitems.push_back((**it).rmpcolitem);
				it=node->i_ExtremePoints[k].erase(it);
				++count;
			} else ++it;
		}
//		out_log << "\t " << " -> " << node->i_ExtremePoints[k].size() << '\t' << sum << endl;
		assert(!node->i_ExtremePoints[k].empty());
		assert(fabs(sum-1)<1E-4);
	}
	if (count) RMP->delete_cols(colitems);
	return count;
}
