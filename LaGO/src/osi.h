// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef OSI_H
#define OSI_H

#include "standard.h"
#include "problem.h"
#include "opt.h"

class OsiSolverInterface;
class CglCutGenerator;

/** Interface to COIN OSI.
*/
class OSISolver : public MIPSolver {
	public:
		class RowItem : public MIPSolver::RowItem {
			friend class OSISolver;
			private:
				int index;
				list<RowItem>::iterator it;
			public:
				RowItem(int index_)
				: index(index_)
				{ }
		};

		class ColItem : public MIPSolver::ColItem {
			friend class OSISolver;
			private:
				int index;
				list<ColItem>::iterator it;
			public:
				ColItem(int index_)
				: index(index_)
				{ }
		};

	private:
		/** The OSI-Solver, we use.
		*/
		Pointer<OsiSolverInterface> osisolver;
		
		Pointer<CglCutGenerator> cutgenerator;

		/** To constant of the objective function.
		*/
		double obj_const;

		/** Indicates, whether the next solve() will be a cold start or not.
		*/
		bool cold_start;

		bool lastpoint_feasible();

		list<RowItem> addedrows;
		list<ColItem> addedcols;
	public:

		OSISolver(const MipProblem& mip);

		int nr_col();
		int nr_row();

		void set_tol(double tol);
		void set_maxiter(int maxiter);

		void reset();

		SolutionStatus solve();
		SolutionStatus solve(const UserVector<double>& x);
		SolutionStatus solveMIP();
		
		SolutionStatus feasible();

		void get_primal(UserVector<double>& x);
		double get_primal(const MIPSolver::ColItem& colitem);
		using MIPSolver::get_primal;
		
		int get_colindex(const MIPSolver::ColItem& colitem) { return ((ColItem*)(&colitem))->index; }

		void get_dual(UserVector<double>& mu);
		double get_dual(const MIPSolver::RowItem& rowitem);

		void get_reducedcosts(UserVector<double>& rc);
		double get_reducedcosts(const MIPSolver::ColItem& colitem);

		void get_rowactivity(UserVector<double>& rowact);
		double get_rowactivity(const MIPSolver::RowItem& rowitem);

		double get_optval();

		int get_iter();

		void set_obj(const UserVector<double>& obj, double obj_const_=0.);
		void modify_obj(int i, double coeff);

		const MIPSolver::RowItem* add_row(const UserVector<double>& row, double low, double up);
		const MIPSolver::RowItem* add_row(const UserVector<double>& row, const ivector& indices, double low, double up);
		void add_rows(const vector<pair<dvector, ivector> >& rows, const dvector& low, const dvector& up);
		void add_rows(list<const MIPSolver::RowItem*>& rowitems, const vector<pair<dvector, ivector> >& rows, const dvector& low, const dvector& up);
		void delete_row(const MIPSolver::RowItem& rowitem) { delete_rows(list<const MIPSolver::RowItem*>(1, (RowItem*)&rowitem)); }
		void delete_rows(const list<const MIPSolver::RowItem*>& rowitems);
		void modify_row(const MIPSolver::RowItem& rowitem, double low, double up);
		void modify_row(int index, double low, double up);

		const MIPSolver::ColItem* add_col(double low, double up, MipProblem::VarType vartype=MipProblem::CONTINUOUS);
		const MIPSolver::ColItem* add_col(const UserVector<double>& col, double obj_coeff, double low, double up, MipProblem::VarType vartype=MipProblem::CONTINUOUS);
		void add_cols(list<const MIPSolver::ColItem*>& colitems, const dvector& low, const dvector& up);
		void add_cols(list<const MIPSolver::ColItem*>& colitems, vector<Pointer<UserVector<double> > >& cols, const vector<double>& obj_coeff, const dvector& low, const dvector& up);
		void delete_col(const MIPSolver::ColItem& colitem) { delete_cols(list<const MIPSolver::ColItem*>(1, (ColItem*)&colitem)); }
		void delete_cols(const list<const MIPSolver::ColItem*>& colitems);
		void modify_col(const MIPSolver::ColItem& colitem, double low, double up, MipProblem::VarType type);
		void modify_col(const MIPSolver::ColItem& colitem, const UserVector<double>& col, double obj_coeff, double low, double up, MipProblem::VarType vartype);
		void modify_col(int index, double low, double up, MipProblem::VarType type);

		double get_collow(int index);
		double get_colup(int index);

		int generate_cuts(list<Pointer<SimpleCut> >& rowcuts);

//		void print();
};

#endif
