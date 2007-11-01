// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "osi.h"

#define COIN_BIG_INDEX 0

#ifdef COIN_HAS_CPX
#define CPLEX_AVAILABLE
#endif

#include "OsiSolverInterface.hpp"
#ifdef CPLEX_AVAILABLE
#include "OsiCpxSolverInterface.hpp"
#else
//#include "coin/OsiVolSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#endif
#include "CoinPackedVector.hpp"

OSISolver::OSISolver(const MipProblem& mip)
: MIPSolver(), obj_const(mip.getObjConst()), cold_start(true)
{
#ifdef CPLEX_AVAILABLE
	OsiCpxSolverInterface* osisolvercpx=new OsiCpxSolverInterface();
	// switch off CPLEX output
	CPXsetintparam(osisolvercpx->getEnvironmentPtr(), CPX_PARAM_SCRIND, CPX_OFF );
	CPXsetlogfile(osisolvercpx->getEnvironmentPtr(), NULL); // maybe not necessary
	CPXsetdblparam(osisolvercpx->getEnvironmentPtr(), CPX_PARAM_TILIM, 60.); // 60 second time limit
	osisolver=osisolvercpx;
#else
	OsiClpSolverInterface* osisolverclp=new OsiClpSolverInterface();
//	osisolverclp->getModelPtr()->setDblParam(ClpMaxSeconds, 60.);
	osisolverclp->getModelPtr()->setIntParam(ClpMaxNumIteration, 10000);
//	osisolverclp->getModelPtr()->messageHandler()->setLogLevel(0);
//	osisolverclp->setHintParam(OsiDoScale, false);
	osisolver=osisolverclp;
//	osisolver=new OsiVolSolverInterface();
#endif
	osisolver->messageHandler()->setLogLevel(0);
	osisolver->loadProblem(mip.dim(), mip.rows(),
		mip.getMatrix().GetColPtr(), mip.getMatrix().GetRowInd(), mip.getMatrix().GetVal(),
		(const Pointer<double>)mip.getColLower(), (const Pointer<double>)mip.getColUpper(),
		(const Pointer<double>)mip.getObj(),
		(const Pointer<double>)mip.getRowLower(), (const Pointer<double>)mip.getRowUpper()
	);
	osisolver->setDblParam(OsiObjOffset, -obj_const);

	for (int i=0; i<nr_col(); i++) {
		switch (mip.getColType(i)) {
			case MipProblem::CONTINUOUS: break;
			case MipProblem::BINARY:
			case MipProblem::INTEGER:
				osisolver->setInteger(i);
//				out_log << "found " << i << " to be integer " << mip.getColLower()(i) << '\t' << mip.getColUpper()(i)
//				<< '\t' << osisolver->getColLower()[i] << '\t' << osisolver->getColUpper()[i] << endl;
				break;
			default: out_err << "OSISolver: Variable type unknown. Aborting." << endl; exit(-1);
		}
		if (mip.getColLower()(i)==-INFINITY) osisolver->setColLower(i, -osisolver->getInfinity());
		if (mip.getColUpper()(i)==INFINITY) osisolver->setColUpper(i, osisolver->getInfinity());
	}

	for (int i=0; i<nr_row(); i++) {
		if (mip.getRowLower()(i)==-INFINITY) osisolver->setRowLower(i, -osisolver->getInfinity());
		if (mip.getRowUpper()(i)==INFINITY) osisolver->setRowUpper(i, osisolver->getInfinity());
	}
};

int OSISolver::nr_col() { return osisolver->getNumCols(); }

int OSISolver::nr_row() { return osisolver->getNumRows(); }

void OSISolver::set_tol(double tol) { osisolver->setDblParam(OsiPrimalTolerance, tol); }

void OSISolver::set_maxiter(int maxiter) { osisolver->setIntParam(OsiMaxNumIteration, maxiter); }

void OSISolver::reset() {
	cold_start=true;
	// remove cols
	int Colnum=addedcols.size(); // slow :(
	int* Ind=NULL;
	if (Colnum) {
		Ind=new int[Colnum];
		int i=0;
		for (list<ColItem>::iterator it(addedcols.begin()); it!=addedcols.end(); ++it, ++i) Ind[i]=it->index;
		osisolver->deleteCols(Colnum, Ind);
		addedcols.clear();
	}

	// remove rows
	int Rownum=addedrows.size(); // slow :(
	if (Rownum) {
		if (Rownum>Colnum) { delete Ind; Ind=new int[Rownum]; }
		int i=0;
		for (list<RowItem>::iterator it(addedrows.begin()); it!=addedrows.end(); ++it, ++i) Ind[i]=it->index;
		osisolver->deleteRows(Rownum, Ind);
		addedrows.clear();
	}

	if (Ind) delete[] Ind;
}

bool OSISolver::lastpoint_feasible() {
	double tol; assert(osisolver->getDblParam(OsiPrimalTolerance, tol));

	// check violation of rows
	const double* rowval=osisolver->getRowActivity();
	const double* rowlow=osisolver->getRowLower();
	const double* rowup=osisolver->getRowUpper();
	for (int i=0; i<nr_row(); i++)
		if (rowlow[i]>rowval[i]+tol || rowup[i]<rowval[i]-tol) return false;

	// check box violations
	const double* colval=osisolver->getColSolution();
	const double* collow=osisolver->getColLower();
	const double* colup=osisolver->getColUpper();
	for (int i=0; i<nr_col(); i++) {
		if (collow[i]>colval[i]+tol || colup[i]<colval[i]-tol) return false;
		if ((osisolver->isBinary(i) || osisolver->isInteger(i)) && fabs(round(colval[i])-round(colval[i]))>tol) return false;
	}
	return true;
}

MIPSolver::SolutionStatus OSISolver::solve() {
#ifndef CPLEX_AVAILABLE
	osisolver->messageHandler()->setLogLevel(0);
#endif
	if (cold_start) {
		osisolver->initialSolve();
		cold_start=false;
//		osisolver->resolve();
/*#ifdef CPLEX_AVAILABLE
		// OsiCPLEX initial solve (primal cplex) cannot distinguish between UNBOUNDED and INFEASIBLE, calling resolve (dual simplex) comes to a decision
		OsiCpxSolverInterface* osicpx=dynamic_cast<OsiCpxSolverInterface*>(&*osisolver);
		if (CPXgetstat(osicpx->getEnvironmentPtr(), osicpx->getLpPtr(OsiCpxSolverInterface::FREECACHED_RESULTS))==CPX_STAT_INForUNBD)
			osisolver->resolve();
#endif*/
	}
	else {
		osisolver->resolve();
#ifndef CPLEX_AVAILABLE
		if (osisolver->isProvenPrimalInfeasible() || osisolver->isProvenDualInfeasible()) {
			osisolver->initialSolve();
			if (osisolver->isProvenOptimal())
				out_log << "Resolve reported infeasible, initialSolve says feasible...";
		}
#endif
	}

// #ifdef CPLEX_AVAILABLE
// 	if (osisolver->isProvenPrimalInfeasible()) {
// 		OsiCpxSolverInterface* osicpx=dynamic_cast<OsiCpxSolverInterface*>(&*osisolver);
// 		CPXiiswrite(osicpx->getEnvironmentPtr(), osicpx->getLpPtr(), "infeas.iis");
// 	}
// #endif

	if (osisolver->isAbandoned()) { cold_start=true; return ABORTED; }
	if (osisolver->isProvenPrimalInfeasible()) return INFEASIBLE;
	if (osisolver->isProvenDualInfeasible()) return UNBOUNDED;
	if (osisolver->isPrimalObjectiveLimitReached() || osisolver->isDualObjectiveLimitReached()) return UNBOUNDED;
	if (osisolver->isIterationLimitReached()) return ITERATIONLIMITEXCEEDED;
	if (osisolver->isProvenOptimal()) return SOLVED;
	if (lastpoint_feasible()) return FEASIBLE;
	return UNKNOWN;
}

MIPSolver::SolutionStatus OSISolver::solve(const UserVector<double>& x) {
	assert(x.dim()==nr_col());
	osisolver->setColSolution((Pointer<double>)x);
	return solve();
}

//#include "coin/CbcModel.hpp"

MIPSolver::SolutionStatus OSISolver::solveMIP() {
	return UNKNOWN;
//	CbcModel cbc(*osisolver);
//	cbc.initialSolve();
//	
//	if (cbc.isInitialSolveAbandoned()) return ABORTED;
//	if (cbc.isInitialSolveProvenPrimalInfeasible()) return INFEASIBLE;
//	if (cbc.isInitialSolveProvenDualInfeasible()) return UNBOUNDED;
//	
////	double* lowerstore=new double[nr_col()];
////	double* upperstore=new double[nr_col()];
////	for (int i=0; i<nr_col(); ++i) {
////		lowerstore[i]=osisolver->getColLower()[i];
////		upperstore[i]=osisolver->getColUpper()[i];
////	}
//	cbc.branchAndBound();
////	cbc.cleanModel(lowerstore, upperstore);
//	
//	if (cbc.isAbandoned()) return ABORTED;
//	if (cbc.isProvenOptimal()) return SOLVED;
//	if (cbc.isProvenInfeasible()) return INFEASIBLE;
//	if (cbc.isNodeLimitReached() || cbc.isSecondsLimitReached()) return ITERATIONLIMITEXCEEDED;
//	if (cbc.isSolutionLimitReached()) return UNBOUNDED;
//	return UNKNOWN;
}
/*

MIPSolver::SolutionStatus OSISolver::solveMIP() {
//	osisolver->messageHandler()->setLogLevel(10);
	
	double* lowerstore=new double[nr_col()];
	double* upperstore=new double[nr_col()];
	for (int i=0; i<nr_col(); ++i) {
		lowerstore[i]=osisolver->getColLower()[i];
		upperstore[i]=osisolver->getColUpper()[i];
//		out_log << i << ": " << lowerstore[i] << '\t' << upperstore[i];
//		if (osisolver->isInteger(i)) out_log << "\t is integer";
//		out_log << endl;
	}
	try {
		osisolver->branchAndBound();
	} catch(CoinError e) {
		out_log << "Catched CoinError: "; e.print(); 
	}
	for (int i=0; i<nr_col(); ++i) 	osisolver->setColBounds(i, lowerstore[i],upperstore[i]);
	delete[] lowerstore;
	delete[] upperstore;
//	osisolver->writeMps("mip");
//	osisolver->messageHandler()->setLogLevel(0);
//out_log << osisolver->isProvenDualInfeasible() << osisolver->isPrimalObjectiveLimitReached() << osisolver->isDualObjectiveLimitReached() << endl;
	if (osisolver->isAbandoned()) return ABORTED;
	if (osisolver->isProvenOptimal()) return SOLVED;
	if (osisolver->isProvenPrimalInfeasible()) return INFEASIBLE;
	if (osisolver->isProvenDualInfeasible()) return UNBOUNDED;
	if (osisolver->isPrimalObjectiveLimitReached() || osisolver->isDualObjectiveLimitReached()) return UNBOUNDED;
	if (osisolver->isIterationLimitReached()) return ITERATIONLIMITEXCEEDED;
	if (lastpoint_feasible()) return FEASIBLE;
	return UNKNOWN;
}
*/
MIPSolver::SolutionStatus OSISolver::feasible() {
	// save original objective
	dvector orig_obj(osisolver->getObjCoefficients(), nr_col());

	// clear objective
	dvector zeroobj(nr_col());
	int* indices=new int[nr_col()];
	for (int i=nr_col()-1; i>=0; --i) indices[i]=i;
	osisolver->setObjCoeffSet(indices, indices+nr_col(), (Pointer<double>)zeroobj);
	osisolver->setDblParam(OsiObjOffset, 0.);

	MIPSolver::SolutionStatus status(solve());

	// restore original objective
	osisolver->setObjCoeffSet(indices, indices+nr_col(), (Pointer<double>)orig_obj);
	osisolver->setDblParam(OsiObjOffset, -obj_const);
	delete indices;

	switch (status) {
		case SOLVED:
		case FEASIBLE:
		case UNBOUNDED: // but how could this happen when the objective is 0?
			return FEASIBLE;
	}
	return status;
}

void OSISolver::get_primal(UserVector<double>& x) {
	assert(x.dim()<=nr_col());
	x.set(osisolver->getColSolution());
}

double OSISolver::get_primal(const MIPSolver::ColItem& colitem) {
	return osisolver->getColSolution()[((const ColItem&)colitem).index];
}

void OSISolver::get_dual(UserVector<double>& mu) {
	assert(mu.dim()<=nr_row());
	mu.set(osisolver->getRowPrice());
}

double OSISolver::get_dual(const MIPSolver::RowItem& rowitem) {
	return osisolver->getRowPrice()[((const RowItem&)rowitem).index];
}

void OSISolver::get_reducedcosts(UserVector<double>& rc) {
	assert(rc.dim()<=nr_col());
	rc.set(osisolver->getReducedCost());
}
double OSISolver::get_reducedcosts(const MIPSolver::ColItem& colitem) {
	return osisolver->getReducedCost()[((const ColItem&)colitem).index];
}

void OSISolver::get_rowactivity(UserVector<double>& rowact) {
	assert(rowact.dim()<=nr_row());
	rowact.set(osisolver->getRowActivity());
}

double OSISolver::get_rowactivity(const MIPSolver::RowItem& rowitem) {
	return osisolver->getRowActivity()[((const RowItem&)rowitem).index];
}

double OSISolver::get_optval() { return osisolver->getObjValue(); }

int OSISolver::get_iter() { return osisolver->getIterationCount(); }

void OSISolver::set_obj(const UserVector<double>& obj, double obj_const_) {
	assert(obj.dim()==nr_col());
	obj_const=obj_const_;

	int* indices=new int[nr_col()];
	for (int i=nr_col()-1; i>=0; --i) indices[i]=i;
	osisolver->setObjCoeffSet(indices, indices+nr_col(), (Pointer<double>)obj);
	osisolver->setDblParam(OsiObjOffset, -obj_const);
	delete[] indices;
}

void OSISolver::modify_obj(int i, double coeff) {
	osisolver->setObjCoeff(i, coeff);
}

const MIPSolver::RowItem* OSISolver::add_row(const UserVector<double>& row, double low, double up) {
	assert(row.dim()==nr_col());
	// to be improved
	osisolver->addRow(CoinPackedVector(row.dim(), (Pointer<double>)row), low==-INFINITY ? -osisolver->getInfinity() : low, up==INFINITY ? osisolver->getInfinity() : up);
	addedrows.push_back(RowItem(nr_row()-1));
	list<RowItem>::iterator it(addedrows.end()); --it;
	it->it=it;
	return &(*it);
}

const MIPSolver::RowItem* OSISolver::add_row(const UserVector<double>& row, const ivector& indices, double low, double up) {
	osisolver->addRow(indices.size(), (Pointer<int>)indices, (Pointer<double>)row, low==-INFINITY ? -osisolver->getInfinity() : low, up==INFINITY ? osisolver->getInfinity() : up);
	addedrows.push_back(RowItem(nr_row()-1));
	list<RowItem>::iterator it(addedrows.end()); --it;
	it->it=it;
	return &(*it);
}

void OSISolver::add_rows(const vector<pair<dvector, ivector> >& rows, const dvector& low, const dvector& up) {
	CoinPackedVectorBase** coin_rows=new CoinPackedVectorBase*[rows.size()];
	double* coin_lb=new double[rows.size()];
	double* coin_ub=new double[rows.size()];
	for (int i=0; i<rows.size(); ++i) {
		coin_rows[i]=new CoinPackedVector(rows[i].first.size(), (Pointer<int>)rows[i].second, (Pointer<double>)rows[i].first);
		coin_lb[i]=low[i]==-INFINITY ? -osisolver->getInfinity() : low[i];
		coin_ub[i]=up[i]==INFINITY ? osisolver->getInfinity() : up[i];
	}
	osisolver->addRows(rows.size(), coin_rows, coin_lb, coin_ub);
	delete[] coin_lb;
	delete[] coin_ub;
	for (int i=0; i<rows.size(); ++i)
		delete coin_rows[i];
	delete[] coin_rows;
}

void OSISolver::add_rows(list<const MIPSolver::RowItem*>& rowitems, const vector<pair<dvector, ivector> >& rows, const dvector& low, const dvector& up) {
	for (int i=0; i<rows.size(); ++i) {
		addedrows.push_back(RowItem(nr_row()-1+i));
		list<RowItem>::iterator it(addedrows.end()); --it;
		it->it=it;
		rowitems.push_back(&(*it));
	}
	add_rows(rows, low, up);	
}


void OSISolver::delete_rows(const list<const MIPSolver::RowItem*>& rowitems) {
	if (rowitems.empty()) return;
	int num=rowitems.size();
	int* RowInd=new int[num];
	int i=0;
	for (list<const MIPSolver::RowItem*>::const_iterator it(rowitems.begin()); it!=rowitems.end(); ++it, ++i)
		RowInd[i]=((const RowItem*)(*it))->index;
	osisolver->deleteRows(num, RowInd);
	delete[] RowInd;

	for (list<const MIPSolver::RowItem*>::const_iterator it(rowitems.begin()); it!=rowitems.end(); ++it) {
		list<RowItem>::iterator it2(((const RowItem*)(*it))->it);
		for (it2=addedrows.erase(it2); it2!=addedrows.end(); ++it2)
			--it2->index;
	}
}

void OSISolver::modify_row(const MIPSolver::RowItem& rowitem, double low, double up) {
	modify_row(((const RowItem&)rowitem).index, low, up);
}

void OSISolver::modify_row(int index, double low, double up) {
	osisolver->setRowBounds(index, low==-INFINITY ? -osisolver->getInfinity() : low, up==INFINITY ? osisolver->getInfinity() : up);
}

const MIPSolver::ColItem* OSISolver::add_col(double low, double up, MipProblem::VarType vartype) {
	osisolver->addCol(CoinPackedVector(), low==-INFINITY ? -osisolver->getInfinity() : low, up==INFINITY ? osisolver->getInfinity() : up, 0.);
	if (vartype==MipProblem::BINARY || vartype==MipProblem::INTEGER) osisolver->setInteger(nr_col()-1);

	addedcols.push_back(ColItem(nr_col()-1));
	list<ColItem>::iterator it(addedcols.end()); --it;
	it->it=it;
	return &(*it);
}

const MIPSolver::ColItem* OSISolver::add_col(const UserVector<double>& col, double objcoeff, double low, double up, MipProblem::VarType vartype) {
	osisolver->addCol(CoinPackedVector(col.dim(), (Pointer<double>)col),
		low==-INFINITY ? -osisolver->getInfinity() : low, up==INFINITY ? osisolver->getInfinity() : up, objcoeff);
	if (vartype==MipProblem::BINARY || vartype==MipProblem::INTEGER) osisolver->setInteger(nr_col()-1);

	addedcols.push_back(ColItem(nr_col()-1));
	list<ColItem>::iterator it(addedcols.end()); --it;
	it->it=it;
	return &(*it);
}

void OSISolver::add_cols(list<const MIPSolver::ColItem*>& colitems, const dvector& low, const dvector& up) {
	if (!low.dim()) return;
	int nrcol=nr_col();
	dvector mylow(low); dvector myup(up);
	for (int i=0; i<low.dim(); i++) {
		if (low(i)==-INFINITY) mylow[i]=-osisolver->getInfinity();
		if (up(i)==INFINITY) myup[i]=osisolver->getInfinity();
		addedcols.push_back(ColItem(nrcol+i));
		list<ColItem>::iterator it(addedcols.end()); --it;
		it->it=it;
		colitems.push_back(&(*it));
	}
	CoinPackedVectorBase** cols=new CoinPackedVectorBase*[low.dim()];
	cols[0]=new CoinPackedVector();
	for (int i=1; i<low.dim(); i++) cols[i]=cols[0];
	osisolver->addCols(low.dim(), cols, (Pointer<double>)mylow, (Pointer<double>)myup, Pointer<double>(dvector(low.dim())));
	delete cols[0];
	delete cols;
}

void OSISolver::add_cols(list<const MIPSolver::ColItem*>& colitems, vector<Pointer<UserVector<double> > >& cols, const vector<double>& objcoeff, const dvector& low, const dvector& up) {
	if (!low.dim()) return;
	dvector mylow(low), myup(up);
	int nrcol=nr_col();
	for (int i=0; i<mylow.dim(); ++i) {
		if (mylow(i)==-INFINITY) mylow[i]=-osisolver->getInfinity();
		if (myup(i)==INFINITY) myup[i]=osisolver->getInfinity();
		addedcols.push_back(ColItem(nrcol+i));
		list<ColItem>::iterator it(addedcols.end()); --it;
		it->it=it;
		colitems.push_back(&(*it));
	}
	CoinPackedVectorBase** mycols=new CoinPackedVectorBase*[low.dim()];
	double* myobjcoeff=new double[low.dim()];
	for (int i=0; i<low.dim(); i++) {
		mycols[i]=new CoinPackedVector(cols[i]->dim(), (Pointer<double>)*cols[i]);
		myobjcoeff[i]=objcoeff[i];
	}
	osisolver->addCols(low.dim(), mycols, (Pointer<double>)mylow, (Pointer<double>)myup, myobjcoeff);

	delete myobjcoeff;
	for (int i=0; i<low.dim(); i++) delete mycols[i];
	delete mycols;
}

void OSISolver::delete_cols(const list<const MIPSolver::ColItem*>& colitems) {
	int num=colitems.size(); // :-(
	if (!num) return;

	int* ColInd=new int[num];
	int i=0;
	list<ColItem>::iterator first(((const ColItem*)(*colitems.begin()))->it); // pointing to the column, which will be removed, and has the smallest index
	for (list<const MIPSolver::ColItem*>::const_iterator it(colitems.begin()); it!=colitems.end(); ++it, ++i) {
		ColInd[i]=((const ColItem*)(*it))->index;
		if (first->index>ColInd[i]) first=((const ColItem*)(*it))->it;
	}
	osisolver->deleteCols(num, ColInd);
	delete ColInd;

	// now, making first pointing to the column, before the first, that will be removed, if possible
	int nr=first->index; // index of first removed column
	if (first==addedcols.begin()) first=addedcols.end();
	else --first;

	for (list<const MIPSolver::ColItem*>::const_iterator it(colitems.begin()); it!=colitems.end(); ++it)
		addedcols.erase(((const ColItem*)(*it))->it);

	// now, make first pointing to the first column, which got a new index
	if (first==addedcols.end()) first=addedcols.begin();
	else ++first;
	for (; first!=addedcols.end(); ++first, ++nr) // update index numbers
		first->index=nr;
	assert(nr==nr_col());
/*
	for (list<const MIPSolver::ColItem*>::const_iterator it(colitems.begin()); it!=colitems.end(); ++it) {
		list<ColItem>::iterator it2(((const ColItem*)(*it))->it);
		for (it2=addedcols.erase(it2); it2!=addedcols.end(); ++it2)
			--it2->index;
	}
*/}

void OSISolver::modify_col(const MIPSolver::ColItem& colitem, double low, double up, MipProblem::VarType vartype) {
	modify_col(((const ColItem&)colitem).index, low, up, vartype);
}

void OSISolver::modify_col(const MIPSolver::ColItem& colitem_, const UserVector<double>& col, double objcoeff, double low, double up, MipProblem::VarType vartype) {
	const ColItem& colitem(*(const ColItem*)&colitem_);

	osisolver->deleteCols(1, &colitem.index);
	osisolver->addCol(CoinPackedVector(col.dim(), (Pointer<double>)col),
		low==-INFINITY ? -osisolver->getInfinity() : low, up==INFINITY ? osisolver->getInfinity() : up, objcoeff);
	if (vartype==MipProblem::BINARY || vartype==MipProblem::INTEGER) osisolver->setInteger(nr_col()-1);
	else osisolver->setContinuous(nr_col()-1);

	for (list<ColItem>::iterator it(colitem.it); it!=addedcols.end(); ++it) --it->index; // update index numbers

	addedcols.splice(addedcols.end(), addedcols, colitem.it); // move colitem.it to end
	addedcols.back().index=nr_col()-1;

	assert(&addedcols.back()==&colitem);
}

void OSISolver::modify_col(int index, double low, double up, MipProblem::VarType type) {
	assert(index<nr_col());
	osisolver->setColBounds(index, low==-INFINITY ? -osisolver->getInfinity() : low, up==INFINITY ? osisolver->getInfinity() : up);
	switch (type) {
		case MipProblem::CONTINUOUS:
			osisolver->setContinuous(index);
			break;
		case MipProblem::BINARY:
		case MipProblem::INTEGER:
			osisolver->setInteger(index);
			break;
		default: out_err << "OSISolver: Variable type unkown. Aborting." << endl; exit(-1);
	}
}

double OSISolver::get_collow(int index) {
	assert(index<nr_col());
	return osisolver->getColLower()[index];
}

double OSISolver::get_colup(int index) {
	assert(index<nr_col());
	return osisolver->getColUpper()[index];
}

// #include "CglGomory.hpp"
#include "CglMixedIntegerRounding/CglMixedIntegerRounding.hpp"
#include "cuts.h"

int OSISolver::generate_cuts(list<Pointer<SimpleCut> >& rowcuts) {
//	if (!cutgenerator) cutgenerator=new CglGomory();
/*	if (!cutgenerator)*/ cutgenerator=new CglMixedIntegerRounding();
//	CglGomory gomory;
//	CglMixedIntegerRounding mir;
	OsiCuts cuts;
	cutgenerator->generateCuts(*osisolver, cuts);
//	mir.generateCuts(*osisolver, cuts);
//	cuts.printCuts();
	int nr=0;
	for (int i=0; i<cuts.sizeRowCuts(); ++i) {
		OsiRowCut& cut(cuts.rowCut(i));
		CoinPackedVector row(cut.row()); row.sortIncrIndex();
		Pointer<UserVector<double> > coeff(new SparseVector<double>(nr_col(), row.getNumElements(), row.getIndices(), row.getElements()));
		*coeff*=.5;
		if (cut.ub()<osisolver->getInfinity()) {
			rowcuts.push_back(new SimpleCut(coeff, -cut.ub()));
			++nr;	
			if (cut.lb()>-osisolver->getInfinity()) {
				Pointer<UserVector<double> > coeff2(coeff->getcopy()); *coeff2*=-1.;
				rowcuts.push_back(new SimpleCut(coeff2, cut.lb()));
				++nr;
			}
		} else if (cut.lb()>-osisolver->getInfinity()) {
			*coeff*=-1.;
			rowcuts.push_back(new SimpleCut(coeff, cut.lb()));
			++nr;	
		}
	}
	if (cuts.sizeRowCuts() || cuts.sizeColCuts()) {
		out_log << "Got " << cuts.sizeRowCuts() << " row cuts and " << cuts.sizeColCuts() << " column cuts. "
			<< "Returned " << nr << " row cuts." << endl;
	}
	
//	OsiSolverInterface::ApplyCutsReturnCode ret=osisolver->applyCuts(cuts);
//	clog << "Applied " << ret.getNumApplied() << " cuts." << endl;
	
//	clog << ret.getNumIneffective() << ' ' << ret.getNumInconsistent() << ' ' << ret.getNumInconsistentWrtIntegerModel() 
//	<< ' ' << ret.getNumInfeasible() << ' ' << ret.getNumApplied() << endl;
	
//	out_log << solve() << '\t' << get_optval() << endl;
	return nr;
}

//void OSISolver::print() { osisolver->writeMps("lp"); }
