// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

// snopt.cc
// interface for SNOPT
//
// C++ front end to SNOPT

#include "snopt.h"
#include <sys/stat.h>
#include "f77filehandle.h"

#ifdef NOUNDERSCORE
#define SNINIT sninit
#define SNOPT snopt
#else
#define SNINIT sninit_
#define SNOPT snopt_
#endif

bool snoptlicenceok=false;

/** Initialize SNOPT.
    @param specsfile The name of a SPECS-file with SNOPT-options.
    @param printfile The name of the file to open for SNOPT-print-messages.
    @param summaryfile The name of the file to open for SNOPT-summary-messages.
    @param inform Hold's exit-level when returns.
    @param cw Character-workspace for SNOPT.
    @param lencw Size of cw.
    @param iw Integer-workspace for SNOPT.
    @param leniw Size of iw.
    @param rw Double-workspace for SNOPT.
    @param lenrw Size of rw.
*/
extern "C" void SNINIT(
      int& Specsfileid, int& Printfileid, int& Summaryfileid,
      int& inform,
      char* cw, int &lencw,
      int* iw, int &leniw,
      double* rw, int &lenrw);

/** The SNOPT-solver.
    @param ... A lot, I don't want to describe them all here. You should have a look at the private-part of snopt.h for it.
*/
extern "C" void SNOPT(
      char* start, int& m, int& n, int& ne, int& nName,
      int& nnCon, int& nnObj, int& nnJac,
      int& iObj, double& ObjAdd, char* Prob,
      void (*fcon)(int&, int&, int&, int&, double*, double*, double*, int&, char*, int&, int*, int&, double*, int&),
      void (*fobj)(int&, int&, double*, double&, double*, int&, char*, int&, int*, int&, double*, int&),
      double *a, int *ha, int *ka, double *bl, double *bu,
      char* Names,
      int* hs, double* xs, double *pi, double* rc,
      int& inform, int& mincw, int& miniw, int& minrw,
      int& nS, int& nInf, double& sInf, double& Obj,
      char *cu, int& lencu, int* iu, int& leniu, double* ru, int& lenru,
      char *cw, int& lencw, int* iw, int& leniw, double* rw, int& lenrw );

SnOpt::SnOpt(const Pointer<MinlpProblem> p_, Pointer<Param> param_, char* param_prefix, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: LocOpt(p_->dim(), out_solver_p_, out_solver_log_p_), minlp(p_), param(param_ ? param_ : Pointer<Param>(new Param())),
  param_prefix_end(param_prefix ? strlen(param_prefix)+1 : 0)
{ param_name=new char[param_prefix_end+100];
  if (param_prefix) {
  	strcpy((char*)param_name, param_prefix);
  	param_name[param_prefix_end-1]=' ';
  }
  param_name[param_prefix_end]=0;
  
	if (!snoptlicenceok) {
		out_err << "No valid SNOPT licence. Aborting." << endl;
		exit(-1);
	}

	init();
};

SnOpt::~SnOpt() {
  if (Specsfileid) closeF77file(Specsfileid);
  if (Printfileid) closeF77file(Printfileid);
  if (Summaryfileid) closeF77file(Summaryfileid);
}

void SnOpt::sort() {
  vector<bool> con_lin(minlp->con.size(), true);
	nonlin_con.clear();
	lin_con.clear();
	nonlin_block_conobj.clear();
	nonlin_block_cononly.clear();
	nonlin_block_objonly.clear();
	lin_block.clear();
	all_block.clear();

	for (int k=minlp->block.size()-1; k>=0; --k) {
		bool linear=true;
		for (int c=minlp->con.size()-1; c>=0; --c) {
			if (minlp->con[c]->A[k] || minlp->con[c]->s[k]) {
				linear=false;
				con_lin[c]=false;
			}
		}
		if (linear)
			if (minlp->obj->A[k] || minlp->obj->s[k]) nonlin_block_objonly.push_back(k);
			else lin_block.push_back(k);
		else
			if (minlp->obj->A[k] || minlp->obj->s[k]) nonlin_block_conobj.push_back(k);
			else nonlin_block_cononly.push_back(k);
	}
	if ((!nonlin_block_cononly.empty()) && (!nonlin_block_objonly.empty()))
		nonlin_block_conobj.splice(nonlin_block_conobj.end(), nonlin_block_cononly);

	all_block=nonlin_block_conobj;
	if (!nonlin_block_cononly.empty()) all_block.insert(all_block.end(), nonlin_block_cononly.begin(), nonlin_block_cononly.end());
	if (!nonlin_block_objonly.empty()) all_block.insert(all_block.end(), nonlin_block_objonly.begin(), nonlin_block_objonly.end());
	if (!lin_block.empty()) all_block.insert(all_block.end(), lin_block.begin(), lin_block.end());

	nnJac=0;
	nnObj=0;
	int size;
	for (list<int>::iterator it(nonlin_block_conobj.begin()); it!=nonlin_block_conobj.end(); ++it) {
		size=minlp->block[*it].size();
		nnJac+=size;
		nnObj+=size;
	}
	iObj=0;
	for (list<int>::iterator it(nonlin_block_cononly.begin()); it!=nonlin_block_cononly.end(); ++it) {
		nnJac+=minlp->block[*it].size();
		if ((!iObj) && minlp->obj->b[*it]) iObj=minlp->con.size()+1;
	}
	for (list<int>::iterator it(nonlin_block_objonly.begin()); it!=nonlin_block_objonly.end(); ++it)
		nnObj+=minlp->block[*it].size();
	if (!iObj)
		for (list<int>::iterator it(lin_block.begin()); it!=lin_block.end(); ++it)
			if (minlp->obj->b[*it]) {
				iObj=minlp->con.size()+1;
				break;
			}

	nnCon=0;
	for (int c=con_lin.size()-1; c>=0; --c)
		if (con_lin[c]) lin_con.push_back(c);
		else {
			nonlin_con.push_back(c);
			++nnCon;
		}

	// set number of constraints
  m = minlp->con.size()+(iObj ? 1 : 0);
  // set number of variables without slacks
  n = minlp->dim();

	if (out_solver_log_p) {
	  out_solver_log << "Number of constraints: " << m << endl;
	  out_solver_log << "Number of variables: " << n << endl;
		out_solver_log << "Block sortation: ";
		for (list<int>::iterator it(nonlin_block_conobj.begin()); it!=nonlin_block_conobj.end(); ++it) out_solver_log << *it << ' ';
		out_solver_log << '|' << ' ';
		for (list<int>::iterator it(nonlin_block_cononly.begin()); it!=nonlin_block_cononly.end(); ++it) out_solver_log << *it << ' ';
		out_solver_log << '|' << ' ';
		for (list<int>::iterator it(nonlin_block_objonly.begin()); it!=nonlin_block_objonly.end(); ++it) out_solver_log << *it << ' ';
		out_solver_log << '|' << ' ';
		for (list<int>::iterator it(lin_block.begin()); it!=lin_block.end(); ++it) out_solver_log << *it << ' ';
		out_solver_log << endl;
		out_solver_log << "Constraint sortation: ";
		for (list<int>::iterator it(nonlin_con.begin()); it!=nonlin_con.end(); ++it) out_solver_log << *it << ' ';
		out_solver_log << '|' << ' ';
		for (list<int>::iterator it(lin_con.begin()); it!=lin_con.end(); ++it) out_solver_log << *it << ' ';
		out_solver_log << endl;

	  // number of nonlinear constraints
  	out_solver_log << "nnCon: " << nnCon << endl;
	  // number of nonlinear objective variables
  	out_solver_log << "nnObj: " << nnObj << endl;
  	// number of nonlinear variables in the jacobian-matrix
  	out_solver_log << "nnJac: " << nnJac << endl;
  	// says, which constraint holds the linear part of the objective
  	out_solver_log << "iObj: " << iObj << endl;
	}
}

void SnOpt::init_jacobian() {
	SparseMatrix jac(m, n);

	int varind=0;
	for (list<int>::iterator itb(nonlin_block_conobj.begin()); itb!=nonlin_block_conobj.end(); ++itb) {
		int c=0;
		for (list<int>::iterator itc(nonlin_con.begin()); itc!=nonlin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->A[*itb] || minlp->con[*itc]->b[*itb] || minlp->con[*itc]->s[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 1.);
		for (list<int>::iterator itc(lin_con.begin()); itc!=lin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->b[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 2.*(*minlp->con[*itc]->b[*itb])(i));
		varind+=minlp->block[*itb].size();
	}

	for (list<int>::iterator itb(nonlin_block_cononly.begin()); itb!=nonlin_block_cononly.end(); ++itb) {
		int c=0;
		for (list<int>::iterator itc(nonlin_con.begin()); itc!=nonlin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->A[*itb] || minlp->con[*itc]->b[*itb] || minlp->con[*itc]->s[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 1.);
		for (list<int>::iterator itc(lin_con.begin()); itc!=lin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->b[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 2.*(*minlp->con[*itc]->b[*itb])(i));
		if (iObj && minlp->obj->b[*itb])
			for (int i=0; i<minlp->block[*itb].size(); ++i)
				jac.AddElement(c, varind+i, 2.*(*minlp->obj->b[*itb])(i));
		varind+=minlp->block[*itb].size();
	}

	for (list<int>::iterator itb(nonlin_block_objonly.begin()); itb!=nonlin_block_objonly.end(); ++itb) {
		int c=0;
		for (list<int>::iterator itc(nonlin_con.begin()); itc!=nonlin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->b[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 2.*(*minlp->con[*itc]->b[*itb])(i));
		for (list<int>::iterator itc(lin_con.begin()); itc!=lin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->b[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 2.*(*minlp->con[*itc]->b[*itb])(i));
		varind+=minlp->block[*itb].size();
	}

	for (list<int>::iterator itb(lin_block.begin()); itb!=lin_block.end(); ++itb) {
		int c=0;
		for (list<int>::iterator itc(nonlin_con.begin()); itc!=nonlin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->b[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 2.*(*minlp->con[*itc]->b[*itb])(i));
		for (list<int>::iterator itc(lin_con.begin()); itc!=lin_con.end(); ++itc, ++c)
			if (minlp->con[*itc]->b[*itb])
				for (int i=0; i<minlp->block[*itb].size(); ++i)
					jac.AddElement(c, varind+i, 2.*(*minlp->con[*itc]->b[*itb])(i));
		if (iObj && minlp->obj->b[*itb])
			for (int i=0; i<minlp->block[*itb].size(); ++i)
				jac.AddElement(c, varind+i, 2.*(*minlp->obj->b[*itb])(i));
		varind+=minlp->block[*itb].size();
	}
	assert(varind==n);

	jac.finish();

	ne=jac.nonzeros();
  a  = new double[ne]; const double* val=jac.GetVal();
  ha = new int[ne]; const int* rowind=jac.GetRowInd();
  for (int i=0; i<ne; ++i, ++val, ++rowind) {
    a[i] = *val;
    ha[i] = *rowind+1;
  }
  ka = new int[n+1]; const int* colptr=jac.GetColPtr();
  for (int i=0; i<n+1; ++i, ++colptr) ka[i] = *colptr+1;
//out_log << dvector(a, ne) << ivector(ha, ne) << ivector(ka, n+1);
  // set maximal number of nonzero-elements
  out_solver_log << "Number of nonzero's in the jacobian: " << ne << endl;
}

void SnOpt::init_bounds() {
  // the lower and upper bounds of the continuous variables and slack-variables
  out_solver_log << "Set lower and upper bounds of variables" << endl;

	int varind=0;
	int i0, j;
	for (list<int>::iterator it(all_block.begin()); it!=all_block.end(); ++it) {
		for (j=0; j<minlp->block[*it].size(); ++j, ++varind) {
			i0=minlp->block[*it][j];
			bl[varind]=minlp->lower[i0];
			if (bl[varind]==-INFINITY) bl[varind]=-1.E+20;
			bu[varind]=minlp->upper[i0];
			if (bu[varind]==INFINITY) bu[varind]=1.E+20;
		}
	}

  // if dual bound is 0, it's an inequality-constraint, if it's not 0 (-INFINITY), it's an equation
  // so lower-bound is -infinity for inequalities and the constant part for equations
  // upper-bound is minus the constant part minus the value of the discrete variables, if not in upper-left part
  out_solver_log << "Set lower and upper bounds of constraints" << endl;
	int c=0;
	for (list<int>::iterator itc(nonlin_con.begin()); itc!=nonlin_con.end(); ++itc, ++c) {
		bu[n+c]=-minlp->con[*itc]->c;
    bl[n+c]=(minlp->con_eq[*itc] ? bu[n+c] : -1.E+20); // -10^20 was recommended in the snopt-manual
	}
	for (list<int>::iterator itc(lin_con.begin()); itc!=lin_con.end(); ++itc, ++c) {
		bu[n+c]=-minlp->con[*itc]->c;
    bl[n+c]=(minlp->con_eq[*itc] ? bu[n+c] : -1.E+20); // -10^20 was recommended in the snopt-manual
	}
	if (iObj) {
		bl[n+m-1]=-1.E+20;
		bu[n+m-1]=1.E+20;
	}

	if (out_solver_log_p) {
  	out_solver_log << "lower bounds of (slack-)var.: "; for (int i=0; i<n+m; i++) out_solver_log << bl[i] << " "; out_solver_log << endl;
  	out_solver_log << "upper bounds of (slack-)var.: "; for (int i=0; i<n+m; i++) out_solver_log << bu[i] << " "; out_solver_log << endl;
	}
}

void SnOpt::init() {
	sort();

  lencw=500;
  cw=new char[lencw*8];

  leniw=100*MAX(5, (n+m));
  iw = new int[leniw];

  lenrw=2*leniw;
  rw = new double[lenrw];

	strcpy(param_name+param_prefix_end, "snopt specs");
  Specsfile=param->get(param_name);
	if (Specsfile) {
		Specsfileid=openF77file(Specsfile);
	  if (Specsfileid<0) {
	  	out_err << "SnOpt: Error opening " << Specsfile << ", maybe no free filedescriptor left." << endl;
	  	Specsfileid=0;
	  }
  } else Specsfileid=0;
  out_solver_log << "Specsfile: " << (Specsfile ? (char*)Specsfile : "none") << endl;

	strcpy(param_name+param_prefix_end, "snopt print");
  Printfile=param->get(param_name);
	if (Printfile) {
	  Printfileid=openF77file(Printfile);
	  if (Printfileid<0) {
	  	out_err << "SnOpt: Error opening " << Printfile << ", maybe no free filedescriptor left." << endl;
	  	Printfileid=0;
	  }
  } else Printfileid=0;

  out_solver_log << "Printfile: " << (Printfile ? (char*)Printfile : "none") << endl;

	strcpy(param_name+param_prefix_end, "snopt summary");
  Summaryfile=param->get(param_name);
	if (Summaryfile) {
		Summaryfileid=openF77file(Summaryfile);
	  if (Summaryfileid<0) {
	  	out_err << "SnOpt: Error opening " << Summaryfile << ", maybe no free filedescriptor left." << endl;
	  	Summaryfileid=0;
	  }
  } else Summaryfileid=0;
  out_solver_log << "Summaryfile: " << (Summaryfile ? (char*)Summaryfile : "none") << endl;

	int inform=0;
  SNINIT(Specsfileid, Printfileid, Summaryfileid,
         inform,
         cw,lencw,
         iw,leniw,
         rw,lenrw);

  out_solver_log << "SNINIT called: " << inform << " " << endl;

  start = strdup("Cold");

	init_jacobian();

	// constant part of the objective
  ObjAdd=minlp->obj->c;
  out_solver_log << "ObjAdd: " << ObjAdd << endl;

	bl=new double[n+m];
  bu=new double[n+m];
	init_bounds();

  xs = new double[n+m];
  hs = new int[n+m];
  pi=new double[m];

  for (int i=0; i<n+m; ++i) {
  	hs[i]=4;
  	xs[i]=bl[i];
  }

  rc=new double[n+m];

  out_solver_log << "SNOPT::Init() finished." << endl;
}

void SnOpt::reinit() {
	init_bounds();

  for (int i=0; i<n+m; ++i) {
  	hs[i]=4;
  	xs[i]=bl[i];
  }

  out_solver_log << "SNOPT::reinit() finished." << endl;
};

bool SnOpt::nnobj(double& val, double* grad, int mode, const double* x) {
	val=0.;
	int size;
	int i, k;

#define nnobjblock(kk) { \
		k=kk; \
		size=minlp->block[k].size(); \
		if ((!minlp->obj->A[k]) && (!minlp->obj->b[k]) && (!minlp->obj->s[k])) { \
			x+=size; \
			if (mode) for (i=0; i<size; ++i) *(grad++)=0.; \
			continue; \
		} \
		dvector xb(x, size); \
		if (mode==0 || mode==2) val+=minlp->obj->eval(xb, k); \
		if (mode) { \
			dvector g(size); \
			minlp->obj->grad(g, xb, k); \
			for (i=0; i<size; ++i) *(grad++)=g(i); \
		} \
		x+=size; \
	};

	list<int>::iterator it(nonlin_block_conobj.begin());
	for (; it!=nonlin_block_conobj.end(); ++it) nnobjblock(*it);
	for (it=nonlin_block_objonly.begin(); it!=nonlin_block_objonly.end(); ++it) nnobjblock(*it);
#undef nnobjblock

	if (!finite(val)) { out_err << "SnOpt::nnobj: infinite functionvalue" << endl; return false; }
	if (mode) for (int i=nnObj; i>0; --i) if (!finite(*(--grad))) { out_err << "SnOpt::nnobj: infinite gradient" << endl; return false; }
	return true;
}

extern "C" void FUNOBJ (
      int& mode, int& nnObj,
      double* x, double& fObj, double* gObj, int& nState,
      char* cu, int& lencu,
      int* iu, int& leniu,
      double* ru, int& lenru )
{ if (!((SnOpt*)cu)->nnobj(fObj, gObj, mode, x)) nState=-1;
}

bool SnOpt::nncon(double* val, double* grad, int mode, const double* x) {
	int size;
	int c=0;
	int k;
	while (c<nnCon) val[c++]=0.;

	vector<dvector> congrad(mode ? nnCon : 0);

#define	nnconblock(kk) { \
    k=kk; \
		size=minlp->block[k].size(); \
		dvector xb(x, size); \
		c=0; \
		for (list<int>::iterator itc(nonlin_con.begin()); itc!=nonlin_con.end(); ++itc, ++c) { \
			if ((!minlp->con[*itc]->A[k]) && (!minlp->con[*itc]->b[k]) && (!minlp->con[*itc]->s[k])) { \
				if (mode) congrad[c].resize(0); \
				continue; \
			} \
			if (mode==0 || mode==2) val[c]+=minlp->con[*itc]->eval(xb, k); \
			if (mode) { \
				congrad[c].resize(size); \
				minlp->con[*itc]->grad(congrad[c], xb, k); \
			} \
		} \
		if (mode) \
			for (int i=0; i<size; ++i) \
				for (c=0; c<nnCon; ++c) \
					if (congrad[c].size()) *(grad++)=congrad[c](i); \
		x+=size; \
	}

	list<int>::iterator itb(nonlin_block_conobj.begin());
	for (; itb!=nonlin_block_conobj.end(); ++itb) nnconblock(*itb);
	for (itb=nonlin_block_cononly.begin(); itb!=nonlin_block_cononly.end(); ++itb) nnconblock(*itb);
#undef nnconblock

	for (c=nnCon; c; --c) if (!finite(*(val++))) { out_err << "SnOpt::nncon: infinite functionvalue" << endl; return false; }
	if (mode) for (c=nnCon*nnJac; c>0; --c) if (!finite(*(--grad))) { out_err << "SnOpt::nncon: infinite gradient" << endl; return false; }
	return true;
}

extern "C" void FUNCON (
      int& mode, int& nnCon, int& nnJac, int& neJac,
      double* x, double* fCon,
      double* gCon,
      int& nState,
      char* cu, int& lencu,
      int* iu, int& leniu,
      double* ru, int& lenru )
{ if (!((SnOpt*)cu)->nncon(fCon, gCon, mode, x)) nState=-1;
}

int SnOpt::solve() {
	timer.start();

	if (iter_max<INF) iw[88]=iter_max;
	iw[77]=-1; // no check on funcon and funobj

	int nS; // Number of superbasic variables.
	int nInf; // Number of infeasible constraints.
	double sInf;
	int inform;
  int mincw,miniw,minrw; //Info about required character/integer/real snopt workspace, if inform between 42 and 44.
	int lenzero=0; // Length for user workspaces: 0.
	int nName=1; //Number of col&row name in char array Names.
	char Prob[8]={' ',' ',' ',' ',' ',' ',' ',' '};

  SNOPT(start, m, n, ne, nName,
	 nnCon, nnObj, nnJac,
	 iObj, ObjAdd, Prob,
	 &FUNCON, &FUNOBJ,
	 a, ha, ka, bl, bu,
	 NULL /* Names */,
	 hs, xs, pi, rc,
	 inform, mincw, miniw, minrw,
	 nS, nInf, sInf, opt_val_,
	 (char*)this, lenzero, NULL, lenzero, NULL, lenzero,
	 cw, lencw, iw, leniw, rw, lenrw);

  if ((inform>=42 && inform <=44) || inform==20) {
    out_solver_log << "SnOpt: Too less workspace for SNOPT (inform=" << inform << "). Trying to increase it." << endl;
    if (mincw>lencw) {
      Pointer<char> oldcw=cw;
      cw=new char[mincw];
      for (int i=0; i<MIN(lencw, mincw); i++) cw[i]=oldcw[i];
      lencw=mincw;
      out_solver_log << "SnOpt: cw length increased to " << lencw << endl;
    }

    if (miniw>leniw) {
      Pointer<int> oldiw=iw;
      miniw=(int)(1.5*miniw);
      iw=new int[miniw];
      for (int i=0; i<MIN(leniw, miniw); i++) iw[i]=oldiw[i];
      leniw=miniw;
      out_solver_log << "SnOpt: iw length increased to " << leniw << endl;
    }

    if (minrw>lenrw) {
      Pointer<double> oldrw=rw;
      minrw=(int)(1.5*minrw);
      rw=new double[minrw];
      for (int i=0; i<MIN(lenrw, minrw); i++) rw[i]=oldrw[i];
      lenrw=minrw;
      out_solver_log << "SnOpt: rw length increased to " << lenrw << endl;
    }

    SNOPT(start, m, n, ne, nName,
          nnCon, nnObj, nnJac,
          iObj, ObjAdd, Prob,
          &FUNCON, &FUNOBJ,
          a, ha, ka, bl, bu,
          NULL,
          hs, xs, pi, rc,
          inform, mincw, miniw, minrw,
          nS, nInf, sInf, opt_val_,
		      (char*)this, lenzero, NULL, lenzero, NULL, lenzero,
		      cw, lencw, iw, leniw, rw, lenrw);
  }

  out_solver << "SnOpt: solve finished: " << inform << endl;

	if (inform>5 && inform!=9)
		out_err << "Warning: SnOpt returned " << inform << endl;

	int j=0;
	int i,size;
	for (list<int>::iterator it(all_block.begin()); it!=all_block.end(); ++it) {
		size=minlp->block[*it].size();
		for (i=0; i<size; ++i, ++j)
			sol_point[minlp->block[*it][i]]=xs[j];
	}

  if (out_solver_log_p) {
	  out_solver_log << "SnOpt: solution-vector: " << endl;
		for (i=0; i<sol_point.size(); i++)
    	out_solver_log << minlp->var_names[i] << "\t" << sol_point(i) << endl;
	}

	opt_val_=minlp->obj->eval(sol_point);
  out_solver << "SnOpt: The objective value is: " << opt_val_ << endl;

  out_solver << "SnOpt: Number of infeasible constraints: " << nInf << endl;
  out_solver << "SnOpt: Sum of infeasibilities: " << sInf << endl;

  start=strdup("Cold"); // next start could be cold

	time_=timer.stop();

  return inform;
};

int SnOpt::solve(dvector& x) {
	int varind=0;
	int i0,size;
	for (list<int>::iterator it(all_block.begin()); it!=all_block.end(); ++it) {
		size=minlp->block[*it].size();
		for (int i=0; i<size; ++i, ++varind) {
			i0=minlp->block[*it][i];
			xs[varind]=x(i0);
			if (minlp->discr[i0]) {
				bl[varind]=bu[varind]=x(i0);
				hs[varind]=3; // basic
			}
			if (fabs(xs[varind]-bl[varind])<rtol) hs[varind]=0; // nonbasic on lower bound
			else if (fabs(xs[varind]-bu[varind])<rtol) hs[varind]=1; // nonbasic on upper bound
			else hs[varind]=2; // superbasic
		}
	}

  start=strdup("Warm");

  return solve();
}

dvector SnOpt::get_lag_multipliers() {
	dvector mypi(minlp->con.size());
	int c=0;
	for (list<int>::iterator it(nonlin_con.begin()); it!=nonlin_con.end(); ++it, ++c) mypi[*it]=-pi[c];
	for (list<int>::iterator it(lin_con.begin()); it!=lin_con.end(); ++it, ++c) mypi[*it]=-pi[c];
	return mypi;
}


dvector SnOpt::get_mu_q() {
	dvector mu_q(n);

	int size, i;
	int j=0;
	for (list<int>::iterator it(all_block.begin()); it!=all_block.end(); ++it) {
		size=minlp->block[*it].size();
		for (i=0; i<size; ++i, ++j) {
			if (bu[j]-bl[j]<rtol) continue;
			if (fabs(xs[j]-bl[j])<rtol) mu_q[minlp->block[*it][i]]=rc[j]/(bu[j]-bl[j]);
			if (fabs(xs[j]-bu[j])<rtol) mu_q[minlp->block[*it][i]]=rc[j]/(bl[j]-bu[j]);
		}
	}

	return mu_q;
}
