// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "ampl.h"

#ifdef COIN_HAS_ASL
#include <setjmp.h>

ASL *asl;

Pointer<MinlpProblem> ampl::get_problem() {
  FILE *nl;
  asl=ASL_alloc(ASL_read_fgh);

  want_xpi0=3; // if a starting guess is provided in the ampl-file, i want it

  nl = jac0dim(stubfile, (fint)strlen(stubfile));
  fgh_read(nl, 0);
	
  asl->i.err_jmp_=new Jmp_buf;
  if (setjmp(asl->i.err_jmp_->jb)) {
  	out_err << "caught error in asl" << endl;
  	exit(1);
  }
/*  asl->i.err_jmp1_=new Jmp_buf;
  if (setjmp(asl->i.err_jmp1_->jb)) {
  	out_err << "caught other error in asl" << endl;
  	exit(1);
  }
*/
  Pointer<MinlpProblem> prob(new MinlpProblem(MinlpProblem::MINLP, obj_name(0)));  // using objective name as problem name; or would the filename be better?
/*
  out_log << "var: " << n_var << "; binary: " << nbv << "; other integer: " << niv << "; nonlinear: " << MAX(nlvc,nlvo) << endl;
  out_log << "obj: " << n_obj << endl;
  out_log << "con: " << n_con << " nonlinear: " << nlc << endl;
*/
  int n=0;  // variable-counter
	
	for (int i=0; i<n_var; i++) {
		if (LUv[2*i]<=-1e+99) LUv[2*i]=-INFINITY;
		if (LUv[2*i]+1>=1e+99) LUv[2*i+1]=INFINITY;
	}
	
	bool integers=false; // indicates, whether we have integer-variables
  // nonlinear variables in constraints and objective
  for (int i=0; i<nlvb; i++, n++) {
		if (i>=nlvb-nlvbi && (LUv[2*n+1]-LUv[2*n])>1+rtol) integers=true;
    prob->add_var(n, 0, (i>=nlvb-nlvbi), LUv[2*n], LUv[2*n+1], var_name(n));
	}
  // nonlinear variables only in constraints
  for (int i=0; i<nlvc-nlvb; i++, n++) {
		if (i>=nlvc-nlvb-nlvci && (LUv[2*n+1]-LUv[2*n])>1+rtol) integers=true;
    prob->add_var(n, 0, (i>=nlvc-nlvb-nlvci), LUv[2*n], LUv[2*n+1], var_name(n));
	}
	
  // nonlinear variables only in objective
  for (int i=0; i<nlvo-nlvc; i++, n++) {
		if (i>=nlvo-nlvc-nlvoi && (LUv[2*n+1]-LUv[2*n])>1+rtol) integers=true;
    prob->add_var(n, 0, (i>=nlvo-nlvc-nlvoi), LUv[2*n], LUv[2*n+1], var_name(n));
	}

  // linear variables
  for (int i=0; i<n_var-MAX(nlvo, nlvc); i++, n++) {
		if (i>=n_var-MAX(nlvo, nlvc)-nbv-niv && (LUv[2*n+1]-LUv[2*n])>1+rtol) integers=true;
    prob->add_var(n, 0, (i>=n_var-MAX(nlvo, nlvc)-nbv-niv), LUv[2*n], LUv[2*n+1], var_name(n));
	}

	if (integers) {
		out_err << "Integer variables in problem. Not supported yet, aborting." << endl;
		exit(-1);
	}

  Pointer<SepQcFunc> func;
  Pointer<UserVector<double> > grad(new SparseVector<double>(n_var));  // for linear objectives or constraints

  // objective
  if (n_obj!=1) {
    if (n_obj<1) out_err << "Objective missing." << endl;
    else out_err << "Too many objective functions." << endl;
    exit(-1);
  }


	func=new SepQcFunc(prob->block);
	if (nlo) { // nonlinear objective
		Pointer<SparsityInfo> sparsity=new SparsityInfo();
		sparsity->linear=new map<int, SparsityInfo::LinearVariable>;
		sparsity->nonlinear=new map<int, SparsityInfo::NonlinearVariable>;
  		for (ograd *og=Ograd[0]; og; og=og->next)
			sparsity->nonlinear->insert(pair<int, SparsityInfo::NonlinearVariable>(og->varno, SparsityInfo::NonlinearVariable()));
		func->s[0]=new amplObj(prob->dim(), sparsity); // if nonlinear objectives
		func->set_curvature(0, Func::UNKNOWN);
  	} else {  // linear objective
		for (ograd *og=Ograd[0]; og; og=og->next) grad->SetElement(og->varno, 0.5*og->coef);
		func->b[0]=grad;
		func->c=objconst(0);
		func->set_curvature(0, Func::LINEAR);
	}
	prob->add_obj(func);

  // constraints
  cgrad *cg;
  for (int i=0; i<n_con; i++) {
    if (LUrhs[2*i+1]>INFINITY && LUrhs[2*i]<-INFINITY) {
      out_err << "Error: unbounded constraint " << con_name(i) << ": -inf <= f <= +inf" << endl;
      exit(-1);
    }

	Pointer<SparsityInfo> sparsity;

    if (i>=nlc) { // linear
		grad=new SparseVector<double>(n_var);
		for (cg=Cgrad[i]; cg; cg=cg->next)
      		grad->SetElement(cg->varno, 0.5*cg->coef);
    } else { // nonlinear
		sparsity=new SparsityInfo();
		sparsity->linear=new map<int, SparsityInfo::LinearVariable>;
		sparsity->nonlinear=new map<int, SparsityInfo::NonlinearVariable>;
		for (cg=Cgrad[i]; cg; cg=cg->next) // cannot distinguish between linear and nonlinear variables :-(
			sparsity->nonlinear->insert(pair<int, SparsityInfo::NonlinearVariable>(cg->varno, SparsityInfo::NonlinearVariable()));
    }

    if (fabs(LUrhs[2*i+1]-LUrhs[2*i])<rtol) { // constraint of form a <= f <= a -> f - a = 0
		func=new SepQcFunc(prob->block);
		if (i<nlc) {
			func->s[0]=new amplCon(i, prob->dim(), sparsity);
			func->set_curvature(0, Func::UNKNOWN);
		} else {
			func->b[0]=grad;
			func->set_curvature(0, Func::LINEAR);
		}
		func->c=-LUrhs[2*i+1];
		prob->add_con(func, true, con_name(i));
//      out_log << "added con " << i << ": " << con_name(i) << " = " << LUrhs[2*i] << endl;
    }
    else { // constraint of form a <= f <= b, a!=b
      if (LUrhs[2*i]>-INFINITY) {   // a <= f  -> -f + a <=0
        func=new SepQcFunc(prob->block);
        if (i<nlc) {
        	func->s[0]=new MinusFunc(new amplCon(i, prob->dim(), sparsity));
			func->set_curvature(0, Func::UNKNOWN);
        } else {
        	func->b[0]=new dvector(*grad * (-1));
			func->set_curvature(0, Func::LINEAR);
        }
        func->c=LUrhs[2*i];
        prob->add_con(func, false, con_name(i));
//        out_log << "added con " << i << ": " << LUrhs[2*i] << " <= " << con_name(i) << endl;
      }
      if (LUrhs[2*i+1]<INFINITY) {   // f <= b -> f - b <= 0
        func=new SepQcFunc(prob->block);
        if (i<nlc) {
			func->s[0]=new amplCon(i, prob->dim(), sparsity);
			func->set_curvature(0, Func::UNKNOWN);
        } else {
        	func->b[0]=grad;
			func->set_curvature(0, Func::LINEAR);
        }
        func->c=-LUrhs[2*i+1];
        prob->add_con(func, false, con_name(i));
//        out_log << "added con " << i << ": " << con_name(i) << " <= " << LUrhs[2*i+1] << endl;
      }
    }
  }

  // starting guess
  if (X0) prob->primal_point.set(X0, n_var);
  else
		for (int i=0; i<prob->primal_point.dim(); i++)
			if (prob->lower(i)<=-INFINITY)
				if (prob->upper(i)>=INFINITY)
					prob->primal_point[i]=random(prob->lower(i), prob->upper(i));
				else prob->primal_point[i]=prob->upper(i);
			else
				if (prob->upper(i)>=INFINITY)
					prob->primal_point[i]=random(prob->lower(i), prob->upper(i));
				else prob->primal_point[i]=.5*(prob->lower(i)+prob->upper(i));

  return prob;
};

#endif // COIN_HAS_ASL
