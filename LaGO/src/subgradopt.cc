// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Hernan Alperin, Stefan Vigerske

#include "subgradopt.h"

extern "C" {
#include "ranlib.h"
}

int SubGradOpt::solve() {
  dvector d(dim());
  dvector grad(dim());
  dvector old_grad(dim());

  int f_eval_count=0;
	int ret;

  double f_value_here;

  double lambda;
  for (iter_=0; iter_<iter_max; iter_++) { // iteration count
    ret=obj.valgrad(f_value_here, grad, sol_point);  //S check for return-code!
    f_eval_count++;
		if (ret) return ret;

    int ret=check(f_value_here); //S check return-code!
    if (ret) break;

    if (grad.sq_norm2() < rtol) break;
    d = grad;

    if (ConjGrad) {
      double alpha=0.;
      if (iter_ % dim() > 0) {  // minor iterations (up to dim() conjugates)
        alpha=grad.sq_norm2()/old_grad.sq_norm2();
        d+=alpha*d_old; // alter direction to make it conjugate with prevs.
      }
      d_old=d;
    }
    old_grad = grad;

    switch (LineSearch) {

    case 2: {
      // find the first increasing point on the direction $d$
      //
      double beta=2*start_beta; // since beta is divided by 2 before
      double f_value_there;

      do {
	beta/=2;
	f_value_there=obj.eval(sol_point + beta * d);
        f_eval_count++;

      } while (f_value_here > f_value_there && beta > rtol);

      lambda=beta;
      if (beta > rtol) start_beta = beta*elasticity;

    } break;

    case 4: {
      // find the max on the direction $d$ with binary search
      // for $\nabla f(x)^Td =0$.
      double beta    = start_beta;
/*S      double f_value_there = obj.eval(sol_point+beta*d);
      dvector grad_there(dim());
      obj.grad(grad_there,sol_point+beta*d);
*/
      double f_value_there;
      dvector grad_there(dim());
      ret=obj.valgrad(f_value_there, grad_there, sol_point+beta*d);
			if (ret) return ret;
/*S  could be (or: it is) faster and you get a return-code! */

      f_eval_count++;
      while (grad_there*d > 0) {
	beta*=2;
	obj.grad(grad_there,sol_point+beta*d);
        f_eval_count++;
      }
      // find a decreasing point on the direction $d$
      // the maximum is between 0 and beta
      double here   = 0;
      double there  = beta;
      double middle = beta/2;
      do {
        dvector new_point(sol_point+middle*d);
	for (int i=0; i<new_point.size(); i++)
           if (new_point(i)<lower_bound[i]) new_point[i]=lower_bound[i];
	obj.grad(grad_there,new_point);
        f_eval_count++;
	if (grad_there*d >= 0)
	  here = middle;
	else
	  there = middle;
	middle = here + (there - here)/2;

      } while (there-here > tol);
      lambda = middle;
    } break;

    default: lambda=-d*obj.HessMult(sol_point,d);
      if (fabs(lambda)>rtol) lambda=old_grad.sq_norm2()/lambda;

    }

/*    out_solver_log <<
      " iter = " << iter_ <<
      " val = " << f_value_here <<
      " evaluation count = " << f_eval_count << endl;
*/    //    out_solver_log << " point = " << sol_point;

    sol_point+=lambda*d;

    for (int i=0; i<sol_point.size(); i++)
      if (sol_point(i)<lower_bound[i]) sol_point[i]=lower_bound[i];



  }
  f_value_here = obj.eval(sol_point);
  f_eval_count++;

  out_solver_log << " opt = " << f_value_here << " point = " << sol_point;
  opt_val_=f_value_here;

//	if (iter()==iter_max) return 3;
	return 0;
}


int Random::solve() {

  setall((long)rand(),(long)rand());     // set random-seed

  for (iter_=0; iter_<iter_max; iter_++) {
    for (int j=0; j<sol_point.size(); j++) sol_point[j] = gennor(0., dispersion);

    double val=obj.eval(sol_point);
    check(val); //S check return code?!
  }

	return 0;
}


