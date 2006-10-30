// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Hernan Alperin, Stefan Vigerske

#ifndef SUBGRADOPT_H
#define SUBGRADOPT_H

#include "standard.h"
#include "opt.h"
#include "func.h"
#include "param.h"

/** Subgradient algorithm */

class SubGradOpt : public DualSolver {
  private:
    /** The type of the line search.
    */
    int LineSearch;

    double start_beta;
    double elasticity;

    /** Old direction.
    */
    dvector d_old;
    /** Indicates, whether to use the conjugate gradient method, or not.
    */
    bool ConjGrad;

  public:
    /** Constructor for a function, the lower bounds and parameters.
        @param f The function, to optimize.
        @param lower_bound_ The lower bounds of the variables.
        @param param Some parameters.
        @param out_solver_p_ A Pointer to an ostream to print solver relevant output to.
        @param out_solver_log_p A Pointer to an ostream to print solver relevant logging output to.
    */
    SubGradOpt(DualFunc& f, vector<double> lower_bound_, Param& param_, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p)
    : DualSolver(f, lower_bound_, param_, out_solver_p_, out_solver_log_p_),
      d_old(f.dim()),
      ConjGrad(param.get_i("Conjugate Gradient", 0)),
      LineSearch(param.get_i("Line Search", 0)),
      start_beta(param.get_d("start beta", 1.)),
      elasticity(param.get_d("elasticity", 4.))
    { iter_max=param.get_i("SubGradOpt max iter",100);
    	tol=param.get_d("SubGradOpt tol", 0.001);
    }

    int solve(dvector& z) {
      sol_point=z;
      return solve();
    }

    int solve();

};


class Random : public DualSolver {
  private:
    double dispersion;

  public:
    /** Constructor for a function, the lower bounds and parameters.
        @param f The function, to optimize.
        @param lower_bound_ The lower bounds of the variables.
        @param param Some parameters.
        @param out_solver_p_ A Pointer to an ostream to print solver relevant output to.
        @param out_solver_log_p A Pointer to an ostream to print solver relevant logging output to.
    */
    Random(DualFunc& f, vector<double> lower_bound_, Param& param_, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p)
    : DualSolver(f, lower_bound_, param_, out_solver_p_, out_solver_log_p_),
      dispersion(param_.get_d("dispersion", 1.))
      { iter_max=param_.get_i("Random max iter",1);
      }
    
    int solve(dvector& z) {
      sol_point=z;
      return solve();
    }

    int solve();

};

#endif
