// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

// snopt.h
// C++ interface to SNOPT

#ifdef SNOPT_AVAILABLE
#ifndef SNOPT_H
#define SNOPT_H

#include "standard.h"
#include "opt.h"
#include "problem.h"
#include "param.h"
#include "usermatrix.h"

#ifdef NOUNDERSCORE
#define FUNCON funcon
#define FUNOBJ funobj
#else
#define FUNCON funcon_
#define FUNOBJ funobj_
#endif

extern bool snoptlicenceok;

/** Subroutine, that computes the obj fn and its gradient.
    This function is called by SNOPT.
    @param mode The mode: 0 to evaluate, 1 to compute the gradient and 2 to evaluate and compute
    @param nnObj The number of nonlinear variables in the objective.
    @param x The vector to compute for. The size is nnObj.
    @param fObj The nonlinear value of the nonlinear objective-variables, only set in mode 0 or 2.
    @param gObj The nonlinear gradient of the nonlinear objective-variables, only set in mode 1 or 2.
    @param nState To return a status.
    @param cu Used as a point to the SnOpt-object.
    @param lencu The size of cu. It's set to 0 by Init().
    @param iu Integer-Workspace, not used.
    @param leniu The size of iu. It's set to 0 by Init().
    @param ru Real-Workspace, not used.
    @param lenru The size of ru. It's set to 0 by Init().
*/
extern "C" void FUNOBJ (
      int& mode, int& nnObj,
      double* x, double& fObj, double* gObj, int& nState,
      char* cu, int& lencu,
      int* iu, int& leniu,
      double* ru, int& lenru );

/** The function, that computes the vector of nonlinear constraints and its Jacobian.
    This function is called by SNOPT.
    @param mode The mode: 0 to evaluate, 1 to compute the gradient or 2 to evaluate and compute.
    @param nnCon The number of nonlinear constraints.
    @param nnJac The number of nonlinear variables.
    @param neJac The number of nonzero elements in the nonlinear jacobian. In a dense case, it's nnCon*nnJac.
    @param x The vector to compute, size is nnJac.
    @param fCon The array of size nnCon to store the evaluated constraints in.
    @param gCon The array of size neJac to store the gradient in (column-wise).
    @param nState Not used.
    @param cu Used as a point to the SnOpt-object.
    @param lencu The size of cu. It's set to 0 by Init().
    @param iu Integer-Workspace, not used.
    @param leniu The size of iu. It's set to 0 by Init().
    @param ru Real-Workspace, not used.
    @param lenru The size of ru. It's set to 0 by Init().
*/
extern "C" void FUNCON (
      int& mode, int& nnCon, int& nnJac, int& neJac,
      double* x, double* fCon,
      double* gCon, int &nState,
      char* cu, int& lencu,
      int* iu, int& leniu,
      double* ru, int& lenru );

/** The interface to the SnOpt-solver.
    @class SnOpt
		@param snopt specs
		%options filename
		%default empty
		%level 1
		Name of a parameterfile for SNOPT.
		@param snopt print
		%options filename
		%default empty
		Name of a file for the printing-output of SNOPT.
		@param snopt summary
		%options filename
		%default empty
		Name of a file for the summary-output of SNOPT.
*/
class SnOpt : public LocOpt {
	friend void FUNOBJ (
      int& mode, int& nnObj,
      double* x, double& fObj, double* gObj, int& nState,
      char* cu, int& lencu,
      int* iu, int& leniu,
      double* ru, int& lenru );
	friend void FUNCON (
      int& mode, int& nnCon, int& nnJac, int& neJac,
      double* x, double* fCon,
      double* gCon, int &nState,
      char* cu, int& lencu,
      int* iu, int& leniu,
      double* ru, int& lenru );
  private:
    /** The optimization problem.
    */
    const Pointer<MinlpProblem> minlp;

    /** The parameters.
    */
    Pointer<Param> param;

    Pointer<char> param_name;
    int param_prefix_end;

		/** The indices of the nonlinear constraints.
		*/
		list<int> nonlin_con;
		/** The indices of the linear constraints.
		*/
		list<int> lin_con;
		/** The indices of the blocks, which involve nonlinear variables in constraints and objective.
		*/
		list<int> nonlin_block_conobj;
		/** The indices of the blocks, which involve nonlinear variables only in the constraints.
		*/
		list<int> nonlin_block_cononly;
		/** The indices of the blocks, which involve nonlinear variables only in the objective.
		*/
		list<int> nonlin_block_objonly;
		/** The indices of the blocks which have exclusively linear variables in constraints and objective.
		*/
		list<int> lin_block;

		/** A list of all blocks.
		    The concenation of nonlin_blcok_conobj, nonlin_block_cononly, nonlin_block_objonly and lin_block.
		*/
		list<int> all_block;

    /** The name of the SPECS-file, which contains options for the SNOPT-solver.
    */
    Pointer<char> Specsfile;
    int Specsfileid; // the fortran-id of the specsfile

    /** The name of the file for print messages.
    */
    Pointer<char> Printfile;
    int Printfileid;  // the fortran-id of the printfile

    /** The name of the file for summary messages.
    */
    Pointer<char> Summaryfile;
    int Summaryfileid; // the fortran-id of the summaryfile

    /** How to obtain a starting basis.
        "Cold", "Basis file" or "Warm".
    */
    Pointer<char> start;

    /* Number of constraints.
    */
    int m;

    /* Number of variables excluding slacks (n>0)
    */
    int n;

    /** The maximal number of non-zero elements in the Jacobian-matrix.
    */
    int ne;

    /** The index of the constraint, where the linear constraints start.
        Or the number of non-linear constraints.
    */
    int nnCon;

    /** The number of nonlinear variables in the objective.
    */
    int nnObj;

    /** The number of nonlinear variables in the constraints.
        It's the number of variables in the block's of nonlin_block
    */
    int nnJac;

    /** The index of the constraint, which is the linear part in the objective.
        In fortran-style.
    */
    int iObj;

    /** Additive constant to obj. for reporting.
    */
    double ObjAdd;

    /** array with the $ne$ elems of A
    */
    Pointer<double> a;
    /** Row indices.
    */
    Pointer<int> ha;
    /** pointer to the begin of each column $n+1$ elems.
        must be $ka(0) = 1$ and $ka(n+1) = ne+1$.
    */
    Pointer<int> ka;

    /** Vector with lower bounds for var and slacks dim $n+m$.
    */
    Pointer<double> bl;
    /** Vector with upper bounds for var and slacks dim $n+m$
    */
    Pointer<double> bu;

    /** Sometimes the initial state for vars and slacks values in ${0,\dots,5}$.
    */
    Pointer<int> hs;
    /** Sometimes the initial value for vars and slacks.
        Depends on hs.
    */
    Pointer<double> xs;

    /** Estimates of dual Lagr. multipliers.
    */
    Pointer<double> pi;

    /** Length for snopt character workspace.
    */
    int lencw;
    /** Length for snopt int workspace.
    */
    int leniw;
    /** Length for snopt real workspace.
    */
    int lenrw;
    /** Snopt character workspace.
    */
    Pointer<char> cw;
    /** Snopt int workspace.
    */
    Pointer<int> iw;
    /** Snopt real workspace.
    */
    Pointer<double> rw;

    /** vector for reduced costs $n+m$.
    */
    Pointer<double> rc;

    /** Sort's the variables and constraints by linearity.
    */
    void sort();

		void init_jacobian();

		void init_bounds();

		void init();

		bool nnobj(double& val, double* grad, int mode, const double* x);

		bool nncon(double* val, double* grad, int mode, const double* x);

  public:
    /** Constructor for a MinlpProblem and Snopt-parameters.
        @param p_ A Pointer to the problem to solve.
        @param param_ A Pointer to the parameters for the solver, default is NULL.
        @param out_solver_p_ A pointer to an ostream to print to, default is out_out_p.
        @param out_solver_log_p_ A pointer to an ostream to print logging information to, default is out_log_p.
    */
    SnOpt(const Pointer<MinlpProblem> p_, Pointer<Param> param_=NULL, char* param_prefix=NULL, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);

    /** Destructor.
        Closes the files.
    */
    ~SnOpt();

    /** Initialization part two.
        Set's start to "Cold" and updates the bounds.
    */
		void reinit();

    /** Solves the problem.
        Calls SNOPT.
        @return The return code from SNOPT: inform.
    */
    int solve();

    /** Solves the problem for a starting point.
        Set's xs and hs, so that x is taken as starting point.
				Fixes the discrete variables to the values in x.
        Set's start to "Warm".
        Call's solve().
        @return The return code from SNOPT: inform.
    */
    int solve(dvector &x);

    /** Returns the Lagrangian multipliers.
        @return A dvector with the dimension equals the number of constraints of the original problem (without the multiplier for the linear part of the objective).
    */
    dvector get_lag_multipliers();

    dvector get_mu_q();
};

#endif // SNOPT_H
#endif // SNOPT_AVAILABLE
