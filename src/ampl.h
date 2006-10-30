// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske


#ifdef AMPL_AVAILABLE
#ifndef AMPL_H
#define AMPL_H

#include "standard.h"
#include "func.h"
#include "problem.h"

#include "asl.h"
#undef filename

/** Global pointer to Ampl-Interface.
    Maybe needed by the ASL-functions.
*/
extern ASL *asl;

/** Class for the AMPL-Interface.
*/
class ampl {
private:
  /** The name of the AMPL .nl-file to read.
  */
  Pointer<char> stubfile;

public:

  /** Constructor for the name of the .nl file.
      If stubfile_ if NULL, gives an error-message and exists.
      @param stubfile_ The name of the .nl-file without ".nl"
  */
  ampl(char* stubfile_) {
    if (! stubfile_) {
      out_err << "ampl::ampl: Given filename was NULL !" << endl;
      exit(-1);
    }
    stubfile = strdup(stubfile_);
  }

  /** Read's the AMPL file and constructs a MINLP-problem.
      @return The problem, which was described in the AMPL-file.
  */
  Pointer<MinlpProblem> get_problem();

  /** Write the AMPL solution-file stub.sol.
      @param sol_point The solution point.
      @param message A message.
  */
  void write_sol_file(const UserVector<double>& sol_point, char* message="") {
    write_sol(message, (double*)(const Pointer<double>)sol_point, NULL, NULL);
  }

};

/** Class for an AMPL-objective-function.
    This class is used to compute a nonlinear objective of a problem, which was described by an ampl-file.
*/
class amplObj : public Func {
private:
	Func::CurvatureType curv_type;
public:
  /** (Standard-)Constructor with optional the dimension.
      @param dim_ The dimension of the function.
      @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
      @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
  */
  amplObj(int dim_=0, Pointer<SparsityInfo> sparsity_=NULL, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
  : Func(dim_, out_func_p_, out_func_log_p_)
  { sparsity=sparsity_;
  }

  /** Evaluates the objective for a UserVector<double>.
      Calls objval(...) from the AMPL/Solver Interface Library (ASL).
      @param x The UserVector<double> to evaluate.
      @return The value of the objective.
  */
  double eval(const UserVector<double>& x) const {
  	delete asl->i.err_jmp_;
  	asl->i.err_jmp_=new Jmp_buf;
	  if (setjmp(asl->i.err_jmp_->jb)) {
  		out_log << "ASL-error evaluating objective" << endl;
  		out_log << "Last point: " << endl;
  		for (int i=0; i<x.dim(); i++) out_log << var_name(i) << ": " << x(i) << endl;
	    write_sol("ASL-error evaluating objective, aborted", (double*)(const Pointer<double>)x, NULL, NULL);
  		exit(0);
	  }
    return objval(0, (double*)(const Pointer<double>)x, NULL);
  }

  /** Computes the gradient of the objective for a dvector.
      Calls objgrd(...) from the AMPL/Solver Interface Library (ASL).
      @param g The dvector to store the result in.
      @param x The dvector to compute the gradient for.
  */
	void grad(dvector& g, const dvector& x) const {
  	delete asl->i.err_jmp_;
  	asl->i.err_jmp_=new Jmp_buf;
	  if (setjmp(asl->i.err_jmp_->jb)) {
  		out_log << "ASL-error evaluating gradient of objective" << endl;
  		out_log << "Last point: " << endl;
  		for (int i=0; i<x.dim(); i++) out_log << var_name(i) << ": " << x(i) << endl;
	    write_sol("ASL-error computing gradient of objective, aborted", (double*)(const Pointer<double>)x, NULL, NULL);
  		exit(0);
	  }
    objgrd(0, (double*)(const Pointer<double>)x, (double*)(Pointer<double>)g, NULL);
	}

  void grad(UserVector<double>& g, const UserVector<double>& x) const {
  	delete asl->i.err_jmp_;
  	asl->i.err_jmp_=new Jmp_buf;
	  if (setjmp(asl->i.err_jmp_->jb)) {
  		out_log << "ASL-error computing gradient of objective" << endl;
  		out_log << "Last point: " << endl;
  		for (int i=0; i<x.dim(); i++) out_log << var_name(i) << ": " << x(i) << endl;
	    write_sol("ASL-error computing gradient of objective, aborted", (double*)(const Pointer<double>)x, NULL, NULL);
  		exit(0);
	  }
    double* gr=new double[g.size()];
    objgrd(0, (double*)(const Pointer<double>)x, gr, NULL);
    g.set(gr, g.dim());
    delete gr;
  }
	
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
	using Func::grad;
#endif

  /** Computes the product of a Hessian and a dvector.
      @param y The dvector to store the result in.
      @param x The dvector to compute the Hessian for.
      @param z The dvector to multiply with the Hessian.
  */
  void HessMult(dvector& y, const dvector& x, const dvector& z) const {
    eval(x);
    hvcomp((double*)(const Pointer<double>)y, (double*)(const Pointer<double>)z, 0, NULL, NULL);
  }

  void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
    eval(x);
    double* y0=new double[y.size()];
    hvcomp(y0, (double*)(const Pointer<double>)z, 0, NULL, NULL);
    y.set(y0, y.dim());
    delete y0;
  }
	
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
	using Func::HessMult;
#endif

	void set_curvature(CurvatureType ct) { curv_type=ct; };
	CurvatureType get_curvature() const { return curv_type; };

  /** Prints some information about this function.
      Prints the dimension.
      @param out The ostream to print to.
  */
  void print(ostream& out) const {
    out << "amplObj: dim=" << dim() << endl;
  }

};

/** Class for an AMPL-constraints-function.
*/
class amplCon : public Func {
private:
  /** The number of the constraint, this function represents.
  */
  fint connr;

  /** Array for Hessain-Multiplication, to select the constraint, which Hessian should be computed.
  */
  Pointer<real> con_select;

	Func::CurvatureType curv_type;
public:
  /** Constructor for the number of the constraint and the dimension of the function.
      @param connr_ The number of the constraint, this function represents.
      @param dim_ The dimension.
      @param conv_ Indicates, whether the function is convex, default is false.
      @param out_func_p_ An ostream to print function-related output to, default is out_out_p.
      @param out_func_log_p_ An ostream to print function-related logging-output to, default is out_func_p.
      @see amplCon(int)
  */
  amplCon(int connr_, int dim_=0, Pointer<SparsityInfo> sparsity_=NULL, Pointer<ostream> out_func_p_=out_out_p, Pointer<ostream> out_func_log_p_=out_log_p)
  : Func(dim_, out_func_p_, out_func_log_p_), connr(connr_), con_select(new real[n_con])
  { for (int i=0; i<n_con; i++) con_select[i]=0;
    con_select[connr]=1;
    sparsity=sparsity_;
  }

  /** Evaluates a constraint for a UserVector<double>.
      Calls conival(connr, ...) from the AMPL/Solver Interface Library (ASL).
      @param x The UserVector<double> to evaluate.
      @return The value of the connr-th constraint.
  */
  double eval(const UserVector<double>& x) const {
  	delete asl->i.err_jmp_;
  	asl->i.err_jmp_=new Jmp_buf;
	  if (setjmp(asl->i.err_jmp_->jb)) {
  		out_log << "ASL-error evaluating constraint " << connr << ": " << con_name(connr) << endl;
  		out_log << "Last point: " << endl;
  		for (int i=0; i<x.dim(); i++) out_log << var_name(i) << ": " << x(i) << endl;
  		char* str=new char[200]; sprintf(str, "ASL-error evaluating constraint %i: %s", connr, con_name(connr));
	    write_sol(str, (double*)(const Pointer<double>)x, NULL, NULL);
  		exit(0);
  }
  	return conival(connr, (double*)(const Pointer<double>)x, NULL);
  }

  /** Computes the gradient of a constraint for a dvector.
      Calls congrd(connr, ...) from the AMPL/Solver Interface Library (ASL).
      @param g The dvector to store the connr-th gradient in.
      @param x The dvector to compute the gradient for.
  */
  void grad(dvector &g, const dvector &x) const {
  	delete asl->i.err_jmp_;
  	asl->i.err_jmp_=new Jmp_buf;
	  if (setjmp(asl->i.err_jmp_->jb)) {
  		out_log << "ASL-error computing gradient of constraint " << connr << ": " << con_name(connr) << endl;
  		out_log << "Last point: " << endl;
  		for (int i=0; i<x.dim(); i++) out_log << var_name(i) << ": " << x(i) << endl;
  		char* str=new char[200]; sprintf(str, "ASL-error computing gradient of constraint %i: %s", connr, con_name(connr));
	    write_sol(str, (double*)(const Pointer<double>)x, NULL, NULL);
  		exit(0);
	  }
    congrd(connr, (double*)(const Pointer<double>)x, (double*)(Pointer<double>)g, NULL);
  }

  void grad(UserVector<double>& g, const UserVector<double>& x) const {
  	delete asl->i.err_jmp_;
  	asl->i.err_jmp_=new Jmp_buf;
	  if (setjmp(asl->i.err_jmp_->jb)) {
  		out_log << "ASL-error computing gradient of constraint " << connr << ": " << con_name(connr) << endl;
  		out_log << "Last point: " << endl;
  		for (int i=0; i<x.dim(); i++) out_log << var_name(i) << ": " << x(i) << endl;
  		char* str=new char[200]; sprintf(str, "ASL-error computing gradient of constraint %i: %s", connr, con_name(connr));
	    write_sol(str, (double*)(const Pointer<double>)x, NULL, NULL);
  		exit(0);
	  }
    double* gr=new double[g.size()];
    congrd(connr, (double*)(const Pointer<double>)x, gr, NULL);
    g.set(gr, g.dim());
    delete gr;
  }
	
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
	using Func::grad;
#endif

  /** Computes the product of a Hessian and a dvector.
      @param y The dvector to store the result in.
      @param x The dvector to compute the Hessian for.
      @param z The dvector to multiply with the Hessian.
  */
  void HessMult(dvector& y, const dvector& x, const dvector& z) const {
    eval(x);
    hvcomp((double*)(Pointer<double>)y, (double*)(const Pointer<double>)z, -1, NULL, con_select);
  }

  void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
    eval(x);
    double* y0=new double[y.size()];
    hvcomp(y0, (double*)(const Pointer<double>)z, -1, NULL, con_select);
    y.set(y0, y.dim());
    delete y0;
  }
	
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
	using Func::HessMult;
#endif

	void set_curvature(CurvatureType ct) { curv_type=ct; };
	CurvatureType get_curvature() const { return curv_type; };
	
  /** Prints some information about this function.
      Prints the dimension and constraint-number.
      @param out The ostream to print to.
  */
  void print(ostream& out) const {
    out << "amplCon " << connr << ": dim=" << dim() << ": " << con_name(connr) << endl;
  }
};


#endif // AMPL_H
#endif // AMPL_AVAILABLE
