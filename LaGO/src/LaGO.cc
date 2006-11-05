// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "standard.h"
#include "param.h"
#include "problem.h"
#include "opt.h"
#include "minlpopt.h"
#include "ampl.h"
#include "gams.h"
#include "snopt.h"

int main(int argc, char** argv) {
	Timer t;
//	out_log_p=NULL;
//	out_out_p=NULL;
//	out_err_p=NULL;

	out_out << "LaGO very-early-version" << endl;
	if (argc<2) {
#ifdef AMPL_AVAILABLE		
		out_out << "usage: " << argv[0] << " stubfile (without .nl) [<parameter-file>]" << endl;
#else
#ifdef GAMS_AVAILABLE
		out_out << "usage: " << argv[0] << " <gams-file> [<parameter-file>]" << endl;
#else
		out_out << "No AMPL or GAMS interface available." << endl;
		exit(-1);
#endif
#endif	
	}

	Pointer<Param> param;
	if (argc>2) {	
		out_out << "Reading parameter file " << argv[2] << endl;
		param=new Param(argv[2], NULL);
		param->read();
	} else {
		param=new Param(NULL, NULL);
	}
	
// 	cout << "Reading problem-file " << argv[1] << endl;
#ifdef AMPL_AVAILABLE
	ampl interface(argv[1]);
	Pointer<MinlpProblem> prob(interface.get_problem());
#else
	gams interface(param);
	Pointer<MinlpProblem> prob(interface.get_problem(argv[1]));
	prob->prob_name=strdup(argv[1]);
#endif

#ifdef SNOPT_AVAILABLE
	snoptlicenceok=true;
#endif	
	
	if (!LocOpt::nlp_solver_available()) cout << "No ";
	cout << "NLP Solver available." << endl;

	cout << "Constructing Solver" << endl;
	Pointer<MinlpOpt> solver(new MinlpOpt(prob, param, true));
//	Pointer<LocOpt> solver=LocOpt::get_solver_origprob(prob, param);
	cout << "Solving ..." << endl;
	int ret=solver->solve();

	cout << "Solved: " << ret << endl;
	int ns=prob->feasible(solver->sol_point, 1.E-4, ret ? NULL : out_log_p);
	cout << "Nonsatisfied constraints: " << ns << endl;
	cout.setf(ios::fixed);
	cout.precision(20);
	cout << "Optimal value: " << prob->obj->eval(solver->sol_point) << endl;
	cout.unsetf(ios::fixed);
	cout.precision(6);

	cout << "Time: " << t.stop() << endl;

#ifdef AMPL_AVAILABLE
  interface.write_sol_file(solver->sol_point);
#else
	int model_status=ret ? 11 : 1;
	int solver_status=1;
	interface.write_sol_file(solver->sol_point, model_status, solver_status, solver->iter(), t, prob);
#endif

	cout << "LaGO finished." << endl;
	return 0;
}
