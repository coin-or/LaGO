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

int main(int argc, char** argv) {
	Timer t;
//	out_log_p=NULL;
//	out_out_p=NULL;
//	out_err_p=NULL;

	out_out << "LaGO very-early-version" << endl;
	if (argc<2) {
		out_out << "usage: " << argv[0] << " stubfile (without .nl) [-AMPL] [parameter-file in resource/-subdir]" << endl;
		exit(-1);
	}
	
	out_out << "Reading parameters." << endl;
	Pointer<Param> param=new Param(argc>2 ? argv[2] : NULL, "");
	param->read();
	
	cout << "Reading problem-file " << argv[1] << endl;
	Pointer<MinlpProblem> prob;
	ampl interf(argv[1]);
	prob=interf.get_problem();
//	prob->prob_name=strdup(argv[1]);
	
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

    interf.write_sol_file(solver->sol_point);
	
	cout << "LaGO finished." << endl;
	return 0;
}
