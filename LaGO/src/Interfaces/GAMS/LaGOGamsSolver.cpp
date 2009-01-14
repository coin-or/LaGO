// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOAlgorithm.hpp 135 2007-10-09 19:37:41Z stefan $

#include "LaGOGamsSolver.hpp"

#include <sys/wait.h>
#include <fstream>

#define CNAMES

#include "iolib.h"
//#include "dict.h"
#include "nliolib.h"
//#include "gcprocs.h"
//#include "g2dexports.h"
//#include "clicelib.h"

namespace LaGO {

pid_t childpid;
void killChild(int) {
	cerr << "Killing child because time is up." << endl;
	kill(childpid, SIGTERM);
	signal(SIGALRM, SIG_DFL);
	alarm(0);
}

int startProcess(const char* name, char*const* const args, int timelimit) {
	int status;

	childpid=vfork();
	if (!childpid) execvp(name, args); // call process in child process
	if (childpid==-1) return -1;
	
	if (timelimit) {
		alarm(timelimit); // set alarm
		signal(SIGALRM, killChild);
	}

	if (waitpid(childpid, &status, 0)==-1) { // waitpid was interrupted by SIGALRM probably
		cerr << "Waiting for process interrupted." << endl;
		alarm(0);
		return 3;
	}
	alarm(0);
	if (!WIFEXITED(status)) return -1;  // some trouble in child process
	return WEXITSTATUS(status);
}


GamsSolver::GamsSolver(const GamsReader& gams_, const vector<int>& discrete_var_)
: gams(gams_), discrete_var(discrete_var_), discr_fix(discrete_var.size()), first_call(true),
  solution_primal(iolib.ncols), solution_conval(iolib.nrows), solution_dualcon(iolib.nrows), solution_dualvar(iolib.ncols)
{
	tmpsolfn=strdup(iolib.flnsol);
	iolibsave=new tiolib;
	*(tiolib*)iolibsave=iolib;
	solvername=strdup("conopt");

	int slen=strlen(solvername);
	for (int i=0; i<iolib.nosolvers; ++i)
		//TODO: should be case insensitive
		if (strncmp(solvername, iolib.line1[i], slen)==0 && isspace(iolib.line1[i][slen])) {
			subsolver=i;
			break;
		}
	if (subsolver<0) {
		cerr << "GamsSolver: Could not find GAMS NLP solver " << solvername << endl;
		exit(EXIT_FAILURE);
	}

	args=new char*[3];
	args[0]=strdup(iolib.line3[subsolver]);
	args[1]=new char[1024];
	sprintf(args[1], "%soqcntr.scr", iolib.gscrdr);
	args[2]=NULL;
	
	optfile=0;

	// Create SBB-control-file
	sprintf(iolib.flnsbbopt, "%soqinfo.scr", iolib.gscrdr);
	ofstream file(iolib.flnsbbopt);
	file << "restart 1" << endl << "rfile oqfl000.scr" << endl;
	file.close();
}

GamsSolver::~GamsSolver() {
	if (tmpsolfn) free(tmpsolfn);
	if (solvername) free(solvername);
	if (args) {
		free(args[0]);
		delete[] args[1];
		delete[] args;
	}
	delete (tiolib*)iolibsave;
}

void GamsSolver::setStartPoint(const DenseVector& x) {
	startpoint=x;
}

void GamsSolver::setBounds(const DenseVector& lower_, const DenseVector& upper_) {
	lower=lower_;
	upper=upper_;	
}


bool GamsSolver::callSolver() {
	sprintf(iolib.flnsbbsav, "oqfl000.scr");
	sprintf(iolib.flnsbbopt, "%soqinfo.scr", iolib.gscrdr);

	// round and fix discrete variables in startpoint
	for (int i=0; i<(int)discrete_var.size(); ++i) {
		int i0=discrete_var[i];
		double val=closestInteger(startpoint[i0]);
		startpoint[i0]=val;
		lower[i0]=val;
		upper[i0]=val;
		discr_fix[i]=val;
	}

	// create start file, fake dual variables and variable and equation status
	double* duals=CoinCopyOfArray((double*)NULL, iolib.nrows);
	int* varequstatus=CoinCopyOfArray((int*)NULL, CoinMax(iolib.ncols, iolib.nrows), 3);
  if (cioSBBSave(1,1,1,1,discr_fix.getElements(), discr_fix.getElements(), startpoint.getElements(), duals, varequstatus, varequstatus)) {
		cerr << "GamsSolver: Could not write to S/R File." << endl;
		exit(EXIT_FAILURE);
	}

	iolib.SBBFlag=first_call ? 1 : 2;
	iolib.ignbas=1;

	sprintf(iolib.flnsta, "%soqsta.scr", iolib.gscrdr);
	sprintf(iolib.flnsol, "%soqsol.scr", iolib.gscrdr);

	if (iolib.ilog!=0 && iolib.ilog!=1 && iolib.ilog!=3) {
		sprintf(iolib.flnlog, "%soqlog.scr", iolib.gscrdr);
		iolib.ilog=2;
	}

	iolib.useopt=optfile;
	sprintf(iolib.flnopt, "%s%s.opt", iolib.gwrkdr, solvername);

	gfWriteCntr(args[1], &iolib, 30); // write control file
	first_call=false;
	
	// calling sub-solver, timelimit: 10 minutes
	int ret=startProcess(iolib.line3[subsolver], args, 600);
	if (ret) {
		cerr << "Spawn of GAMS NLP solver " << solvername << " failed. exit code: " << ret << endl;
		return false;
	}

  if (iolib.ilog!=0 && iolib.ilog!=1 && iolib.ilog!=3) {
		ifstream log(iolib.flnlog);
		if (!log.good()) {
			cerr << "Could not open logfile " << iolib.flnlog << " for reading." << endl;
			exit(EXIT_FAILURE);
		}
		char buf[1024];
		while (!log.eof()) {
			//TODO: shouldn't we be more careful and read at most 1024 char.?
			log >> buf;
			fprintf(gfiolog, buf);
		}
		if (remove(iolib.flnlog))
			cerr << "Could not remove file " << iolib.flnlog << endl;
	}
	if (remove(iolib.flnsta)) {
		cerr << "Could not remove file " << iolib.flnsta << endl;
		cerr << "Spawn of GAMS NLP solver " << solvername << " failed." << endl;
		return false;
	}

	int dom_violations; // domain violations
	gfrsol(iolib.flnsol, &model_status, &solver_status, &iteration_number, &time, &solution_value, &dom_violations,
		solution_conval.getElements(), solution_dualcon.getElements(), NULL, varequstatus,
		gams.con_types, gams.con_rhs,
		solution_primal.getElements(), solution_dualvar.getElements(), NULL, varequstatus,
		lower.getElements(), upper.getElements(), iolib.nrows, iolib.ncols);
	if (gams.reformed) solution_primal[gams.objvar]=0.;

	strcpy(iolib.flnsol, tmpsolfn);
	iolib=*(tiolib*)iolibsave;

	delete[] duals;
	delete[] varequstatus;
	
	return true;
}

} // namespace LaGO
