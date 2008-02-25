// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske


#include "tools.h"
#include <unistd.h>
#include <stdlib.h>

extern "C" {
#include "ranlib.h"
}

#include <sys/wait.h>
//#include <signal.h>

map<const void*, int> refcount;

void print_all_Pointer(ostream& out) {
  for (map<const void*, int>::iterator p=refcount.begin(); p!=refcount.end(); p++) {
    out << p->first << ": " << p->second << endl;
  }
}

// --------------------------------------------------------

Pointer<ostream> out_out_p(&cout, false);
Pointer<ostream> out_log_p(&clog, false);
Pointer<ostream> out_err_p(&cerr, false);

// --------------------------------------------------------

// initalize random-number-generator
//void* initstatereturn=initstate(time(NULL), new char[32], 32);
/*
int init_random() {
  srand(time(NULL));
  return 0;
}
int init_random_dummy=init_random();
*/

/* uniformly distributed random number from (0,1)
   * using the C system random number generator
*/
/*inline double uniform() {
   return (double)rand() / (RAND_MAX - 1);
}

int random(int lb, int ub) {
   return lb+(int)((ub-lb+0.99) * uniform()); // as suggested in the rand()-manual
}

double random(double lb, double ub) {
   return lb + (ub - lb) * uniform();
}
*/

int random(int lb, int ub) {
	return (int)random((double)lb, (ub+0.99));
}

double random(double lb, double ub) {
	if (lb<=-INFINITY && ub>=INFINITY) return gennor(0., 1.);
	if (lb<=-INFINITY) return ub-genexp(1.);
	if (ub>= INFINITY) return lb+genexp(1.);
	return genunf(lb, ub);
}

// --------------------------------------------------------

void Timer::start() {
//	starttime=clock();
	struct rusage r;

	getrusage(RUSAGE_SELF, &r);
//	starttime=r.ru_stime;
	starttime=r.ru_utime;
//	starttime.tv_sec+=r.ru_stime.tv_sec;
//	starttime.tv_usec+=r.ru_stime.tv_usec;

	getrusage(RUSAGE_CHILDREN, &r);
	starttime.tv_sec+=r.ru_utime.tv_sec;
	starttime.tv_usec+=r.ru_utime.tv_usec;
//	starttime.tv_sec+=r.ru_stime.tv_sec;
//	starttime.tv_usec+=r.ru_stime.tv_usec;
}

Timer& Timer::stop() {
//	endtime=clock();
	struct rusage r;

	getrusage(RUSAGE_SELF, &r);
//	endtime=r.ru_stime;
	endtime=r.ru_utime;
//	endtime.tv_sec+=r.ru_stime.tv_sec;
//	endtime.tv_usec+=r.ru_stime.tv_usec;

	getrusage(RUSAGE_CHILDREN, &r);
	endtime.tv_sec+=r.ru_utime.tv_sec;
	endtime.tv_usec+=r.ru_utime.tv_usec;
//	endtime.tv_sec+=r.ru_stime.tv_sec;
//	endtime.tv_usec+=r.ru_stime.tv_usec;

	return *this;
}

Timer::operator double() const {
/*	if (endtime>starttime) return (endtime-starttime)/(double)CLOCKS_PER_SEC;
	return (endtime+(((clock_t)-1)-starttime))/(double)CLOCKS_PER_SEC;
*/	return (double)endtime.tv_sec+(double)endtime.tv_usec/1.E+6-(double)starttime.tv_sec-(double)starttime.tv_usec/1.E+6;
}


void Timer::print(ostream& out, double time) {
	long sec=(long)time;
	long usec=(long)((time-(double)(long)time)*1000000);

  usec/=10000;
  int hours = sec/3600;
  int mins  = (sec%3600)/60;
  int secs  = sec%60;
  if (hours)     out << hours << ":";
  if (mins > 9)  out << mins  << ":";
  else if (mins && hours)
                 out << "0" << mins << ":";
  else if (mins) out << mins << ":";
  else if (hours) out << "00:";   // hours and not minutes.
  if (secs > 9)  out << secs;
  else if (secs && (hours || mins))
                 out << "0" << secs;
  else if (secs) out << secs;    // nor hours nor minutes
  else if (hours || mins) out << "00"; // not secs, but higher.
  else out << "0";  // nor secs, hours, nor mins.
  if (usec > 9)  out << "." << usec;
  else if (usec) out << ".0" << usec;
  else           out << ".00";

}

#include <malloc.h>

unsigned int get_mem() {
//struct mallinfo {
//  int arena;    /* non-mmapped space allocated from system */
//  int ordblks;  /* number of free chunks */
//  int smblks;   /* number of fastbin blocks */
//  int hblks;    /* number of mmapped regions */
//  int hblkhd;   /* space in mmapped regions */
//  int usmblks;  /* maximum total allocated space */
//  int fsmblks;  /* space available in freed fastbin blocks */
//  int uordblks; /* total allocated space */
//  int fordblks; /* total free space */
//  int keepcost; /* top-most, releasable (via malloc_trim) space */
//};

	struct mallinfo meminfo=mallinfo();
//	out_log << meminfo.arena << '\t' << meminfo.ordblks << '\t' << meminfo.smblks << '\t'
//		<< meminfo.hblks << '\t' << meminfo.hblkhd << '\t' << meminfo.usmblks << '\t' << meminfo.fsmblks << '\t'
//		<< meminfo.uordblks << '\t' << meminfo.fordblks << '\t' << meminfo.keepcost << endl;

	return meminfo.arena+meminfo.uordblks;
}

pid_t childpid;
void kill_child(int) {
	out_err << "Killing child because time is up." << endl;
	kill(childpid, SIGTERM);
	signal(SIGALRM, SIG_DFL);
	alarm(0);
}

int start_process(const char* name, char*const* const args, int timelimit, const char* envvarname, const char* envvarvalue) {
	int status;

	childpid=vfork();
	if (!childpid) execvp(name, args); // call process in child process
	if (childpid==-1) return -1;
	
	if (envvarname) {
		assert(envvarvalue);
		setenv(envvarname, envvarvalue, 1);	
	}

	if (timelimit) {
		alarm(timelimit); // set alarm
		signal(SIGALRM, kill_child);
	}

	if (waitpid(childpid, &status, 0)==-1) { // waitpid was interrupted by SIGALRM probably
		out_err << "Waiting for process interrupted." << endl;
		alarm(0);
		return 3;
	}
	alarm(0);
	if (!WIFEXITED(status)) return -1;  // some trouble in child process
	return WEXITSTATUS(status);
}

