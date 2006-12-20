// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef GAMS_H
#define GAMS_H

#include "standard.h"
#ifdef COIN_HAS_GAMSIO
#define GDX_AVAILABLE

#include "param.h"
#include "problem.h"
#include "opt.h"
#include "preprocessapi.h"

class gamsLocOpt;
class SolCandidate;

extern "C" struct dictRec;

/** Gams interface.
		@class gams
		@param GAMS write solution candidates
		%options $\geq 0$
		%default 0
		%level 2
		If 0, no solution candidates are written to gdx files.
		If greater 0, every solution candidate is written to a file solcand$<$Value$>$\_$<$Randomnumber$>$.gdx, where Value is the objective value.
		If more files were written than the number, which is specified by this parameter, the file with the worst solution candidate is removed again.
		@param GAMS write startpoint
		%options $\geq 0$
		%default 0
		%level 2
		If greater 0, every starting point of a local search is written to a file startpoint\_$<$Randomnumber$>$.gdx.
*/
class gams {
	friend class gamsLocOpt;
	friend class gamsFunc;

	private:
		Pointer<Param> param;

		// original types of constraints
		ivector con_type;

		dvector rhs;
		dvector lower, upper;
		double obj_sign;
		int objcon, objvar;
		Pointer<char> objcon_name;

		bool reformed;

		/** The list of the written gdx files.
		*/
		multimap<double, Pointer<char> > written_gdx;
		int written_gdx_limit;

		/** get name of row i
		    @param dict
		    @param gi row index, [0..nRows)
 		    @param bufLen size of target buffer
				@param name target buffer for row name
        @return target on success, NULL on failure
		*/
		char* getRowName (int i, char *name, int bufLen);
		
		/** get name of column j
		    @param dict
		    @param gj column index, [0..nCols)
 		    @param bufLen size of target buffer
				@param name target buffer for column name
        @return target on success, NULL on failure
		*/
		char* getColName (int j, char *name, int bufLen);
#ifdef GDX_AVAILABLE
		void gdx_error(int n);
#endif
	public:
		static void init_cplex_licence(int connr=0, int varnr=0, int nnz=0, int nlnz=0, int ndisc=0);
		
		void init_snopt_licence();

		struct dictRec* dict;

		gams(Pointer<Param> param_=NULL);

		~gams();

		Pointer<MinlpProblem> get_problem(char* gamsfile);

		void write_sol_file(const dvector& sol_point, int model_status, int solver_status, int iter, double time, Pointer<MinlpProblem> prob);
		void write_sol_set(const set<SolCandidate>& sol_set);

		void write_matlab(const dvector& x, const char* filename, vector<Pointer<char> >& var_names);

		void write_gams(const dvector& x, const char* filename, const vector<bool>& discr);
#ifdef GDX_AVAILABLE
		/** Writes a gdx file which contains the variables and the values from a given vector.
		*/
		void write_gdx(const dvector& x, char* filename, double val);
		void write_gdx(const dvector& x, char* filename);
		void read_gdx(dvector& x, char* filename);
#endif

		void write_box(const dvector& lower, const dvector& upper);
};

// points to gams-interface; needed for gamsLocOpt
extern gams* gamsptr;

class gamsNLData {
public: 
	/** Stripped NL Instruction data. */
	unsigned int* instr;
	double* nlCons;
	int* startIdx;
	int* numInstr;
	
	/** Sum of Instruction lengths. */
	int lenins;
	/** Max instruction length. */
	int maxins;
	
	/** Used for function evaluation. */
	double* s;
	double* sbar;
	double* resstack;
	int resstacksize;
		
	int* jacNX;
	int* jacVR;
	double* jacVL;
	int* hesLagCL;
	int* hesLagRW;
	int* hesLagNX;
	double* hesLagVL;
	
	int jacOneLenMax;
	int hesLagHeadPtr;
	int hesOneLenMax;
	int hesLagLenMax;
	
	/** Number of domain violations so far.
	*/
	int domain_violations;
	
	gamsNLData()
	: instr(NULL), nlCons(NULL), startIdx(NULL), numInstr(NULL),
		lenins(1), maxins(1), resstacksize(1),
		s(NULL), sbar(NULL),
		jacNX(NULL), jacVR(NULL), jacVL(NULL),
		hesLagCL(NULL), hesLagRW(NULL), hesLagNX(NULL), hesLagVL(NULL),
		jacOneLenMax(-1),
		hesLagHeadPtr(-1), hesOneLenMax(-1), hesLagLenMax(-1),
		domain_violations(0)
	{ }
	
	~gamsNLData() {
		if (instr) delete[] instr;
		if (nlCons) delete[] nlCons;
		if (startIdx) delete[] startIdx;
		if (numInstr) delete[] numInstr;
		if (s) delete[] s;
		if (sbar) delete[] sbar;
		if (resstack) delete[] resstack;
		if (jacNX) delete jacNX;
		if (jacVR) delete jacVR;
		if (jacVL) delete jacVL;
		if (hesLagCL) delete hesLagCL;
		if (hesLagRW) delete hesLagRW;
		if (hesLagNX) delete hesLagNX;
		if (hesLagVL) delete hesLagVL;
	}

};
	
class gamsFunc : public Func {
	private:
		Pointer<gamsNLData> data;
	
		int connr;

		Func::CurvatureType curv_type;

	public:
		gamsFunc(int n, int connr_, Pointer<gamsNLData> data_, Pointer<SparsityInfo> sparsity_=NULL)
		: Func(n), connr(connr_), data(data_), curv_type(Func::UNKNOWN)
		{ sparsity=sparsity_; }
		
		double eval(const UserVector<double>& x) const {
			return eval((Pointer<double>)x);
		}
		
		double eval(const double* x) const;
		
		int valgrad(double& val, UserVector<double>& g, const UserVector<double>& x) const {
			dvector g0(g.dim());
			int ret=valgrad(val, (Pointer<double>)g0, (Pointer<double>)x);
			g=g0;
			return ret;
		}

		int valgrad(double& val, double* g, const double* x) const;

		void grad(UserVector<double>& g, const UserVector<double>& x) const {
			dvector g0(g.dim());
			double val;
			valgrad(val, (Pointer<double>)g0, (Pointer<double>)x);
			g=g0;
		}

		void grad(dvector& g, const dvector& x) const {
			double val;
			valgrad(val, (Pointer<double>)g, (Pointer<double>)x);
		}

		void grad(double* g, const double* x) const {
			double val;
			valgrad(val, g, x);
		}
		using Func::grad;

		void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
			dvector y0(y.dim());
			HessMult((Pointer<double>)y0, (Pointer<double>)x, (Pointer<double>)z);
			y=y0;
		}
		
		void HessMult(dvector& y, const dvector& x, const dvector& z) const {
			HessMult((Pointer<double>)y, (Pointer<double>)x, (Pointer<double>)z);
		}
		
		void HessMult(double* y, const double* x, const double* z) const;
		using Func::HessMult;

#ifdef FILIB_AVAILABLE
		bool is_interval_compliant() const { return true; } // :-)

		interval<double> eval(const IntervalVector& x) const;

		int valgrad(interval<double>& val, IntervalVector& y, const IntervalVector& x) const;
#endif

		void set_curvature(CurvatureType ct) { curv_type=ct; };
		CurvatureType get_curvature() const { return curv_type; };

		void print(ostream& out) const {
			out << "GamsFunc " << connr << endl;
		}
};

/** Calls a GAMS-solver to solve the original problem with different bounds, fixed variables or another starting point.
		@class gamsLocOpt
    @param GAMS LocOpt solver
    %options name of GAMS solver
    %default conopt
		%level 2
    Name of the GAMS solver to use to solve the original MINLP with fixed binaries.
		@param GAMS LocOpt optionfile
		%options $\geq 0$
		%default 0
		%level 2
		The number of the optionalfile for the local optimizer. 0 for no optionfile.
		@param GAMS LocOpt preprocessing
		%options name of a script
		%default none
		%level 2
		The library to link which defines some algorithm to do some preprocessing for a given startpoint before giving it to the local optimizer.
		@param GAMS write solution candidates
		See gams-section.
*/
class gamsLocOpt : public LocOpt {
	friend class gams;
	private:
		Pointer<Param> param;
		Pointer<MinlpProblem> prob;
		bool second_run;

		dvector lower_discr;
		dvector upper_discr;
		dvector lower;
		dvector upper;

		dvector con_val;
		dvector duals_con;
		dvector duals_var;
		ivector basind_con;
		ivector basind_var;

		int model_status, solver_status;

		int recent_calls;
		Pointer<char> tmpsolfn; // temporary space to store solutionfilename
		char** args; // arguments for subsolver-call

		int optfile;
		int subsolver;
		Pointer<char> solvername;
		void* iolibsave;

		/** Skript name for a preprocessing step.
		*/
		char** preprocessargs;
		bool preprocess_keepgdx;

		bool write_solcand;
		bool write_startpoint;

		/** A random value to identify points, we process.
		*/
		double rd;

		/** Calls a preprocessing program.
		    The improved point will be written to sol_point.
		*/
		int do_preprocessing(const dvector& start);

		static void* preprocess_handle;

		LocOptPreprocessing* preprocessing;

	public:
		gamsLocOpt(Pointer<MinlpProblem> prob_, Pointer<Param> param_, Pointer<ostream> out_solver_p_=out_out_p, Pointer<ostream> out_solver_log_p_=out_log_p);
	
	  ~gamsLocOpt();
	
		int solve() { return solve(prob->primal_point); }
	
		int solve(dvector& start);
	
		dvector get_lag_multipliers();
};

#endif // COIN_HAS_GAMSIO
#endif
