#include "preprocessapi.h"

#ifdef _WIN32
#define VIS
#else
#ifdef __sun
#define SOL
#define NOSTRNICMP
#else
#define LXY
#endif
#endif
#define CNAMES
#ifdef NOUNDERSCORE
#define FNAME_LCASE_NODECOR
#else
#define FNAME_LCASE_DECOR
#endif

#include "iolib.h"
#include "dict.h"
#include "nliolib.h"
#include "gcprocs.h"
#include "g2dexports.h"

/** Preprocessing by calling another model.
		@class GAMSPreprocessing
    Let I be the set of variables, which are shared between the model in LaGO and the model for preprocessing.
		For a given point x, we start the preprocessing, where the startingvalues of the variables in I are set to the values in x.
		If it found a feasible point, the values in x for the variables in I are set to the values from the feasible point.
		@param GAMS Preprocessing solver
    %options A GAMS solver
		%default conopt
		%level 1
		The solver to use to solve the preprocessing model.
		@param GAMS Preprocessing optionfile
		%options integer $\geq 0$
		%default 0
		%level 2
		The optionfile for the solver to use for the preprocessing model.
*/
class GAMSPreprocessing : public LocOptPreprocessing {
	private:
		tiolib preprocess_iolib, LaGO_iolib;
		struct dictRec* pp_dict;

		int preprocess_dim;

		char** preprocess_solver_args;
		int subsolver;

		/** Maps variable indices from LaGO to preprocessing and from preprocessing to LaGO.
		*/
		map<int,int> var_LaGO2pp, var_pp2LaGO;
		/** The variable indices from LaGO of variables, which do not appear in the preprocessing. And 0, projected onto the box of this variable.
		*/
		list<pair<int,double> > var_LaGO_notpp;

		int create_varmapping(struct dictRec* LaGO_dict);

		vector<int> i_discr;
		dvector fix_discr;

		dvector start, sol_point;

		bool recent_calls;

		dvector con_val, con_duals, rhs, var_duals, lower, upper;
		ivector con_type, con_basind, var_basind;

	public:
		GAMSPreprocessing(const Pointer<Param>& param_)
		: LocOptPreprocessing(param_), pp_dict(NULL), subsolver(-1), preprocess_solver_args(NULL),
		  recent_calls(false)
		{ }

		~GAMSPreprocessing() {
			iolib=preprocess_iolib;
			gfclos();
			iolib=LaGO_iolib;
			if (preprocess_solver_args) {
				delete preprocess_solver_args[0];
				delete preprocess_solver_args[1];
				delete preprocess_solver_args;
			}
		}

		int init(struct dictRec* LaGO_dict, const MinlpProblem& LaGO_prob);

		int run(UserVector<double>& x);
};

int GAMSPreprocessing::create_varmapping(struct dictRec* LaGO_dict) {
  char quote, *targ, *end, *s, tbuf[32];
  int uelIndices[10], nIndices, pp_symIndex, oldSym=-1, LaGO_symIndex, LaGO_index;

	for (int i=0; i<preprocess_dim; ++i) {
  	if (gcdColUels(pp_dict, i, &pp_symIndex, uelIndices, &nIndices) != 0) return 1;

		if (oldSym!=pp_symIndex) { // starting a new symbol
		  if ((s=gcdSymName(pp_dict, pp_symIndex, tbuf, 32)) == NULL) return 1; // get symbol name

			// get symbol index in LaGO
			LaGO_symIndex=gcdSymIndex(LaGO_dict, s);
			if (LaGO_symIndex<0) {
				out_err << "Cannot find symbol " << s << " in LaGO problem." << endl;
				continue;
			}

			oldSym=pp_symIndex;
		}

		// translate uelIndices from pp to LaGO
		for (int k = 0;  k < nIndices;  ++k) {
    	if ((gcdUelLabel(pp_dict, uelIndices[k], tbuf, 32, &quote))==NULL) return 1;
			uelIndices[k]=gcdUelIndex(LaGO_dict, tbuf);
		}

		// get variable index in LaGO
		if ((LaGO_index=gcdColIndex(LaGO_dict, LaGO_symIndex, nIndices, uelIndices, LaGO_index))<0) {
			out_err << "Cannot find variable " << s << '(';
	  	gcdColUels(pp_dict, i, &pp_symIndex, uelIndices, &nIndices);
			for (int k=0; k<nIndices; ++k) out_err << gcdUelLabel(pp_dict, uelIndices[k], tbuf, 32, &quote) << ',';
			out_err << ") in LaGO problem." << endl;
			continue;
		}

		var_pp2LaGO.insert(pair<int,int>(i, LaGO_index));
		var_LaGO2pp.insert(pair<int,int>(LaGO_index, i));

//		out_log << "Match preprocess " << i << " to LaGO " << LaGO_index << endl;
  }

	return 0;
}

int GAMSPreprocessing::init(struct dictRec* LaGO_dict, const MinlpProblem& LaGO_prob) {
	char* controlfile=getenv("GAMSPREPROCESS_CNTRFILE");
	if (!controlfile) {
		out_err << "Environmentvariable GAMSPREPROCESS_CNTRFILE not set." << endl;
		return 4;
	}

	LaGO_iolib=iolib; // save LaGO's iolib

	gfinit();

	// Read control file
	cntrec info;
	gfrcnt(true, true, &info, controlfile);
	preprocess_dim=info.kgv[2];

	start.resize(preprocess_dim);
	sol_point.resize(preprocess_dim);

	con_val.resize(info.kgv[1]);
	con_duals.resize(info.kgv[1]);
	rhs.resize(info.kgv[1]);
	con_type.resize(info.kgv[1]);
	con_basind.resize(info.kgv[1]);
	var_duals.resize(preprocess_dim);
	lower.resize(preprocess_dim);
	upper.resize(preprocess_dim);
	var_basind.resize(preprocess_dim);
	i_discr.reserve(iolib.ndisc);

	gfopst(); // open status file

	// to allow reading of columns
//	rowrec rowdata;
//	for (int c=0; c<info.kgv[1]; c++) gfrrow(&rowdata);

	// Read startvalues and discrete variables
/*	colrec coldata;
	for (int i=0; i<preprocess_dim; ++i) {
		gfrcol(&coldata);
		if (coldata.idata[3]) i_discr.push_back(i);
		start[i]=coldata.cdata[2];
	}
	lower_discr.resize(i_discr.size());
	upper_discr.resize(i_discr.size());
*/
	// load dictionary
  if (!iolib.dictFileWritten) {
		out_err << "Do not have dictionary file." << endl;
		return 5;
	}
  if (gcdLoad(&pp_dict, iolib.flndic, iolib.dictVersion)) {
		out_err << "Loading dictionary file failed" << endl;
		return 5;
	}
	if (create_varmapping(LaGO_dict)) return 6;

	int i=0;
	for (map<int,int>::iterator it(var_LaGO2pp.begin()); it!=var_LaGO2pp.end(); ++it, ++i)
		while (i<it->first) {
			if (LaGO_prob.lower(i)<=0. && LaGO_prob.upper(i)>=0.) var_LaGO_notpp.push_back(pair<int,double>(i, 0.));
			else if (LaGO_prob.lower(i)>0) var_LaGO_notpp.push_back(pair<int,double>(i, LaGO_prob.lower(i)));
			else var_LaGO_notpp.push_back(pair<int,double>(i, LaGO_prob.upper(i)));
			++i;
		}
	while (i<LaGO_prob.dim()) {
		if (LaGO_prob.lower(i)<=0. && LaGO_prob.upper(i)>=0.) var_LaGO_notpp.push_back(pair<int,double>(i, 0.));
		else if (LaGO_prob.lower(i)>0) var_LaGO_notpp.push_back(pair<int,double>(i, LaGO_prob.lower(i)));
		else var_LaGO_notpp.push_back(pair<int,double>(i, LaGO_prob.upper(i)));
		++i;
	}

	// get indices of discrete variables and copy their bounds
	for (int i=0; i<LaGO_prob.i_discr.size(); ++i) {
		map<int,int>::iterator it(var_LaGO2pp.find(LaGO_prob.i_discr[i]));
		if (it!=var_LaGO2pp.end()) i_discr.push_back(it->second);
	}
	fix_discr.resize(i_discr.size());
/*	upper_discr.resize(i_discr.size());
	for (int i=0; i<i_discr.size(); ++i) {
		lower_discr[i]=LaGO_prob.lower(var_pp2LaGO[i_discr[i]]);
		upper_discr[i]=LaGO_prob.upper(var_pp2LaGO[i_discr[i]]);
	}
*/
	// copy bounds of all variables
	for (map<int,int>::iterator it(var_LaGO2pp.begin()); it!=var_LaGO2pp.end(); ++it) {
		lower[it->second]=LaGO_prob.lower(it->first);
		upper[it->second]=LaGO_prob.upper(it->first);
	}

	// get solver index
	char* preprocess_solver=param ? param->get("GAMS Preprocessing solver", "conopt") : strdup("conopt");
	int slen=strlen(preprocess_solver);
	for (int i=0; i<iolib.nosolvers; ++i)
		if (strnicmp(preprocess_solver, iolib.line1[i], slen)==0 && isspace(iolib.line1[i][slen])) {
			subsolver=i;
			break;
		}
	if (subsolver<0) {
		out_err << "Could not find GAMS NLP solver " << preprocess_solver << endl;
		return 2;
	}

	preprocess_solver_args=new char*[3];
	preprocess_solver_args[0]=strdup(iolib.line3[subsolver]);
	preprocess_solver_args[1]=new char[1024];
	sprintf(preprocess_solver_args[1], "%sppcntr.scr", iolib.gscrdr);
	preprocess_solver_args[2]=NULL;

	if (param) iolib.useopt=param->get_i("GAMS Preprocessing optionfile", 0); // number of optionfile

	// Create SBB-control-file
	sprintf(iolib.flnsbbopt, "%sppinfo.scr", iolib.gscrdr);
	ofstream file(iolib.flnsbbopt);
	file << "restart 1" << endl << "rfile ppfl000.scr" << endl;
	file.close();

	sprintf(iolib.flnsbbsav, "ppfl000.scr");
	sprintf(iolib.flnsta, "%sppsta.scr", iolib.gscrdr);
	sprintf(iolib.flnsol, "%sppsol.scr", iolib.gscrdr);
	if (iolib.useopt<2) sprintf(iolib.flnopt, "%s%s.opt", iolib.gwrkdr, (char*)preprocess_solver);
	else if (iolib.useopt<10) sprintf(iolib.flnopt, "%s%s.op%i", iolib.gwrkdr, (char*)preprocess_solver, iolib.useopt);
	else sprintf(iolib.flnopt, "%s%s.o%i", iolib.gwrkdr, (char*)preprocess_solver, iolib.useopt);
/*
	if (iolib.ilog==0 || iolib.ilog==1 || iolib.iolg==3) { }
	else {
		sprintf(iolib.flnlog, "%spplog.scr", iolib.gscrdr);
		iolib.ilog=2;
	}
*/
	iolib.ignbas=1;

	out_log << "GAMSPreprocessing initialized: " << endl
		<< "\t Solver: " << preprocess_solver << endl
		<< "\t Optionfile: " << iolib.flnopt << endl
		<< "\t Shared variables: " << var_pp2LaGO.size() << endl
    << "\t Variables, LaGO has more: " << var_LaGO_notpp.size() << endl
		<< "\t discrete variables: " << i_discr.size() << endl;

	delete preprocess_solver;
	preprocess_iolib=iolib;
	iolib=LaGO_iolib;
	return 0;
}

int GAMSPreprocessing::run(UserVector<double>& x) {
	for (map<int,int>::iterator it(var_pp2LaGO.begin()); it!=var_pp2LaGO.end(); ++it)
		start[it->first]=x(it->second);

	for (int i=0; i<i_discr.size(); ++i)
		fix_discr[i]=lower[i_discr[i]]=upper[i_discr[i]]=start(i_discr[i]);

	iolib=preprocess_iolib;
	if (cioSBBSave(1,1,0,0,(Pointer<double>)fix_discr,(Pointer<double>)fix_discr,(Pointer<double>)start,NULL,NULL,NULL)) {
		out_err << "GAMSPreprocessing: Could not write S/R file." << endl;
		return -1;
	}

	iolib.SBBFlag=recent_calls ? 2 : 1;
	recent_calls=true;

	gfWriteCntr(preprocess_solver_args[1], &iolib, 30);

	int ret=start_process(iolib.line3[subsolver], preprocess_solver_args, 300);
	if (ret) {
		out_err << "Call of GAMS NLP solver " << iolib.line3[subsolver] << " failed with return code " << ret << endl;
		iolib=LaGO_iolib;
		return -1;
	}
	if (remove(iolib.flnsta)) {
		out_err << "Could not remove file " << iolib.flnsta << endl;
		out_err << "Call of GAMS NLP solver " << iolib.line3[subsolver] << " failed with return code " << ret << endl;
		iolib=LaGO_iolib;
		return -1;
	}

	double res, opt_val;
	int dom_violations, iter;
	int model_status, solver_status;
	gfrsol(iolib.flnsol, &model_status, &solver_status, &iter, &res, &opt_val, &dom_violations,
		(Pointer<double>)con_val, (Pointer<double>)con_duals, NULL, (Pointer<int>)con_basind,
		(Pointer<int>)con_type, (Pointer<double>)rhs,
		(Pointer<double>)sol_point, (Pointer<double>)var_duals, NULL, (Pointer<int>)var_basind,
		(Pointer<double>)lower, (Pointer<double>)upper, iolib.nrows, iolib.ncols);

	out_log << "Model status: " << model_status << "\t Solver status: " << solver_status << endl;

	for (map<int,int>::iterator it(var_pp2LaGO.begin()); it!=var_pp2LaGO.end(); ++it)
		x[it->second]=sol_point(it->first);

	// set variables of x, not appearing in preprocessing model, to 0 projected onto box.
	for (list<pair<int,double> >::iterator it(var_LaGO_notpp.begin()); it!=var_LaGO_notpp.end(); ++it)
		x[it->first]=it->second;

	iolib=LaGO_iolib;

	if (solver_status==1 && model_status==4) return 1; // infeasible (in triangular part)
//if (solver_status==1 && model_status==5) return 1; // locally infeasible
	return 0;
}


extern "C" LocOptPreprocessing* create(const Pointer<Param>& param, struct dictRec* LaGO_dict, const MinlpProblem& LaGO_prob) {
	GAMSPreprocessing* pp=new GAMSPreprocessing(param);
	int ret=pp->init(LaGO_dict, LaGO_prob);
	if (ret) {
		out_err << "GAMSPreprocessing::init returned " << ret << endl;
		delete pp;
		return NULL;
	}
	return pp;
}

extern "C" void destroy(LocOptPreprocessing* pp) {
	if (pp) delete (GAMSPreprocessing*)pp;
}
