// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "gams.h"

#ifdef COIN_HAS_GAMSIO

#include "minlpopt.h"
#include <dlfcn.h>

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
#include "clicelib.h"

#ifdef FNAME_LCASE_DECOR
#define G2DINTERVAL0X g2dinterval0x_
#define G2DINTERVAL1X g2dinterval1x_
#else
#define G2DINTERVAL0X g2dinterval0x
#define G2DINTERVAL1X g2dinterval1x
#endif

extern "C" {
void G2D_CALLCONV G2DINTERVAL0X(double xmin[], double xmax[], double* gmin, double* gmax, double vmin[], double vmax[], unsigned int instr[], const int *numIns, double nlCons[]);
void G2D_CALLCONV G2DINTERVAL1X(double xmin[], double xmax[], double dfdxmin[], double dfdxmax[], double vmin[], double vmax[], double vbarmin[], double vbarmax[], unsigned int instr[], const int *numIns, double nlCons[]);

#ifdef COIN_HAS_GAMSCPLEXLICE
#include "gamscplexlice.h"
#endif

//GDX client
#ifdef GDX_AVAILABLE
#include "gdxwrap.h"
#endif
}

extern bool snoptlicenceok;

char* gams::getRowName(int gi, char *name, int bufLen) {
  char
    quote,
    *targ,
    *end,
    *s,
    tbuf[32];
  int
    uelIndices[10],
    nIndices,
    lSym;

  if (gi<0 || gi>=gcdNRows(dict)) return NULL;

  if (gcdRowUels (dict, gi, &lSym, uelIndices, &nIndices) != 0) return NULL;
  if ((s=gcdSymName (dict, lSym, tbuf, 32)) == NULL) return NULL;

  targ = name;
  end = targ + bufLen;

  while (targ < end && *s != '\0') *targ++ = *s++;
  if (targ >= end) return NULL;

  if (0 == nIndices) {
    *targ = '\0';
    return name;
  }

  *targ++ = '(';
  if (targ >= end) return NULL;
  for (int k=0;  k<nIndices;  k++) {
    s = gcdUelLabel (dict, uelIndices[k], tbuf, 32, &quote);
    if (NULL == s) return NULL;

    if (' ' != quote) {
      *targ++ = quote;
      if (targ >= end) return NULL;
    }
    while (targ < end && *s != '\0') *targ++ = *s++;
    if (targ >= end) return NULL;
    if (' ' != quote) {
      *targ++ = quote;
      if (targ >= end) return NULL;
    }
    *targ++ = ',';
    if (targ >= end) return NULL;
  }
  *(targ-1) = ')';
  *targ = '\0';

  return name;
}

char* gams::getColName (int gj, char *name, int bufLen) {
  char
    quote,
    *targ,
    *end,
    *s,
    tbuf[32];
  int
    uelIndices[10],
    nIndices,
    lSym;

  if (gj < 0 || gj >= gcdNCols(dict)) return NULL;

  if (gcdColUels (dict, gj, &lSym, uelIndices, &nIndices) != 0) return NULL;
  if ((s=gcdSymName (dict, lSym, tbuf, 32)) == NULL) return NULL;

  targ = name;
  end = targ + bufLen;

  while (targ < end && *s != '\0') *targ++ = *s++;
  if (targ >= end) return NULL;

  if (0 == nIndices) {
		*targ = '\0';
    return name;
  }

  *targ++ = '(';
  if (targ >= end) return NULL;
  for (int k = 0;  k < nIndices;  k++) {
    s = gcdUelLabel (dict, uelIndices[k], tbuf, 32, &quote);
    if (NULL == s) return NULL;

    if (' ' != quote) {
      *targ++ = quote;
      if (targ >= end) return NULL;
    }
    while (targ < end && *s != '\0') *targ++ = *s++;
    if (targ >= end) return NULL;
    if (' ' != quote) {
      *targ++ = quote;
      if (targ >= end) return NULL;
    }
    *targ++ = ',';
    if (targ >= end) return NULL;
  }
  *(targ-1) = ')';
  *targ = '\0';

  return name;
}

gams* gamsptr=NULL;

void gams::init_cplex_licence(int connr, int varnr, int nnz, int nlnz, int ndisc) {
#ifdef COIN_HAS_GAMSCPLEXLICE
	licenseInit_t initType;
	if (gamscplexlice(connr, varnr, nnz, nlnz, ndisc, 1, &initType, NULL, NULL, NULL, NULL, NULL, NULL)) {
		out_err << "Could not initialize CPLEX license" << endl;
//		exit(-1);
  } else out_log << "GAMS/CPLEX license accepted: " << initType << endl;
#endif
}

void gams::init_snopt_licence() {
#ifdef SNOPT_AVAILABLE
	int expiredSince, bestIdx, mDemo, nDemo, nzDemo, nlnzDemo, discDemo;
	char *sn="SN";
	expiredSince = gcxLiceExp(1, sn, iolib.age, &bestIdx);
//expiredSince=0;
	if (expiredSince >=0 && expiredSince<=29) {
		out_out << "Found GAMS/SNOPT license" << endl;
		snoptlicenceok=true;
	} else {
//		out_log << "SNOPT demo license" << endl;
		snoptlicenceok=false;
/*		gcxGetDemoLimits (&mDemo, &nDemo, &nzDemo, &nlnzDemo, &discDemo);
		if (iolib.nrows>mDemo || iolib.ncols>nDemo || iolib.nnz>nzDemo || iolib.nlnz>nlnzDemo || iolib.ndisc>discDemo)
			out_log << "Model too big for SNOPT demo license" << endl;
*/    }
#endif	
}

gams::gams(Pointer<Param> param_)
: param(param_), obj_sign(1), is_minimization(true), reformed(true), dict(NULL)
{ assert(!gamsptr);
	gamsptr=this;
#ifdef GDX_AVAILABLE
	char* errormsg=new char[128];
	if (gdxWrapInit(NULL, errormsg, 128)) {
	 	cerr << "Could not load GDX I/O library: " << errormsg << endl;
		exit(-1);
	}
	delete[] errormsg;
#endif
}

Pointer<MinlpProblem> gams::get_problem(char* gamsfile) {
	assert(gamsfile);

	Pointer<MinlpProblem> prob(new MinlpProblem());

// 	out_log << "Calling gfinit" << endl;
	gfinit();
	
	out_log << "Reading control file " << gamsfile << endl;
	cntrec info;
	gfrcnt(true, /* because we can handle discrete variables */
		true, /* because we can handle nonlinearities */
		&info,
		gamsfile
	);
	iolib.instrStartIdx=1;
	
// 	out_log << "Opening status file" << endl;
	gfopst();

	if (iolib.useopt) {
		out_log << "Reading parameter file " << iolib.flnopt << endl;
		param->add_file(iolib.flnopt);
		param->read();
	}
	written_gdx_limit=param ? param->get_i("GAMS write solution candidates", 0) : 0;

	init_cplex_licence(iolib.nrows, iolib.ncols, iolib.nnz, iolib.nlnz, iolib.ndisc);
	init_snopt_licence();

	int numrow=info.kgv[1];
	int numcol=info.kgv[2];
	if (info.lgv[3]) obj_sign=1;
	else {
		obj_sign=-1;
		is_minimization=false;
	}
	objcon=info.kgv[15]-1;
	objvar=info.kgv[7];

	out_out << "Rows: " << numrow << endl;
	out_out << "Columns: " << numcol << endl;
/*	out_log << "Objective variable: " << objvar << endl;
	out_log << "Objective row: " << objcon << endl;
	out_log << "Nonzeros in objective: " << info.kgv[14] << endl;
	out_log << "Objective sign: " << obj_sign << endl;*/
  
  	if (numrow>5000 || numcol>5000) {
  		out_err << "This version of LaGO is limited to max. 5000 constraints and max. 5000 variables. Aborting..." << endl;
  		exit(-1);
  	}

	reformed&=(info.kgv[14]==1); // reformable, when objective var appears nonlinear in objective

// 	out_log << "Reading row informations" << endl;
	dvector constants(numrow); // constants of constraints
	con_type.resize(numrow); // type of constraint: eq, leq, geq, free
	rhs.resize(numrow);
	rowrec rowdata;
	for (int c=0; c<numrow; c++) {
		gfrrow(&rowdata);
		constants[c]=-rowdata.rdata[2];
		rhs[c]=rowdata.rdata[2];
		con_type[c]=rowdata.idata[1];
		if (con_type[c]==1) constants[c]*=-1;
		if (con_type[c]==3) {
			out_err << "Free constraints not supported" << endl;
			exit(-1);
		}
	}

	reformed&=(con_type[objcon]==0); // objective-row is an equality

	if (iolib.dictFileWritten) {
// 		out_log << "Reading dictionary file" << endl;
    int ret = gcdLoad(&dict, iolib.flndic, iolib.dictVersion);
		if (ret) {
			out_err << "Error reading dictionary file." << endl;
			dict=NULL;
		}
	}

	vector<Pointer<UserVector<double> > > linear(numrow);
	for (int c=0; c<numrow; c++) linear[c]=new SparseVector<double>(numcol);
	vector<Pointer<SparsityInfo> > sparsity(numrow);
	for (int c=0; c<numrow; c++) {
		sparsity[c]=new SparsityInfo();
		sparsity[c]->linear=new map<int, SparsityInfo::LinearVariable>;
		sparsity[c]->nonlinear=new map<int, SparsityInfo::NonlinearVariable>;
	}

// 	out_log << "Reading columns" << endl;
	lower.resize(numcol);
	upper.resize(numcol);
	colrec coldata;
	char* namebuf=new char[50];
	char* name=NULL;
	for (int i=0; i<numcol; i++) {
		gfrcol(&coldata);
		if (coldata.idata[3]>2) {
			out_err << "Variable " << i << " has type " << coldata.idata[3] << endl;
			out_err << "This kind of variables is not supported, aborting." << endl;
			exit(-1);
		}
		double low=coldata.cdata[1]; lower[i]=low;
		double up=coldata.cdata[3]; upper[i]=up;
		if (low==info.xgv[11]) low=-INFINITY;
		if (up==info.xgv[10]) up=INFINITY;
		if (i==objvar) { // objective variable; substituted later
			reformed&=(low==-INFINITY) && (up==INFINITY); // objective variables is free
		}

		if (dict && !(name=getColName(i, namebuf, 50)))
			out_err << "Couldn't retrieve name for variable " << i << endl;
		prob->add_var(i, 0, coldata.idata[3], low, up, name);
		prob->primal_point[i]=coldata.cdata[2];
//		out_out << "Added variable " << i << " type " << coldata.idata[3] << " low: " << low << " up: " << up << " start: " << prob->primal_point(i) << endl;

		double coeff;
		int index;
		int nltyp;
		for (int r=0; r<coldata.idata[1]; r++) {
			gfrcof(&coeff, &index, &nltyp);
			if (con_type[index-1]==1) coeff*=-1; // geq-constraint -> change to leq-constraint
//			out_log << "Row " << index << " coeff=" << coeff << " nltyp=" << nltyp << endl;
			if (nltyp==0) {// linear variable
				linear[index-1]->SetElement(i, .5*coeff);
			} else
				sparsity[index-1]->nonlinear->insert(pair<int, SparsityInfo::NonlinearVariable>(i, SparsityInfo::NonlinearVariable()));
			if (index-1==objcon && i==objvar) reformed&=(nltyp==0); // objective variable appears linear
		}
	}
	out_out << "Discrete variables: " << prob->i_discr.size() << endl;

	out_log << "Reformable: " << (reformed ? "yes" : "no") << endl;

	if (reformed) {
		prob->lower[objvar]=prob->upper[objvar]=prob->primal_point[objvar]=0.;
		if (prob->i_discr.size()==prob->dim()-1 && (!prob->discr[objvar])) {
			prob->discr[objvar]=true;
			prob->i_discr.push_back(objvar);
			prob->i_cont.clear();
		}
	}

	for (int c=0; c<numrow; c++) {
		if (dict && !(name=getRowName(c, namebuf, 50)))
			out_err << "Couldn't retrieve name for row " << c << endl;
		if ((c==objcon) && name) objcon_name=strdup(name);
		if ((c==objcon) && reformed) { // reforming
			prob->add_obj(new SepQcFunc(NULL, linear[c], NULL, constants[c]));
			obj_sign*=-2*(*prob->obj->b[0])(objvar);
			(*prob->obj->b[0])[objvar]=0.;
			if (obj_sign!=1) {
				assert(obj_sign!=0);
				*prob->obj->b[0]/=obj_sign;
				prob->obj->c/=obj_sign;
			}
		} else
			prob->add_con(new SepQcFunc(NULL, linear[c], NULL, constants[c]), con_type[c]==0, name);
	}
	if (!reformed) {
		Pointer<SparsityInfo> si(new SparsityInfo);
		si->linear=new map<int, SparsityInfo::LinearVariable>;
		si->linear->insert(pair<int, SparsityInfo::LinearVariable>(objvar, SparsityInfo::LinearVariable(obj_sign)));
		prob->add_obj(new SepQcFunc(NULL, new SparseVector<double>(prob->dim(), objvar, .5*obj_sign), NULL, 0, si));
	}

	delete[] namebuf;

// 	out_log << "G2DINIT()" << endl;
	G2DINIT();
	G2DSTATSINIT(gci2DData);
//	G2DSETMSGCALLBACK(MSGCB);

	Pointer<gamsNLData> data(new gamsNLData);
// 	out_log << "Calling gfrins" << endl;
	unsigned int* instr=new unsigned int[4*info.knv[4]];
	double* nlwork=new double[info.knv[5]];
	int* startIdx=new int[numrow+1];
	int* numInstr=new int[numrow+1];
	gfrins(instr, nlwork);

// 	out_log << "Calling gciGetInstrStarts" << endl;
	gciGetInstrStarts(instr, startIdx, 1, numInstr);

// 	out_log << "Calling G2DFUNCEXTRACT and G2DADDRESSPATCH" << endl;
	data->instr=new unsigned int[iolib.lenins];
	data->startIdx=new int[numrow+1];
	data->nlCons=new double[iolib.ncons];
	data->numInstr=new int[numrow+1];
	memcpy(data->nlCons, nlwork, iolib.ncons*sizeof(double));

	int check=1;
	int nr_nonlincon=0;
	for (int i=0; i<numrow; i++)
		if (startIdx[i]) {
			data->startIdx[i]=data->lenins;
			G2DFUNCEXTRACT(instr+startIdx[i]-1, data->instr+data->startIdx[i]-1, &(numInstr[i]), &(data->numInstr[i]));
			data->lenins+=data->numInstr[i];
			G2DADDRESSPATCH(data->instr+data->startIdx[i]-1, &(data->numInstr[i]), &check);
			data->maxins=MAX(data->maxins, data->numInstr[i]);
			int rss=G2DGETRESSTACKSIZE(data->instr+data->startIdx[i]-1, &(data->numInstr[i]));
			data->resstacksize=MAX(data->resstacksize, rss);
			if ((i==objcon) && reformed) {
				if (obj_sign==1) prob->obj->s[0]=new gamsFunc(numcol, i, data, sparsity[i]);
				else if (obj_sign==-1) prob->obj->s[0]=new MinusFunc(new gamsFunc(numcol, i, data, sparsity[i]));
				else prob->obj->s[0]=new SumFunc(new gamsFunc(numcol, i, data, sparsity[i]), NULL, 1/	obj_sign, 0);
				prob->obj->set_curvature(0, Func::UNKNOWN);
			} else {
				int index=i-((i>objcon && reformed) ? 1 : 0);
				if (con_type[i]==1) prob->con[index]->s[0]=new MinusFunc(new gamsFunc(numcol, i, data, sparsity[i]));
				else prob->con[index]->s[0]=new gamsFunc(numcol, i, data, sparsity[i]);
				prob->con[index]->set_curvature(0, Func::UNKNOWN);
			}
			++nr_nonlincon;
		} else {
			data->startIdx[i]=0;
			data->numInstr[i]=0;
		}
	out_log << "Nonlinear constraints: " << nr_nonlincon << endl;

	data->resstack=new double[data->resstacksize];
	data->s=new double[data->lenins];
	data->sbar=new double[data->lenins];

	delete[] instr;
	delete[] nlwork;
	delete[] startIdx;
	delete[] numInstr;

	return prob;
}

gams::~gams() {
	if (dict) gcdFree(dict);
	gfclos();
	gamsptr=NULL;
}

void gams::write_matlab(const dvector& x, const char* filename, vector<Pointer<char> >& var_names) {
	ofstream out(filename);
	out.precision(20);
	for (int i=0; i<x.dim(); i++)
		out << var_names[i] << " = \t" << x(i) << ";" << endl;
	out.close();
}

void gams::write_gams(const dvector& x, const char* filename, const vector<bool>& discr) {
  char quote, *targ, *end, *s, tbuf[32];
  int uelIndices[10], nIndices, lSym, oldSym=-1;
	ofstream out(filename);
	out.precision(20);
	out << "$offdigit" << endl;
	for (int i=0; i<x.dim(); ++i) {
  	if (gcdColUels(dict, i, &lSym, uelIndices, &nIndices) != 0) {
			out_err << "Error in gcdColUels variable " << i << endl;
			continue;
		}
//		out_log << i << " " << lSym << " " << nIndices << "\t "; for (int k=0; k<nIndices; ++k) out_log << uelIndices[k] << " ";
//		out_log << "\t --> " << gcdColIndex(dict, lSym, nIndices, uelIndices, 0) << endl;
	  if ((s=gcdSymName(dict, lSym, tbuf, 32)) == NULL) {
			out_err << "Error in gcdSymName variable " << i << endl;
			continue;
		}

		out << s << (discr[i] ? ".fx" : ".l");

		if (nIndices>0) out << '(';
		for (int k = 0;  k < nIndices;  k++) {
    	if ((s=gcdUelLabel(dict, uelIndices[k], tbuf, 32, &quote))==NULL) {
				out_err << "Error in gcdUelLabel variable " << i << " uel " << k << endl;
				continue;
			}
			out << '\'' << s << '\'';
//			if (' ' != quote) out_log << quote; out << s;	if (' ' != quote) out_log << quote;
	    if (k<nIndices-1) out << ',';
		}
		if (nIndices>0) out << ')';
		out << " = " << x(i) << ";" << endl;
  }
	out.close();
}

#ifdef GDX_AVAILABLE
void gams::gdx_error(int n) {
//	int n = GDXGetLastError(gdxio);
	char* mess=new char[256];
	if (!GDXErrorStr(n,&mess)) cout << "Could not retrieve GDXIO error message" << endl;
	else cout << mess << endl;

//	if (mess) delete mess;
	exit(1);
}

void gams::write_gdx(const dvector& x, char* filename, double val) {
	if (written_gdx_limit) {
		if (written_gdx.size()>=written_gdx_limit) {
			if (is_minimization) { // minimization problem, remove from end
				remove(written_gdx.rbegin()->second);
				written_gdx.erase(--written_gdx.end());
			} else { // maximization problem, remove from front
				remove(written_gdx.begin()->second);
				written_gdx.erase(written_gdx.begin());
			}
		}
		written_gdx.insert(pair<double, Pointer<char> >(val, strdup(filename)));
	}

	write_gdx(x, filename);	
}

void gams::write_gdx(const dvector& x, char* filename) {
	PGXFile gdxhandle=NULL;
	int errornr;
	if (errornr=GDXOpenWrite(&gdxhandle, filename, "LaGO")) gdx_error(errornr);
//	out_log << "GDX file open." << endl;

  char quote, *targ, *end, *s, tbuf[32];
  int uelIndices[10], nIndices, lSym, oldSym=-1;
	TgdxStrIndex Elements;
	for (int i=0; i<10; ++i) Elements[i]=new char[32];
	TgdxValues Values;
	for (int i=0; i<5; ++i) Values[i]=0.;
	for (int i=0; i<x.dim(); ++i) {
  	if (gcdColUels(dict, i, &lSym, uelIndices, &nIndices) != 0) return;
	  if ((s=gcdSymName(dict, lSym, tbuf, 32)) == NULL) return;

		if (oldSym!=lSym) { // we start a new symbol
			if (oldSym!=-1) { // close old symbol, if we are not at the first; test i==0 should do the same
				if (!GDXDataWriteDone(gdxhandle)) {
					out_err << "Error in GDXDataWriteDone for symbol nr. " << oldSym << endl;
					exit(-1);
				}

//				out_log << endl;
			}
			if (!GDXDataWriteStrStart(gdxhandle, s, NULL, nIndices, dt_var, 0)) {
				out_err << "Error in GDXDataWriteStrStart for symbol " << s << endl;
				exit(-1);
			}
//			out_log << "Start writing symbol " << s;
			oldSym=lSym;
		}

		for (int k = 0;  k < nIndices;  k++) {
    	if ((s=gcdUelLabel(dict, uelIndices[k], tbuf, 32, &quote))==NULL) {
				out_err << "Error in gcdUelLabel." << endl;
				exit(-1);
			}
			strcpy(Elements[k], tbuf);
//			out_log << " " << tbuf;
		}
		Values[0]=x(i);
//		out_log << " = " << x(i) << "\t ";

		if (!GDXDataWriteStr(gdxhandle, Elements, Values)) {
			out_err << "Error in GDXDataWriteStr for symbol " << s << endl;
			exit(-1);
		}
  }
//	out_log << endl;
	if (!GDXDataWriteDone(gdxhandle)) { // finish last symbol
		out_err << "Error in GDXDataWriteDone for symbol nr. " << oldSym << endl;
		exit(-1);
	}

  if (errornr=GDXClose(&gdxhandle)) gdx_error(errornr);

	for (int i=0; i<10; ++i) delete Elements[i];
}

void gams::read_gdx(dvector& x, char* filename) {
	PGXFile gdxhandle=NULL;
	int errornr;
	if (errornr=GDXOpenRead(&gdxhandle, filename)) gdx_error(errornr);
//	out_log << "GDX file open." << endl;

	int nrsy=0, nruel=0;
	if (GDXSystemInfo(gdxhandle, &nrsy, &nruel)==0) {
		out_err << "Error retrieving symbol infos." << endl;
		exit(-1);
	}
//	out_log << "Symbols: " << nrsy << "\t Unique elements: " << nruel << endl;

  int uelIndices[10]; // local uel indices
	char* name=new char[255]; // buffer for symbol name
	int dim, type;
	int recCount; // entries per symbol
	TgdxStrIndex Elements; // uel names
	for (int i=0; i<10; ++i) Elements[i]=new char[32];
	TgdxValues Values;
	int afdim; // index of first changed uel
	int symIndex; // symbol index in our dictionary
	int varIndex=-1; // last variable index
	for (int i=1; i<=nrsy; ++i) {
		GDXSymbolInfo(gdxhandle, i, name, &dim, &type); // get symbol from gdx file
		if (type!=dt_var) continue; // skip non-variables

		symIndex=gcdSymIndex(dict, name); // get local symbol index
		if (symIndex<0) {
//			out_log << "skipping symbol " << name << endl;
			continue;
		}

		GDXDataReadStrStart(gdxhandle, i, &recCount); // start reading symbol from gdx file

//		out_log << "Symbol " << i << ": " << name << "\t dim: " << dim << "\t Index: " << symIndex << "\t recCount: " << recCount << endl;

		for (int j=1;  j<=recCount; j++) {
			GDXDataReadStr(gdxhandle, Elements, Values, &afdim); // read one entry from gdx file
			for (int k=afdim-1; k<dim; k++) {
				uelIndices[k]=gcdUelIndex(dict, Elements[k]); // get appropriate local uel indices
//				out_log << Elements[k] << " -> " << uelIndices[k] << "\t ";
			}
			varIndex=gcdColIndex(dict, symIndex, dim, uelIndices, varIndex); // get variable index
//			out_log << " ---> " << varIndex << "\t " << Values[0] << endl;
			if (varIndex<0) {
//				out_log << "skipping variable " << name << '(';
//				for (int k=0; k<dim; ++k) out_log << Elements[k] << ",";
//				out_log << ')' << endl;
				continue;
			}
			x[varIndex]=Values[0];
		}

		GDXDataReadDone(gdxhandle);
//		out_log << endl;
	}
  if (errornr=GDXClose(&gdxhandle)) gdx_error(errornr);
	for (int i=0; i<10; ++i) delete Elements[i];
	delete name;
}
#endif // GDX_AVAILABLE

void gams::write_sol_file(const dvector& sol_point, int model_status, int solver_status, int iter, double time, Pointer<MinlpProblem> prob) {
	if (model_status>2) { // not solved
		gfwsolflag(0); // not writing a solution
		gfwsta(model_status, solver_status, iter, time, 0, 0);
		return;
	}

	out_out << "Calling GAMS-NLP-solver to get dual values and nice return codes" << endl;
	Timer t;
	Pointer<gamsLocOpt> solver(new gamsLocOpt(prob, param, NULL, NULL));
	dvector start(sol_point);
	start[objvar]=prob->obj->eval(sol_point);
	int ret=solver->solve(start);
	time+=t.stop();
	bool feasible=!prob->feasible(solver->sol_point, 1E-4);
	out_out << "GAMS solver return: " << ret << " Feasible point: " << (feasible ? "yes" : "no") << endl;
	solver->sol_point[objvar]=prob->obj->eval(solver->sol_point);
	if (!feasible) {
		gfwsolflag(0);
		gfwsta(solver->model_status, solver->solver_status, iter, time, 0, 0);
/*
		if (!ret) {
#ifdef GDX_AVAILABE
			write_gdx(solver->sol_point, "point.gdx", solver->opt_val());
#endif
			write_gams(solver->sol_point, "point.gms", prob->discr);
			write_matlab(solver->sol_point, "point.m", prob->var_names);
		}*/

		return;
	}
	out_out << "Objective value: " << solver->opt_val() << endl;
	gfwsta(solver->model_status, solver->solver_status, iter, time, solver->opt_val(), 0);

	for (int c=0; c<iolib.nrows; c++)
		gfwrow(solver->con_val(c), solver->duals_con(c), solver->basind_con(c), 0);

	for (int i=0; i<iolib.ncols; i++) {
		if (prob->discr[i]) solver->basind_var[i]=3; // superbasic for discrete variables
		gfwcol(solver->sol_point(i), solver->duals_var(i), solver->basind_var(i), 0);
	}
/*
#ifdef GDX_AVAILABLE
	write_gdx(solver->sol_point, "point.gdx", solver->opt_val());
#endif
	write_gams(solver->sol_point, "point.gms", prob->discr);
	write_matlab(solver->sol_point, "point.m", prob->var_names);*/
}

void gams::write_sol_set(const set<SolCandidate>& sol_set) {
	char* name=new char[50];
	int i=0;
	for(set<SolCandidate>::const_iterator it(sol_set.begin()); it!=sol_set.end(); ++it, ++i) {
#ifdef GDX_AVAILABLE
		sprintf(name, "solcand%.4i_%f.gdx", i, it->first);
		write_gdx(it->second, name, it->first);
#else
		sprintf(name, "solcand%.4i_%f.gms", i, it->first);
		write_gams(it->second, name, vector<bool>(it->second.dim(), false));
#endif
	}
	delete name;
}

void gams::write_box(const dvector& lower, const dvector& upper) {
  char quote, *targ, *end, *s, tbuf[32], *name;
  int uelIndices[10], nIndices, lSym, oldSym=-1;
	ofstream out("box.gms");
	out.precision(20);
	out << "$offdigit" << endl;
	for (int i=0; i<lower.dim(); ++i) {
  	if (gcdColUels(dict, i, &lSym, uelIndices, &nIndices) != 0) {
			out_err << "Error in gcdColUels variable " << i << endl;
			continue;
		}
//		out_log << i << " " << lSym << " " << nIndices << "\t "; for (int k=0; k<nIndices; ++k) out_log << uelIndices[k] << " ";
//		out_log << "\t --> " << gcdColIndex(dict, lSym, nIndices, uelIndices, 0) << endl;
	  if ((s=gcdSymName(dict, lSym, tbuf, 32)) == NULL) {
			out_err << "Error in gcdSymName variable " << i << endl;
			continue;
		}
		name=strdup(s);

		out << name << (lower(i)==upper(i) ? ".fx" : ".lo");

		if (nIndices>0) out << '(';
		for (int k = 0;  k < nIndices;  k++) {
    	if ((s=gcdUelLabel(dict, uelIndices[k], tbuf, 32, &quote))==NULL) {
				out_err << "Error in gcdUelLabel variable " << i << " uel " << k << endl;
				continue;
			}
			out << '\'' << s << '\'';
	    if (k<nIndices-1) out << ',';
		}
		if (nIndices>0) out << ')';
		out << " = " << lower(i) << ";" << endl;

		if (lower(i)!=upper(i)) {
			out << name << ".up";

			if (nIndices>0) out << '(';
			for (int k = 0;  k < nIndices;  k++) {
    		if ((s=gcdUelLabel(dict, uelIndices[k], tbuf, 32, &quote))==NULL) {
					out_err << "Error in gcdUelLabel variable " << i << " uel " << k << endl;
					continue;
				}
				out << '\'' << s << '\'';
		    if (k<nIndices-1) out << ',';
			}
			if (nIndices>0) out << ')';
			out << " = " << upper(i) << ";" << endl;
		}
		delete name;
	}
	out.close();

}

// ---------------------------------- gamsFunc ---------------------------------

double gamsFunc::eval(const double* x) const {
	double val;
	int numerr=0;
	double dt;

	G2DFUNCEVAL0((double*)x, &val, data->s,
		data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons,
		&dt, &numerr);

	if (numerr) {
		data->domain_violations+=numerr;
		out_err << "Error evaluating constraint ";
		if (gamsptr->dict) {
			char* namebuf=new char[50];
			out_err << gamsptr->getRowName(connr, namebuf, 50) << endl;
			delete namebuf;
		} else out_err << connr << endl;
		return log(-1.);
	}

	return val;
}

int gamsFunc::valgrad(double& val, double* y, const double* x) const {
	int numerr=0;
	double dt;
//	out_log << "Calling G2DFUNCEVAL0" << endl;
	G2DFUNCEVAL0((double*)x, &val, data->s,
		data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons,
		&dt, &numerr);
	if (numerr) {
		data->domain_violations+=numerr;
		out_err << "Error evaluating constraint ";
		if (gamsptr->dict) {
			char* namebuf=new char[50];
			out_err << gamsptr->getRowName(connr, namebuf, 50) << endl;
			delete namebuf;
		} else out_err << connr << endl;
//		char* namebuf=new char[50];
//		for (int i=0; i<dim(); ++i)
//			out_err << gamsptr->getColName(i, namebuf, 50) << "\t " << x[i] << endl;
		return numerr;
	}

//	out_log << "Calling G2DREVERSEEVAL1" << endl;
	numerr=0;
	G2DREVERSEEVAL1((double*)x, y, data->s, data->sbar,
		data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons,
		&dt, &numerr);
	data->domain_violations+=numerr;
	if (numerr) {
		out_err << "Error computing gradient of constraint ";
		if (gamsptr->dict) {
			char* namebuf=new char[50];
			out_err << gamsptr->getRowName(connr, namebuf, 50) << endl;
			delete namebuf;
		} else out_err << connr << endl;
	}

	return numerr;
}

void gamsFunc::HessMult(double* y, const double* x, const double* z) const {
	double val, gradvecprod; // to store value and gradient*z
	int numerr=0;
	G2DDIR2DX((double*)x, (double*)z, y, &val, &gradvecprod,
		data->resstack, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons, &numerr);
	if (numerr) {
		data->domain_violations+=numerr;
//		out_err << "Error computing hessian-vector-multiplication for constraint ";
//		if (gamsptr->dict) {
//			char* namebuf=new char[50];
//			out_err << gamsptr->getRowName(connr, namebuf, 50) << endl;
//			delete namebuf;
//		} else out_err << connr << endl;
		return;
	}
}

#ifdef FILIB_AVAILABLE
interval<double> gamsFunc::eval(const IntervalVector& x) const {
	double* xmin=new double[x.dim()];
	double* xmax=new double[x.dim()];
	double valmin, valmax;
	for (int i=0; i<x.dim(); i++) {
		xmin[i]=x(i).inf();
		xmax[i]=x(i).sup();
	}

//	out_log << "Calling G2DINTERVAL0X " << connr << endl;
	G2DINTERVAL0X(xmin, xmax, &valmin, &valmax, data->s, data->sbar, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons);

	delete xmin;
	delete xmax;
	
	return interval<double>(valmin, valmax);
}

int gamsFunc::valgrad(interval<double>& val, IntervalVector& y, const IntervalVector& x) const {
	double* xmin=new double[x.dim()];
	double* xmax=new double[x.dim()];
	double valmin, valmax;
	double* ymin=new double[y.dim()];
	double* ymax=new double[y.dim()];
	for (int i=0; i<x.dim(); i++) {
		xmin[i]=x(i).inf();
		xmax[i]=x(i).sup();
		ymin[i]=0.;
		ymax[i]=0.;
	}

//	double inf=filib::fp_traits<double>::infinity();
//	G2DSETSLINFY(&inf);

//	out_log << "Calling G2DINTERVAL0X " << connr << endl;
	G2DINTERVAL0X(xmin, xmax, &valmin, &valmax, data->s, data->sbar, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons);

	double* s2=new double[data->lenins];
	double* sbar2=new double[data->lenins];

	G2DINTERVAL1X(xmin, xmax, ymin, ymax, data->s, data->sbar, s2, sbar2, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons);

	val=interval<double>(valmin, valmax);
	for (int i=0; i<y.dim(); i++) {
		if (ymin[i]==-1E+20) ymin[i]=filib::fp_traits<double>::ninfinity();
		if (ymax[i]==1E+20) ymax[i]=filib::fp_traits<double>::infinity();
		y.SetElement(i, interval<double>(ymin[i], ymax[i]));
	}

	delete xmin;
	delete xmax;
	delete ymin;
	delete ymax;
	delete s2;
	delete sbar2;

	return 0;
}
#endif

// ------------------------------------ gamsLocOpt -----------------------------

void* gamsLocOpt::preprocess_handle=NULL;

gamsLocOpt::gamsLocOpt(Pointer<MinlpProblem> prob_, Pointer<Param> param_, Pointer<ostream> out_solver_p_, Pointer<ostream> out_solver_log_p_)
: LocOpt(prob_->dim(), out_solver_p_, out_solver_log_p_), prob(prob_), param(param_),
  lower_discr(prob_->i_discr.size()), upper_discr(prob_->i_discr.size()),
  lower(gamsptr->lower), upper(gamsptr->upper),
	con_val(iolib.nrows), duals_con(iolib.nrows),
	duals_var(iolib.ncols), basind_con(iolib.nrows, 2), basind_var(iolib.ncols, 2),
  recent_calls(0), model_status(0), solver_status(0), second_run(false), optfile(0), subsolver(-1), iolibsave((void*)new tiolib),
	preprocessargs(NULL),
	preprocessing(NULL)
{	tol=1E-4;
	assert(gamsptr); // a bit late
	tmpsolfn=strdup(iolib.flnsol);
	*(tiolib*)iolibsave=iolib;
	args=new char*[3];
	solvername=param ? param->get("GAMS LocOpt solver", "conopt") : Pointer<char>(strdup("conopt"));
	int slen=strlen(solvername);
	for (int i=0; i<iolib.nosolvers; ++i)
		if (strnicmp(solvername, iolib.line1[i], slen)==0 && isspace(iolib.line1[i][slen])) {
			subsolver=i;
			break;
		}
	if (subsolver<0) {
		out_err << "Could not find GAMS NLP solver " << solvername << endl;
		exit(-1);
	}
	args[0]=strdup(iolib.line3[subsolver]);
	args[1]=new char[1024];
	sprintf(args[1], "%soqcntr.scr", iolib.gscrdr);
	args[2]=NULL;

	if (param) optfile=param->get_i("GAMS LocOpt optionfile", 0); // number of optionfile

	// Create SBB-control-file
	sprintf(iolib.flnsbbopt, "%soqinfo.scr", iolib.gscrdr);
	ofstream file(iolib.flnsbbopt);
	file << "restart 1" << endl << "rfile oqfl000.scr" << endl;
	file.close();

	if (param) {
		Pointer<char> preprocess_lib=param->get("GAMS LocOpt preprocessing");
		if (preprocess_lib) {
			if (!preprocess_handle) preprocess_handle=dlopen(preprocess_lib, RTLD_NOW);
			if (!preprocess_handle) {
				out_err << "Error loading preprocess library: " << dlerror() << endl;

				preprocessargs=new char*[3];
				preprocessargs[0]=preprocess_lib;
				preprocessargs[1]=strcpy(new char[40], "lagostart");
				preprocessargs[2]=NULL;
				preprocess_keepgdx=param->get_i("GAMS LocOpt preprocessing keep gdx",0);
			} else {
				preprocess_create_t* create;
				create=(preprocess_create_t*)dlsym(preprocess_handle, "create");
				preprocessing=create(param, gamsptr->dict, *prob);
			}
		}
	}

	write_solcand=param && param->get_i("GAMS write solution candidates", 0);
	write_startpoint=param && param->get_i("GAMS write startpoint", 0);
}

gamsLocOpt::~gamsLocOpt() {
	free(args[0]);
	delete[] args[1];
	delete[] args;
	delete (tiolib*)iolibsave;
	if (preprocessargs) {
		delete preprocessargs[0];
		delete preprocessargs[1];
		delete preprocessargs;
	}
	if (preprocessing) {
		preprocess_destroy_t* destroy=(preprocess_destroy_t*)dlsym(preprocess_handle, "destroy");
		destroy(preprocessing);
	}
}

int gamsLocOpt::do_preprocessing(const dvector& start) {
	sprintf(preprocessargs[1]+9,"%.10f.gms", rd);
	out_solver_log << "write startpoint to " << preprocessargs[1] << endl;
	gamsptr->write_gams(start, preprocessargs[1], prob->discr);

	preprocessargs[1][21]=0; // skip file extension
	int ret=start_process(preprocessargs[0], preprocessargs, 300); // calling sub-solver

	preprocessargs[1][21]='.'; // add '.' again
	remove(preprocessargs[1]);

	if (ret) return ret;

	for (int i=0; i<prob->dim(); ++i) { // set sol_point to 0, projected onto box
		if (prob->lower(i)<=0. && prob->upper(i)>=0.) sol_point[i]=0.;
		else if (prob->lower(i)>0) sol_point[i]=prob->lower(i);
		else sol_point[i]=prob->upper(i);
	}

	// lagostart....gdx
	preprocessargs[1][23]='d';
	preprocessargs[1][24]='x';
	out_solver_log << "read improved startpoint from " << preprocessargs[1] << endl;
#ifdef GDX_AVAILABLE
	gamsptr->read_gdx(sol_point, preprocessargs[1]);
#else
	out_err << "No GDX support. Cannot read startpoint :-(." << endl;
#endif
	if (!preprocess_keepgdx) remove(preprocessargs[1]);

	return 0;
}

int gamsLocOpt::solve(dvector& start) {
	timer.start();

#ifdef PROJECT_DEBUG
	{
		double maxjac=0.;
		int maxjacconnr=0;
		int maxjacvarnr=0;
		dvector grad(prob->dim());
		for (int c=0; c<prob->con.size(); ++c) {
			prob->con[c]->grad(grad, start);
			for (int i=0; i<grad.dim(); ++i)
				if (fabs(grad(i))>maxjac) {
					maxjac=fabs(grad(i));
					maxjacconnr=c;
					maxjacvarnr=i;
				}
		}
		cout << "Start CONOPT. Biggest absolut element in jacobian for con " << prob->con_names[maxjacconnr] << "\t variable " << prob->var_names[maxjacvarnr]
		<< "\t value of variable: " << start(maxjacvarnr) << "\t value of jacobian: " << maxjac << endl;
	}
#endif
	sprintf(iolib.flnsbbsav, "oqfl000.scr");
	sprintf(iolib.flnsbbopt, "%soqinfo.scr", iolib.gscrdr);

	sol_point=start;

	if (write_startpoint) {
#ifdef GDX_AVAILABLE
			char* name=new char[50];
			double rd=random(0.,1.);
			sprintf(name, "startpoint_%.10f.gdx", rd);
			gamsptr->write_gdx(sol_point, name);
			sprintf(name+24, "viol");
			ofstream violfile(name);
			prob->print_most_violated_constraints(sol_point, violfile, 20);
			delete[] name;
#endif
	}


	if (preprocessing || preprocessargs) {
		rd=random(0.,1.);
		if (prob->feasible(sol_point, tol)) {
			int ret=preprocessing ? preprocessing->run(sol_point) : do_preprocessing(start);
			if (ret) {
				out_solver << "Preprocessing returned " << ret << endl;
				return 1; // preprocessing claims infeasible
			}
		}
	}
	// Create start-file
	int i0;
	for (int i=0; i<prob->i_discr.size(); i++) { // store new bounds for discrete variables
		i0=prob->i_discr[i];
		lower[i0]=lower_discr[i]=upper[i0]=upper_discr[i]=sol_point[i0]=closestint(sol_point(i0));
	}

//	duals_con=0; basind_con=2; basind_var=2;
  // fake dual variables and variable and equation status
	dvector duals(prob->con.size()+1);
	ivector varstatus(prob->dim(),3); 
	ivector equstatus(prob->con.size()+1,3); 
  if (cioSBBSave(1,1,1,1,(Pointer<double>)lower_discr, (Pointer<double>)upper_discr, (Pointer<double>)sol_point, (Pointer<double>)duals, (Pointer<int>)varstatus, (Pointer<int>)equstatus)) {
//  if (cioSBBSave(1,1,0,0,(Pointer<double>)lower_discr, (Pointer<double>)upper_discr, (Pointer<double>)sol_point, NULL, NULL, NULL)) {
		out_err << "Could not write to S/R File." << endl;
		exit(-1);
	}

	iolib.SBBFlag=recent_calls ? 2 : 1;
	iolib.ignbas=1;

	sprintf(iolib.flnsta, "%soqsta.scr", iolib.gscrdr);
	sprintf(iolib.flnsol, "%soqsol.scr", iolib.gscrdr);

	if (iolib.ilog==0 || iolib.ilog==1 || iolib.ilog==3) {
/*		fprintf(gfiolog, "\n");
		fflush(gfiolog);
*/	} else {
		sprintf(iolib.flnlog, "%soqlog.scr", iolib.gscrdr);
		iolib.ilog=2;
	}

	iolib.useopt=optfile;
	sprintf(iolib.flnopt, "%s%s.opt", iolib.gwrkdr, (char*)solvername);

	gfWriteCntr(args[1], &iolib, 30); // write control file

	recent_calls++;
	int ret;
	if (ret=start_process(iolib.line3[subsolver], args, 600)) { // calling sub-solver, timelimit: 10 minutes
		out_err << "Spawn of GAMS NLP solver " << solvername << " failed. exit code: " << ret << endl;
		return -1;
//		exit(-1);
	}

  if (iolib.ilog!=0 && iolib.ilog!=1 && iolib.ilog!=3) {
		ifstream log(iolib.flnlog);
		if (!log.good()) {
			out_err << "Could not open " << iolib.flnlog << " for reading." << endl;
			exit(-1);
		}
		char buf[1024];
		while (!log.eof()) {
			log >> buf;
			fprintf(gfiolog, buf);
		}
		if (remove(iolib.flnlog))
			out_err << "Could not remove file " << iolib.flnlog << endl;
	}
	if (remove(iolib.flnsta)) {
		out_err << "Could not remove file " << iolib.flnsta << endl;
		out_err << "Spawn of GAMS NLP solver " << solvername << " failed." << endl;
		return -1;
//		exit(-1);
	}

	double res; // local resource counter ??
	int dom_violations; // domain violations
	gfrsol(iolib.flnsol, &model_status, &solver_status, &iter_, &res, &opt_val_, &dom_violations,
		(Pointer<double>)con_val, (Pointer<double>)duals_con, NULL, (Pointer<int>)basind_con,
		(Pointer<int>)gamsptr->con_type, (Pointer<double>)gamsptr->rhs,
		(Pointer<double>)sol_point, (Pointer<double>)duals_var, NULL, (Pointer<int>)basind_var,
		(Pointer<double>)lower, (Pointer<double>)upper, iolib.nrows, iolib.ncols);
	if (gamsptr->reformed) sol_point[gamsptr->objvar]=0.;
//	opt_val_*=gamsptr->obj_sign;

	strcpy(iolib.flnsol, tmpsolfn);
	iolib=*(tiolib*)iolibsave;

	out_solver << "Solver status: " << solver_status << "\t Model  status: " << model_status << "\t Objective value: " << opt_val_ << endl;

	timer.stop();

// translating return codes from GAMS to SNOPT
	if (solver_status==1 && (model_status==1 || model_status==2)) {
		if (write_solcand) {
#ifdef GDX_AVAILABLE
			/*if (!preprocessargs) */rd=random(0.,1.);
			char* name=new char[50];
			sprintf(name, "solcand%f_%.10f.gdx", opt_val_, rd);
			gamsptr->write_gdx(sol_point, name, opt_val_);
			delete[] name;
#endif
		}
		out_solver << "feasiblity: "; prob->feasible(sol_point, 1E-4, out_solver_p);
		if ((!second_run) && prob->feasible(sol_point, tol, out_solver_log_p)) {
			out_solver << "Point infeasible in our measure of feasiblity. Starting again from this point." << endl;
			second_run=true;
			return solve(sol_point);
		}
		if (second_run) second_run=false;
		return 0; // (locally) optimal
	}
	if (solver_status==1 && model_status==3)
		return 2; // unbounded
	if (solver_status==1 && (model_status==4 || model_status==5))
		return 1; // infeasible
	if (solver_status==1 && model_status==6)
		return 9; // point cannot be improved (intermediate infeasible)
	if (solver_status==1 && model_status==7)
		return 4; // intermediate nonoptimal
	if (solver_status==2)
		return 3; // iteration limit

	return -1; // other error
}

dvector gamsLocOpt::get_lag_multipliers() {
	if (!gamsptr->reformed) return duals_con;

	dvector lag_mult(prob->con.size());
	for (int c=0; c<gamsptr->objcon; c++) lag_mult[c]=duals_con[c];
	for (int c=gamsptr->objcon; c<prob->con.size(); c++) lag_mult[c]=duals_con[c+1];
	return lag_mult;
}

#endif // COIN_HAS_GAMSIO
