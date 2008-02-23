// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOGamsReader.hpp"
#include "LaGOGamsFunction.hpp"
#include "LaGOMINLPData.hpp"
#include "LaGOScaledFunction.hpp"

//#include <dlfcn.h>

namespace LaGO {

#define CNAMES
#ifdef NOUNDERSCORE
#define FNAME_LCASE_NODECOR
#else
#define FNAME_LCASE_DECOR
#endif

#include "iolib.h"
#include "dict.h"
#include "nliolib.h"
//#include "gcprocs.h"
#include "g2dexports.h"
#include "clicelib.h"

extern "C" {
#if defined(COIN_HAS_CPX) && defined(COIN_HAS_GAMSCPLEXLICE)
#include "gamscplexlice.h"
#endif

//GDX client
//#ifdef GDX_AVAILABLE
//#include "gdxwrap.h"
//#endif
}

GamsReader::GamsReader()
: obj_factor(1), is_minimization(true), reformed(true), con_types(NULL), con_rhs(NULL)
{ data=new Data();
//#ifdef GDX_AVAILABLE
//	char* errormsg=new char[128];
//	if (gdxWrapInit(NULL, errormsg, 128)) {
//	 	cerr << "Could not load GDX I/O library: " << errormsg << endl;
//		exit(-1);
//	}
//	delete[] errormsg;
//#endif
}

GamsReader::~GamsReader() {
	gfclos();
	delete[] con_types;
	delete[] con_rhs;
}

char* GamsReader::getRowName(struct dictRec* dict, int gi, char *name, int bufLen) {
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

char* GamsReader::getColName (struct dictRec* dict, int gj, char *name, int bufLen) {
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

#if defined(COIN_HAS_CPX) && defined(COIN_HAS_GAMSCPLEXLICE)
void GamsReader::initCPLEXLicense(int connr, int varnr, int nnz, int nlnz, int ndisc) const {
	licenseInit_t initType;
	if (gamscplexlice(connr, varnr, nnz, nlnz, ndisc, 1, &initType, NULL, NULL, NULL, NULL, NULL, NULL)) {
		cerr << "Could not initialize CPLEX license" << endl;
  } else cout << "GAMS/CPLEX license accepted: " << initType << endl;
}
#endif

SmartPtr<MINLPData> GamsReader::getProblem(char* cntr_file) {
	assert(cntr_file);

	SmartPtr<MINLPData> prob(new MINLPData());

	gfinit();
	
	clog << "Reading control file " << cntr_file << endl;
	cntrec info;
	gfrcnt(true, /* because we can handle discrete variables */
		true, /* because we can handle nonlinearities */
		&info,
		cntr_file
	);
	iolib.instrStartIdx=1;
	
	// open status file
	gfopst();

#if defined(COIN_HAS_CPX) && defined(COIN_HAS_GAMSCPLEXLICE)
	initCPLEXLicense(iolib.nrows, iolib.ncols, iolib.nnz, iolib.nlnz, iolib.ndisc);
#endif

	int numrow=info.kgv[1];
	int numcol=info.kgv[2];
	if (info.lgv[3]) obj_factor=1;
	else {
		obj_factor=-1;
		is_minimization=false;
	}
	objcon=info.kgv[15]-1;
	objvar=info.kgv[7];
	
	cout << "Rows: " << numrow << endl;
	cout << "Columns: " << numcol << endl;
  
	reformed&=(info.kgv[14]==1); // reformable, when objective var appears linear in objective
	prob->var.reserve(numcol);
	prob->con.reserve(numrow);

	// reading row informations
	vector<double> conlhs(numrow, -getInfinity()); // lower bounds of constraints
	vector<double> conrhs(numrow,  getInfinity()); // upper bounds of constraints
	con_types=new int[numrow];
	con_rhs=new double[numrow];
	rowRec_t rowdata;
	for (int c=0; c<numrow; c++) {
		cioReadRow(&rowdata);
		con_rhs[c]=rowdata.rdata[2];
		con_types[c]=rowdata.idata[1];
		switch (rowdata.idata[1]) {
			case 0: // equality
				conlhs[c]=rowdata.rdata[2];
				conrhs[c]=rowdata.rdata[2];
				break;
			case 1: // greater-equal
				conlhs[c]=rowdata.rdata[2];
				break;
			case 2: // lower-equal
				conrhs[c]=rowdata.rdata[2];
				break;
			case 3: // free row
			default:
				break;  
		}
	}

	reformed&=(conlhs[objcon]==conrhs[objcon]); // objective-row is an equality

	// reading dictionary file
	if (iolib.dictFileWritten) {
    int ret = gcdLoad(&data->dict, iolib.flndic, iolib.dictVersion);
		if (ret) {
			cerr << "Error reading dictionary file." << endl;
			data->dict=NULL;
		}
	}

	// coefficients of linear part of constraints
	vector<SparseVectorCreator> conLin(numrow);
	// nonzero variables in constraints
	vector<set<int> > conSparsity(numrow);

	prob->start_points.resize(1);
	DenseVector& start(prob->start_points.front());
	start.resize(numcol);
	
	// reading columns
	colRec_t coldata;
	char* namebuf=new char[50];
	char* name=NULL;
	double low, up;
	for (int i=0; i<numcol; i++) {
		cioReadCol(&coldata);
		if (coldata.idata[3]>2) {
			cerr << "Variable " << i << " has type " << coldata.idata[3] << endl;
			cerr << "This kind of variables is not supported, aborting." << endl;
			exit(EXIT_FAILURE);
		}
		if (coldata.cdata[1]!=info.xgv[11]) low=coldata.cdata[1];
		else low=-getInfinity();
		if (coldata.cdata[3]!=info.xgv[10]) up=coldata.cdata[3];
		else up=getInfinity();
		if (i==objvar) // objective variable; substituted later
			reformed&=(low==-getInfinity()) && (up==getInfinity()); // objective variables is free

		if (data->dict && !(name=getColName(data->dict, i, namebuf, 50)))
			cerr << "Couldn't retrieve name for variable " << i << endl;
		start[i]=coldata.cdata[2];

		double coeff;
		int index;
		int nltyp;
		bool isnonlinear=false; // whether the variable appears nonlinear in one of the rows
		for (int r=0; r<coldata.idata[1]; r++) {
			gfrcof(&coeff, &index, &nltyp);
			if (nltyp==0) // linear variable
				conLin[index-1].insert(i, coeff);
			else {
				isnonlinear=true;
				conSparsity[index-1].insert(i);
			}
			if (index-1==objcon && i==objvar) reformed&=(nltyp==0); // objective variable appears linear
		}

		prob->var.push_back(MINLPData::Variable(i, low, up, coldata.idata[3], isnonlinear, name));
		if (coldata.idata[3]) prob->discrete_var.push_back(i);  
	}
	cout << "Discrete variables: " << prob->discrete_var.size() << endl;

	cout << "Reformable: " << (reformed ? "yes" : "no") << endl;

	if (reformed) {
		prob->var[objvar].lower=prob->var[objvar].upper=start[objvar]=0.;
		if (prob->numDiscrVariables()==prob->numVariables()-1 && (!prob->var[objvar].isDiscrete())) {
			prob->var[objvar].discrete=true;
			prob->discrete_var.push_back(objvar); // TODO: insert such that vector of indices stays sorted?
		}
	}

	char* objcon_name;
	for (int c=0; c<numrow; c++) {
		if (data->dict && !(name=getRowName(data->dict, c, namebuf, 50)))
			cerr << "Couldn't retrieve name for row " << c << endl;
		if ((c==objcon) && name) objcon_name=strdup(name);
		if ((c==objcon) && reformed) { // reforming
			double& objvar_coeff(conLin[c][objvar]);
			obj_factor*=-objvar_coeff;
			objvar_coeff=0.;
			prob->obj=MINLPData::Objective(NULL, new SparseVector(conLin[c]), -conrhs[c], name);
			
			if (obj_factor!=1) {
				assert(obj_factor!=0);
				*prob->obj.origfuncLin/=obj_factor;
				prob->obj.origfuncConstant/=obj_factor;
			}
		} else
			prob->con.push_back(MINLPData::Constraint(prob->con.size(), conlhs[c], conrhs[c], NULL, new SparseVector(conLin[c]), 0., name));
	}
	if (!reformed) {
		prob->obj=MINLPData::Objective(NULL, new SparseVector(objvar, obj_factor), 0, objcon_name);
		if (objcon_name) free(objcon_name);
	}

	delete[] namebuf;

	G2DINIT();
	G2DSTATSINIT(gci2DData);
//	G2DSETMSGCALLBACK(MSGCB);

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
			data->maxins=CoinMax(data->maxins, data->numInstr[i]);
			int rss=G2DGETRESSTACKSIZE(data->instr+data->startIdx[i]-1, &(data->numInstr[i]));
			data->resstacksize=CoinMax(data->resstacksize, rss);
			if ((i==objcon) && reformed) {
				if (obj_factor==1) prob->obj.origfuncNL=new GamsFunction(i, data, conSparsity[i]);
				else prob->obj.origfuncNL=new ScaledFunction(new GamsFunction(i, data, conSparsity[i]), 1/obj_factor);
			} else {
				int index=i-((i>objcon && reformed) ? 1 : 0);
				prob->con[index].origfuncNL=new GamsFunction(i, data, conSparsity[i]);
			}
			++nr_nonlincon;
		} else {
			data->startIdx[i]=0;
			data->numInstr[i]=0;
		}
	cout << "Nonlinear constraints: " << nr_nonlincon << endl;

	data->resstack=new double[data->resstacksize];
	data->s=new double[data->lenins];
	data->sbar=new double[data->lenins];

	delete[] instr;
	delete[] nlwork;
	delete[] startIdx;
	delete[] numInstr;

	return prob;
}

string GamsReader::getParameterfilename() {
	if (iolib.useopt) return iolib.flnopt;
	return "";
}


void GamsReader::writeSolutionFile(const DenseVector& sol_point, double obj_value, int model_status, int solver_status, int iter, double time) {
	if (model_status>2) { // not solved
		gfwsolflag(0); // not writing a solution
		gfwsta(model_status, solver_status, iter, time, 0, 0);
		return;
	}

	if (!is_minimization) obj_value*=-1.;
	gfwsta(model_status, solver_status, iter, time, obj_value, 0); 

	for (int c=0; c<iolib.nrows; c++)
		gfwrow(0., 0., 3, 0); // TODO: give at least correct row values

	for (int i=0; i<iolib.ncols; i++) {
//		if (prob->discr[i]) solver->basind_var[i]=3; // superbasic for discrete variables
		if (i==objvar && reformed)
			gfwcol(obj_value, 0., 3, 0);
		else
			gfwcol(sol_point(i), 0., 3, 0);
	}
}

GamsReader::Data::Data()
: dict(NULL), instr(NULL), nlCons(NULL), startIdx(NULL), numInstr(NULL),
	lenins(1), maxins(1),
	s(NULL), sbar(NULL), resstack(NULL), resstacksize(1),
	domain_violations(0)
{ }

GamsReader::Data::~Data() {
	delete[] instr;
	delete[] nlCons;
	delete[] startIdx;
	delete[] numInstr;
	delete[] resstack;
	delete[] s;
	delete[] sbar;
	if (dict) gcdFree(dict);
}

} // namespace LaGO
