// Copyright (C) Stefan Vigerske 2009
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

// have to be above LaGO because LaGO gets the std namespace, giving a conflict with Couenne's unary_function
#include "LaGOConfig.h"
#ifdef COIN_HAS_COUENNE
#include "exprAbs.hpp"
#include "exprConst.hpp"
#include "exprCos.hpp"
#include "exprDiv.hpp"
#include "exprExp.hpp"
#include "exprInv.hpp"
#include "exprLog.hpp"
#include "exprMax.hpp"
#include "exprMin.hpp"
#include "exprMul.hpp"
#include "exprOpp.hpp"
#include "exprPow.hpp"
#include "exprSin.hpp"
#include "exprSub.hpp"
#include "exprSum.hpp"
#include "exprVar.hpp"
#endif

#include "LaGOG2DFunction.hpp"
#include "LaGOSymSparseMatrix.hpp"

#include <sstream>

#include "g2dexports.h"

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
}


namespace LaGO {

/** The opcodes of GAMS nonlinear expressions.
 */
typedef enum GamsOpCode_ {
	nlNoOp     =  0, /* no operation */
	nlPushV    =  1, /* push variable */
	nlPushI    =  2, /* push immediate (constant) */
	nlStore    =  3, /* store row */
	nlAdd      =  4, /* add */
	nlAddV     =  5, /* add variable */
	nlAddI     =  6, /* add immediate */
	nlAddL     =  7, /* add local */
	nlSub      =  8, /* minus */
	nlSubV     =  9, /* subtract variable */
	nlSubI     = 10, /* subtract immediate */
	nlSubL     = 11, /* subtract local */
	nlMul      = 12, /* multiply */
	nlMulV     = 13, /* multiply variable */
	nlMulI     = 14, /* multiply immediate */
	nlMulL     = 15, /* multiply local */
	nlDiv      = 16, /* divide */
	nlDivV     = 17, /* divide variable */
	nlDivI     = 18, /* divide immediate */
	nlDivL     = 19, /* divide local */
	nlUMin     = 20, /* unary minus */
	nlUMinV    = 21, /* unary minus variable */
	nlSwap     = 22, /* swap two positions on stack top */
	nlPushL    = 23, /* push local */
	nlPopL     = 24, /* pop local */
	nlPopDeriv = 25, /* pop derivative */
	nlHeader   = 26, /* header */
	nlUMinL    = 27, /* push umin local */
	nlStoreS   = 28, /* store scaled row */
	nlPopDerivS= 29, /* store scaled gradient */
	nlEquScale = 30, /* equation scale */
	nlEnd      = 31, /* end of instruction list */
	nlCallArg1 = 32,
	nlCallArg2 = 33,
	nlCallArgN = 34,
	nlFuncArgN = 35,
	nlPushS    = 36,
	nlPopup    = 37,
	nlArg      = 38,
	nlMulIAdd  = 39,
	nlPushZero = 40,
	nlMulPop1  = 41,
	nlMulPop2  = 42,
	nlMulPop   = 43,
	nlAddPop   = 44, 
	nlSubPop   = 45, 
	nlGetConst = 46, 
	nlMulConst1= 47, 
	nlMulConst2= 48, 
	nlMulConst = 49, 
	nlNegLocal = 50, 
	nlGetLocal = 51, 
	nlSetLocal1= 52, 
	nlSetLocal2= 53, 
	nlSetLocal = 54, 
	nlGetGrad  = 55, 
	nlPushIGrad= 56, 
	nlChk      = 57, 
	nlAddO     = 58, 
	nlPushO    = 59,
	nlInvoc    = 60, 
	nlStackIn  = 61,
	MAXINS     = 62
} GamsOpCode;

typedef enum GamsFuncCode_ {fnmapval=0,fnceil,fnfloor,fnround,
    fnmod,fntrunc,fnsign,fnmin,
    fnmax,fnsqr,fnexp,fnlog,
    fnlog10,fnsqrt,fnabs,fncos,
    fnsin,fnarctan,fnerrf,fndunfm,
    fndnorm,fnpower,fnjdate,fnjtime,
    fnjstart,fnjnow,fnerror,fngyear,
    fngmonth,fngday,fngdow,fngleap,
    fnghour,fngminute,fngsecond,
    fncurseed,fntimest,fntimeco,
    fntimeex,fntimecl,fnfrac,fnerrorl,
    fnheaps,fnfact,fnunfmi,fnpi,
    fnncpf,fnncpcm,fnentropy,fnsigmoid,
    fnlog2,fnboolnot,fnbooland,
    fnboolor,fnboolxor,fnboolimp,
    fnbooleqv,fnrelopeq,fnrelopgt,
    fnrelopge,fnreloplt,fnrelople,
    fnrelopne,fnifthen,fnrpower,
    fnedist,fndiv,fndiv0,fnsllog10,
    fnsqlog10,fnslexp,fnsqexp,fnslrec,
    fnsqrec,fncvpower,fnvcpower,
    fncentropy,fngmillisec,fnmaxerror,
    fntimeel,fngamma,fnloggamma,fnbeta,
    fnlogbeta,fngammareg,fnbetareg,
    fnsinh,fncosh,fntanh,fnmathlastrc,
    fnmathlastec,fnmathoval,fnsignpower,
    fnhandle,fnncpvusin,fnncpvupow,
    fnbinomial,fnrehandle,fngamsver,
    fndelhandle,fntan,fnarccos,
    fnarcsin,fnarctan2,fnsleep,fnheapf,
    fncohandle,fngamsrel,fnpoly,
    fnlicensestatus,fnlicenselevel,fnheaplimit,
    fndummy} GamsFuncCode;

const char* GamsOpCodeName[MAXINS] = {
	"nlNoOp",
	"nlPushV",
	"nlPushI",
	"nlStore",
	"nlAdd",
	"nlAddV",
	"nlAddI",
	"nlAddL",
	"nlSub",
	"nlSubV",
	"nlSubI",
	"nlSubL",
	"nlMul",
	"nlMulV",
	"nlMulI",
	"nlMulL",
	"nlDiv",
	"nlDivV",
	"nlDivI",
	"nlDivL",
	"nlUMin",
	"nlUMinV",
	"nlSwap",
	"nlPushL",
	"nlPopL",
	"nlPopDeriv",
	"nlHeader",
	"nlUMinl",
	"nlStoreS",
	"nlPopDerivS",
	"nlEquScale",
	"nlEnd",
	"nlCallArg1",
	"nlCallArg2",
	"nlCallArgN",
	"nlFuncArgN",
	"nlPushS",
	"nlPopup",
	"nlArg",
	"nlMulIAdd",
	"nlPushZero",
	"nlMulPop1",
	"nlMulPop2",
	"nlMulPop",
	"nlAddPop", 
	"nlSubPop", 
	"nlGetConst", 
	"nlMulConst1", 
	"nlMulConst2", 
	"nlMulConst", 
	"nlNegLocal", 
	"nlGetLocal", 
	"nlSetLocal1", 
	"nlSetLocal2", 
	"nlSetLocal", 
	"nlGetGrad", 
	"nlPushIGrad", 
	"nlChk", 
	"nlAddO", 
	"nlPushO",
	"nlInvoc", 
	"nlStackIn"
};

/** Gives the opcode of a GAMS nonlinear instruction.
 */
GamsOpCode getInstrOpCode(unsigned int instr) {
	int iinstr = instr>>26;
/*	assert(iinstr < MAXINS); */
	return (GamsOpCode)iinstr;
}

/** Gives the address in a GAMS nonlinear instruction.
 * The address will be 0-based.
 */
int getInstrAddress(unsigned int instr) {
	return (instr & 67108863)-1;
}

G2DFunction::G2DFunction(int n_instr_, unsigned int* instr_, double* constants_)
: n_instr(n_instr_), constants(constants_), workfactor(1.)
{
	instr = CoinCopyOfArray(instr_, n_instr);	
	
	GamsOpCode opcode;
	set<int> vars;
	for (int i = 0; i < n_instr; ++i)
	{
		opcode = getInstrOpCode(instr[i]);
		switch (opcode)
		{
			case nlPushV:
			case nlAddV:
			case nlSubV:
			case nlMulV:
			case nlDivV:
			case nlUMinV:
				vars.insert(getInstrAddress(instr[i]));
				break;
			default: ;
		}
	}
	
	sparsity.reserve(vars.size());
	for (set<int>::iterator it(vars.begin()); it!=vars.end(); ++it)
		sparsity.push_back(*it);
	
	// get scratch space requirement for directional 2nd derivative evaluation
	int rss = G2DGETRESSTACKSIZE(instr, &n_instr);
	if (rss < n_instr)
		rss = n_instr;
		
	v    = new double[rss];
	vbar = new double[n_instr];
}

G2DFunction::~G2DFunction()
{
	delete[] instr;
	delete[] v;
	delete[] vbar;
}

double G2DFunction::eval(const DenseVector& x) const
{
	double value;
	int numerr = 0;
	double scale;

//	for (int i=0; i<n_instr; ++i)
//		printf("%d ", instr[i]);
//	printf("\n");
		
	G2DFUNCEVAL0(x.getElements(), &value, v, instr, &n_instr, constants, &scale, &numerr);

	if (numerr) {
//		data->domain_violations+=numerr;
		string message("Error in function evaluation.");
		throw FunctionEvaluationError(message, "G2DFunction", "eval");
	}

	return value;
}

void G2DFunction::gradient(DenseVector& grad, const DenseVector& x) const
{
	double val;
	evalAndGradient(val, grad, x);
}

void G2DFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const
{
	int numerr=0;
	double scale;
	
	G2DFUNCEVAL0(x.getElements(), &value, v, instr, &n_instr, constants, &scale, &numerr);

	if (numerr) {
		string message("Error in function evaluation.");
		throw FunctionEvaluationError(message, "G2DFunction", "evalAndGradient");
	}

	G2DREVERSEEVAL1(const_cast<double*>(x.getElements()), grad.getElements(), v, vbar,
		instr, &n_instr, constants, &scale, &numerr);

//	clog << "Grad. of gams con. " << connr << " at " << x << " is " << grad << endl;

	if (numerr) {
		string message("Error in gradient evaluation.");
		throw FunctionEvaluationError(message, "G2DFunction", "evalAndGradient");
	}
}

void G2DFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const
{
	double val, gradvecprod; // to store value and gradient*z
	int numerr = 0;

	G2DDIR2DX(const_cast<double*>(x.getElements()), const_cast<double*>(factor.getElements()), product.getElements(),
			&val, &gradvecprod, v, instr, &n_instr, constants, &numerr);
	
	if (numerr) {
		stringstream message;
		message << "Error " << numerr << " in computation of Hessian-Vector product." << ends;
		throw FunctionEvaluationError(message.str(), "G2DFunction", "hessianVectorProduct");
	}
}

// TODO: implement fullHessian in LaGO::RestrictedFunction
// TODO: improve fullHessian to not reallocate memory every time
void G2DFunction::fullHessian(SymSparseMatrixCreator& hessian, const DenseVector& x) const
{
	double val, scale;
	int numerr = 0;
	hessian.clear();
	hessian.setDim(x.size());
#if 0
	SymSparseMatrixCreator test;
	Function::fullHessian(test, x);
#endif
	// check, if evaluating the function does not produce error
	G2DFUNCEVAL0(x.getElements(), &val, v, instr, &n_instr, constants, &scale, &numerr);
	if (numerr) {
		string message("Error in function evaluation.");
		throw FunctionEvaluationError(message, "G2DFunction", "fullHessian");
	}
	
	int hesHeadPtr, hesTailPtr, hesNumElem;
	int jacHeadPtr, jacTailPtr, jacNumElem;
	int* hesRowCL; int* hesRowRW; int* hesRowNX; double* hesRowVal;
	int* jacVR; int* jacNX; double* jacVal;

	int hesLagSiz = (int)(workfactor * 10 * sparsity.size()); //TODO
	hesRowCL = new int[hesLagSiz];
	hesRowRW = new int[hesLagSiz];
	hesRowNX = new int[hesLagSiz];
	hesRowVal = new double[hesLagSiz];
	int jacSiz = n_instr;
	jacNX = new int[jacSiz];
	jacVR = new int[jacSiz];
	jacVal = new double[jacSiz];
	
	G2DINITJACLIST(jacNX, &jacSiz);
	G2DINITHESLIST(hesRowNX, &hesLagSiz);
	
	G2DHESROWVAL(x.getElements(), &val, instr, constants, &numerr,
		&hesHeadPtr, &hesTailPtr, &hesNumElem,
		&jacHeadPtr, &jacTailPtr, &jacNumElem,
		jacVR, jacNX, jacVal,
		hesRowCL, hesRowRW, hesRowNX, hesRowVal);

	if (numerr) {
		string message("Error in Hessian calculation.");
		throw FunctionEvaluationError(message, "G2DFunction", "fullHessian");
	}

	
	for (int j = hesHeadPtr; j; j = hesRowNX[j-1])
	{ /* hesRowNX is 1-based, so now is also j */
		hessian.insert(hesRowCL[j-1]-1, hesRowRW[j-1]-1, hesRowVal[j-1]);
	}
	
	delete[] hesRowCL;
	delete[] hesRowRW;
	delete[] hesRowNX;
	delete[] hesRowVal;
	delete[] jacNX;
	delete[] jacVR;
	delete[] jacVal;

#if 0
	assert(test==hessian);
#endif
}

#ifdef COIN_HAS_FILIB
interval<double> G2DFunction::eval(const IntervalVector& x) const
{
	double* xmin=new double[x.getNumElements()];
	double* xmax=new double[x.getNumElements()];
	double* xmin_=xmin;
	double* xmax_=xmax;
	const interval<double>* x_=x.getElements();
	for (int i = x.getNumElements(); i > 0; --i) {
		*xmin_ = x_->inf() == filib::fp_traits<double>::ninfinity() ? -1E+20 : x_->inf();
		*xmax_ = x_->sup() == filib::fp_traits<double>::infinity()  ?  1E+20 : x_->sup();
		++xmin_;
		++xmax_;
		++x_;
	}

	double valmin, valmax;
	G2DINTERVAL0X(xmin, xmax, &valmin, &valmax, v, vbar, instr, &n_instr, constants);
	if (valmin <= -1E+20) valmin = filib::fp_traits<double>::ninfinity();
	if (valmax >=  1E+20) valmax = filib::fp_traits<double>::infinity();

	delete[] xmin;
	delete[] xmax;
	
	return interval<double>(valmin, valmax);	
}

void G2DFunction::evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const
{
	assert(x.getNumElements() == grad.getNumElements());
	double* xmin=new double[x.getNumElements()];
	double* xmax=new double[x.getNumElements()];
	double* ymin=new double[x.getNumElements()];
	double* ymax=new double[x.getNumElements()];
	double* xmin_=xmin;
	double* xmax_=xmax;
	double* ymin_=ymin;
	double* ymax_=ymax;
	const interval<double>* x_ = x.getElements();
	for (int i = x.getNumElements(); i > 0; --i) {
		*xmin_=x_->inf() == filib::fp_traits<double>::ninfinity() ? -1E+20 : x_->inf();
		*xmax_=x_->sup() == filib::fp_traits<double>::infinity()  ?  1E+20 : x_->sup();
		*ymin_=0.;
		*ymax_=0.;
		++xmin_;
		++xmax_;
		++ymin_;
		++ymax_;
		++x_;
	}

	double valmin, valmax;
	G2DINTERVAL0X(xmin, xmax, &valmin, &valmax, v, vbar, instr, &n_instr, constants);
	if (valmin <= -1E+20) valmin=filib::fp_traits<double>::ninfinity();
	if (valmax >=  1E+20) valmax=filib::fp_traits<double>::infinity();
	value = interval<double>(valmin, valmax);

	double* v2    = new double[n_instr];
	double* vbar2 = new double[n_instr];

	G2DINTERVAL1X(xmin, xmax, ymin, ymax, v, vbar, v2, vbar2, instr, &n_instr, constants);

	ymin_=ymin;
	ymax_=ymax;
	interval<double>* grad_=grad.getElements();
	for (int i=x.getNumElements(); i>0; --i, ymin_++, ymax_++, ++grad_) {
		if (*ymin_ <= -1E+20) *ymin_=filib::fp_traits<double>::ninfinity();
		if (*ymax_ >=  1E+20) *ymax_=filib::fp_traits<double>::infinity();
		*grad_=interval<double>(*ymin_, *ymax_);
	}

	delete[] xmin;
	delete[] xmax;
	delete[] ymin;
	delete[] ymax;
	delete[] v2;
	delete[] vbar2;	
}
#endif

#ifdef COIN_HAS_COUENNE
/** Reorders instructions such that they do not contain PushS, Popup, or Swap anymore.
 */
void swapInstr(unsigned int* instr, int len1, int len2) {
	int moves = len1+len2;
	int first = 0;
	unsigned int tmp, dest, source;
	
	do {
		++first;
		tmp = instr[first];
		dest = first;
		source = first + len1;
		do {
			instr[dest] = instr[source];
			--moves;
			dest = source;
			if ((int)source <= len2)
				source += len1;
			else
				source -= len2;
		} while((int)source!=first);
		instr[dest] = tmp;
		--moves;
	} while(moves>0);
}

/** Reorders instructions such that they do not contain PushS, Popup, or Swap anymore.
 */
void reorderInstr(unsigned int* instr, int num_instr) {
	int stacklen = 0;
	int j,k;
	int* instrpos = (int*)malloc(sizeof(int)*num_instr);
	instrpos[0] = 0;

	for (k=1; k<num_instr; ++k) {
		GamsOpCode opcode = getInstrOpCode(instr[k]);
		switch (opcode) {
			case nlUMinV:
			case nlPushV:
			case nlPushI:
			case nlPushZero:
				++stacklen;
				break;
			case nlCallArgN:
				stacklen -= getInstrAddress(instr[k-1]); /* getInstrAdress already decreases by 1 */
				break;
			case nlMulIAdd:
			case nlCallArg2:
			case nlStore:
			case nlStoreS:
			case nlAdd:
			case nlSub:
			case nlDiv:
			case nlMul:
				--stacklen;
				break;
			case nlSwap: {
				swapInstr(instr+instrpos[stacklen-2], instrpos[stacklen-1]-instrpos[stacklen-2], instrpos[stacklen]-instrpos[stacklen-1]);
				instrpos[stacklen-1] = instrpos[stacklen-2] - instrpos[stacklen-1] + instrpos[stacklen];
				instr[k] = ((unsigned int)nlArg)<<26; /* nlNoop should be just zero, but confused g2dHesRowStruct, so we take nlArg */
			} break;
			case nlPushS: {
				int len = getInstrAddress(instr[k])+1;
				int pushshift = instrpos[stacklen-len+1] - instrpos[stacklen-len];
/*				for (int i=0; i<=stacklen; ++i) std::clog << i << '\t' << instrpos[i] << std::endl;
				std::clog << "stacklen " << stacklen << " len " << len << std::endl;
				std::clog << "swap " << instrpos[stacklen-len] << ' ' << instrpos[stacklen-len+1] << ' ' << k-2 << std::endl;
				std::clog << "pushshift " << pushshift << std::endl;
*/				swapInstr(instr+instrpos[stacklen-len], pushshift, k-2-instrpos[stacklen-len+1]);
				for (j = stacklen-len+1; j<=stacklen; ++j)
					instrpos[j] -= pushshift;
				--instrpos[stacklen];
				instr[k-1] = ((unsigned int)nlArg)<<26;
				instr[k] = ((unsigned int)nlArg)<<26;
/*				for (int i=0; i<=stacklen; ++i) std::clog << i << '\t' << instrpos[i] << std::endl;
				for (int i=0; i<num_instr; ++i) std::clog << i << '\t' << GamsOpCodeName[getInstrOpCode(instr[i])] << std::endl;
*/				++stacklen;
			} break;
			case nlPopup: {
				stacklen -= getInstrAddress(instr[k])+1;
				instr[k] = ((unsigned int)nlArg)<<26;
			} break;
			default: ;
		}
		instrpos[stacklen] = k;
	}
	
	free(instrpos);
}


expression* G2DFunction::getAsCouenneExpression(std::vector<exprVar*>& vars, Domain* domain) const
{
	unsigned int* instr_ = new unsigned int[n_instr];
	memcpy(instr_, instr, n_instr*sizeof(unsigned int));
	
#define debugoutput 0
	
	// reorder instructions such that there are no PushS, Popup, or Swap left
	reorderInstr(instr_, n_instr);
	
	list<expression*> stack;

	for (int pos = 0; pos < n_instr; ++pos)
	{	
//		std::cout << "stack size: " << stack.size() << std::endl;
		GamsOpCode opcode = getInstrOpCode(instr_[pos]);
		int address = getInstrAddress(instr_[pos]);
		
		expression* exp = NULL;

		if (debugoutput) std::clog << '\t' << GamsOpCodeName[opcode] << ": ";
		switch(opcode) {
			case nlNoOp : { // no operation
				if (debugoutput) std::clog << "ignored" << std::endl;
			} break;
			case nlPushV : { // push variable
//				address = smag->colMapG2S[address];
				if (debugoutput) std::clog << "push variable " << address << std::endl;
				exp = vars[address]->clone(domain);
			} break;
			case nlPushI : { // push constant
				if (debugoutput) std::clog << "push constant " << constants[address] << std::endl;
				exp = new exprConst(constants[address]);
			} break;
			case nlStore: { // store row
				if (debugoutput) std::clog << "ignored" << std::endl;
			} break;
			case nlAdd : { // add
				if (debugoutput) std::clog << "add" << std::endl;
				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = stack.back(); stack.pop_back();
				exp = new exprSum(term1, term2);
			} break;
			case nlAddV: { // add variable
//				address = smag->colMapG2S[address];
				if (debugoutput) std::clog << "add variable " << address << std::endl;
				
				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = vars[address]->clone(domain);
				exp = new exprSum(term1, term2);
			} break;
			case nlAddI: { // add immediate
				if (debugoutput) std::clog << "add constant " << constants[address] << std::endl;
				
				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = new exprConst(constants[address]);
				exp = new exprSum(term1, term2);
			} break;
			case nlSub: { // minus
				if (debugoutput) std::clog << "minus" << std::endl;
				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = stack.back(); stack.pop_back();
				exp = new exprSub(term2, term1);
			} break;
			case nlSubV: { // subtract variable
//				address = smag->colMapG2S[address];
				if (debugoutput) std::clog << "substract variable " << address << std::endl;

				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = vars[address]->clone(domain);
				exp = new exprSub(term1, term2);
			} break;
			case nlSubI: { // subtract immediate
				if (debugoutput) std::clog << "substract constant " << constants[address] << std::endl;
				
				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = new exprConst(constants[address]);
				exp = new exprSub(term1, term2);
			} break;
			case nlMul: { // multiply
				if (debugoutput) std::clog << "multiply" << std::endl;

				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = stack.back(); stack.pop_back();
				exp = new exprMul(term1, term2);
			} break;
			case nlMulV: { // multiply variable
//				address = smag->colMapG2S[address];
				if (debugoutput) std::clog << "multiply variable " << address << std::endl;
				
				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = vars[address]->clone(domain);
				exp = new exprMul(term1, term2);
			} break;
			case nlMulI: { // multiply immediate
				if (debugoutput) std::clog << "multiply constant " << constants[address] << std::endl;

				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = new exprConst(constants[address]);
				exp = new exprMul(term1, term2);
			} break;
			case nlDiv: { // divide
				if (debugoutput) std::clog << "divide" << std::endl;

				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = stack.back(); stack.pop_back();
				if (term2->Type() == CONST)
					exp = new exprMul(term2, new exprInv(term1));
				else
					exp = new exprDiv(term2, term1);
			} break;
			case nlDivV: { // divide variable
//				address = smag->colMapG2S[address];
				if (debugoutput) std::clog << "divide variable " << address << std::endl;
				
				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = vars[address]->clone(domain);
				if (term1->Type() == CONST)
					exp = new exprMul(term1, new exprInv(term2));
				else
					exp = new exprDiv(term1, term2);
			} break;
			case nlDivI: { // divide immediate
				if (debugoutput) std::clog << "divide constant " << constants[address] << std::endl;

				expression* term1 = stack.back(); stack.pop_back();
				expression* term2 = new exprConst(constants[address]);
				exp = new exprDiv(term1, term2);
			} break;
			case nlUMin: { // unary minus
				if (debugoutput) std::clog << "negate" << std::endl;
				
				expression* term = stack.back(); stack.pop_back();
				exp = new exprOpp(term);
			} break;
			case nlUMinV: { // unary minus variable
//				address = smag->colMapG2S[address];
				if (debugoutput) std::clog << "push negated variable " << address << std::endl;

				exp = new exprOpp(vars[address]->clone(domain));
			} break;
			case nlCallArg1 :
			case nlCallArg2 :
			case nlCallArgN : {
				if (debugoutput) std::clog << "call function ";
				GamsFuncCode func = GamsFuncCode(address+1); // here the shift by one was not a good idea
				
				switch (func) {
					case fnmin : {
						if (debugoutput) std::clog << "min" << std::endl;
						
						expression* term1 = stack.back(); stack.pop_back();
						expression* term2 = stack.back(); stack.pop_back();
						exp = new exprMin(term1, term2);
					} break;
					case fnmax : {
						if (debugoutput) std::clog << "max" << std::endl;

						expression* term1 = stack.back(); stack.pop_back();
						expression* term2 = stack.back(); stack.pop_back();
						exp = new exprMax(term1, term2);
					} break;
					case fnsqr : {
						if (debugoutput) std::clog << "square" << std::endl;
						
						expression* term = stack.back(); stack.pop_back();
						exp = new exprPow(term, new exprConst(2.));
					} break;
					case fnexp:
					case fnslexp:
					case fnsqexp: {
						if (debugoutput) std::clog << "exp" << std::endl;

						expression* term = stack.back(); stack.pop_back();
						exp = new exprExp(term);
					} break;
					case fnlog : {
						if (debugoutput) std::clog << "ln" << std::endl;

						expression* term = stack.back(); stack.pop_back();
						exp = new exprLog(term);
					} break;
					case fnlog10:
					case fnsllog10:
					case fnsqlog10: {
						if (debugoutput) std::clog << "log10 = ln * 1/ln(10)" << std::endl;
						
						expression* term = stack.back(); stack.pop_back();
						exp = new exprMul(term, new exprConst(1./log(10.)));
					} break;
					case fnlog2 : {
						if (debugoutput) std::clog << "log2 = ln * 1/ln(2)" << std::endl;

						expression* term = stack.back(); stack.pop_back();
						exp = new exprMul(term, new exprConst(1./log(2.)));
					} break;
					case fnsqrt: {
						if (debugoutput) std::clog << "sqrt" << std::endl;
						
						expression* term = stack.back(); stack.pop_back();
						exp = new exprPow(term, new exprConst(.5));
					} break;
					case fnabs: {
						if (debugoutput) std::clog << "abs" << std::endl;

						expression* term = stack.back(); stack.pop_back();
						exp = new exprAbs(term);
					} break;
					case fncos: {
						if (debugoutput) std::clog << "cos" << std::endl;
						
						expression* term = stack.back(); stack.pop_back();
						exp = new exprCos(term);
					} break;
					case fnsin: {
						if (debugoutput) std::clog << "sin" << std::endl;

						expression* term = stack.back(); stack.pop_back();
						exp = new exprSin(term);
					} break;
					case fnpower:
					case fnrpower: { // x ^ y
						if (debugoutput) std::clog << "power" << std::endl;
						
						expression* term1 = stack.back(); stack.pop_back();
						expression* term2 = stack.back(); stack.pop_back();
						exp = new exprExp(new exprMul(new exprLog(term2), term1));
					} break;
					case fncvpower: { // constant ^ x
						if (debugoutput) std::clog << "power" << std::endl;
						
						expression* term1 = stack.back(); stack.pop_back();
						expression* term2 = stack.back(); stack.pop_back();

						assert(term2->Type() == CONST);
						exp = new exprExp(new exprMul(new exprConst(log(((exprConst*)term2)->Value())), term1));
						delete term2;
					} break;
					case fnvcpower: { // x ^ constant
						if (debugoutput) std::clog << "power" << std::endl;
						
						expression* term1 = stack.back(); stack.pop_back();
						expression* term2 = stack.back(); stack.pop_back();
						assert(term1->Type() == CONST);
						exp = new exprPow(term2, term1);
					} break;
					case fnpi: {
						if (debugoutput) std::clog << "pi" << std::endl;
						//TODO
						assert(false);
					} break;
					case fndiv:
					case fndiv0: {
						expression* term1 = stack.back(); stack.pop_back();
						expression* term2 = stack.back(); stack.pop_back();
						if (term2->Type() == CONST)
							exp = new exprMul(term2, new exprInv(term1));
						else
							exp = new exprDiv(term2, term1);
					} break;
					case fnslrec: // 1/x
					case fnsqrec: { // 1/x
						if (debugoutput) std::clog << "divide" << std::endl;
						
						expression* term = stack.back(); stack.pop_back();
						exp = new exprInv(term);
					} break;
					case fnceil: case fnfloor: case fnround:
					case fnmod: case fntrunc: case fnsign:
					case fnarctan: case fnerrf: case fndunfm:
					case fndnorm: case fnerror: case fnfrac: case fnerrorl:
			    case fnfact /* factorial */: 
			    case fnunfmi /* uniform random number */:
			    case fnncpf /* fischer: sqrt(x1^2+x2^2+2*x3) */:
			    case fnncpcm /* chen-mangasarian: x1-x3*ln(1+exp((x1-x2)/x3))*/:
			    case fnentropy /* x*ln(x) */: case fnsigmoid /* 1/(1+exp(-x)) */:
			    case fnboolnot: case fnbooland:
			    case fnboolor: case fnboolxor: case fnboolimp:
			    case fnbooleqv: case fnrelopeq: case fnrelopgt:
			    case fnrelopge: case fnreloplt: case fnrelople:
			    case fnrelopne: case fnifthen:
			    case fnedist /* euclidian distance */:
			    case fncentropy /* x*ln((x+d)/(y+d))*/:
			    case fngamma: case fnloggamma: case fnbeta:
			    case fnlogbeta: case fngammareg: case fnbetareg:
			    case fnsinh: case fncosh: case fntanh:
			    case fnsignpower /* sign(x)*abs(x)^c */:
			    case fnncpvusin /* veelken-ulbrich */:
			    case fnncpvupow /* veelken-ulbrich */:
			    case fnbinomial:
			    case fntan: case fnarccos:
			    case fnarcsin: case fnarctan2 /* arctan(x2/x1) */:
			    case fnpoly: /* simple polynomial */
					default : {
						if (debugoutput) std::cerr << "nr. " << (int)func << " - unsuppored. Error." << std::endl;
						return NULL;
					}				
				}
			} break;
			case nlMulIAdd: {
				if (debugoutput) std::clog << "multiply constant " << constants[address] << " and add " << std::endl;

				expression* term1 = stack.back(); stack.pop_back();
				term1 = new exprMul(term1, new exprConst(constants[address]));
				expression* term2 = stack.back(); stack.pop_back();
				
				exp = new exprSum(term1, term2);				
			} break;
			case nlFuncArgN : {
				if (debugoutput) std::clog << "ignored" << std::endl;
			} break;
			case nlArg: {
				if (debugoutput) std::clog << "ignored" << std::endl;
			} break;
			case nlHeader: { // header
				if (debugoutput) std::clog << "ignored" << std::endl;
			} break;
			case nlPushZero: {
				if (debugoutput) std::clog << "push constant zero" << std::endl;
				exp = new exprConst(0.);
			} break;
			case nlStoreS: { // store scaled row
				if (debugoutput) std::clog << "ignored" << std::endl;
			} break;
			// the following three should have been taken out by reorderInstr above; the remaining ones seem to be unused by now
			case nlPushS: // duplicate value from address levels down on top of stack
			case nlPopup: // duplicate value from this level to at address levels down and pop entries in between
			case nlSwap: // swap two positions on top of stack
			case nlAddL: // add local
			case nlSubL: // subtract local
			case nlMulL: // multiply local
			case nlDivL: // divide local
			case nlPushL: // push local
			case nlPopL: // pop local
			case nlPopDeriv: // pop derivative
			case nlUMinL: // push umin local
			case nlPopDerivS: // store scaled gradient
			case nlEquScale: // equation scale
			case nlEnd: // end of instruction list
			default: {
				std::cerr << "not supported - Error." << std::endl;
				return NULL;
			}
		}
		
		if (exp)
			stack.push_back(exp);
	}

	delete[] instr_;

	assert(stack.size() == 1);	
	return stack.back();
}
#endif

void G2DFunction::print(ostream& out) const
{
	out << "G2DFunction on variables";
	for (size_t i = 0; i < sparsity.size(); ++i)
		out << ' ' << sparsity[i];
	out << endl;
}

} /* namespace minlp */
