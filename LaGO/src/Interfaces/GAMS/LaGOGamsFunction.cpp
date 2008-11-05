
#include "LaGOGamsFunction.hpp"

namespace LaGO {

#define CNAMES
#ifdef NOUNDERSCORE
#define FNAME_LCASE_NODECOR
#else
#define FNAME_LCASE_DECOR
#endif

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

GamsFunction::GamsFunction(int connr_, const SmartPtr<GamsReader::Data>& data_, const set<int>& sparsity_)
: data(data_), connr(connr_)
{ sparsity.reserve(sparsity_.size());
	for (set<int>::iterator it(sparsity_.begin()); it!=sparsity_.end(); ++it)
		sparsity.push_back(*it);
}

void GamsFunction::addRowName(string& message) const {
	char namebuf[50];
	char* name=NULL;
	if (data->dict)
		name=GamsReader::getRowName(data->dict, connr, namebuf, 50);
	if (!name) {
		snprintf(namebuf,50,"%d",connr);
		name=namebuf;
	}
	message+=name;
}

double GamsFunction::eval(const DenseVector& x) const {
	double value;
	int numerr=0;
	double dt;

	G2DFUNCEVAL0(x.getElements(), &value, data->s,
	data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons,
	&dt, &numerr);

//	clog << "Value of gams con. " << connr << " at " << x << " is " << value << endl;
	
	if (numerr) {
		data->domain_violations+=numerr;
		string message("Error evaluating constraint ");
		addRowName(message);
		throw FunctionEvaluationError(message, "GamsFunction", "eval");
	}

	return value;
}

void GamsFunction::evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const {
	int numerr=0;
	double dt;
	
	G2DFUNCEVAL0(x.getElements(), &value, data->s,
		data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons,
		&dt, &numerr);

	if (numerr) {
		data->domain_violations+=numerr;
		string message("Error evaluating constraint ");
		addRowName(message);
		throw FunctionEvaluationError(message, "GamsFunction", "evalAndGradient");
	}

	G2DREVERSEEVAL1(const_cast<double*>(x.getElements()), grad.getElements(), data->s, data->sbar,
		data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons,
		&dt, &numerr);

//	clog << "Grad. of gams con. " << connr << " at " << x << " is " << grad << endl;

	if (numerr) {
		data->domain_violations+=numerr;
		string message("Error computing gradient of constraint ");
		addRowName(message);
		throw FunctionEvaluationError(message, "GamsFunction", "evalAndGradient");
	}
}

void GamsFunction::gradient(DenseVector& grad, const DenseVector& x) const {
	double value;
	evalAndGradient(value, grad, x);
}

void GamsFunction::hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const {
	double val, gradvecprod; // to store value and gradient*z
	int numerr=0;

	G2DDIR2DX(const_cast<double*>(x.getElements()), const_cast<double*>(factor.getElements()), product.getElements(), &val, &gradvecprod,
		data->resstack, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons, &numerr);
	
	if (numerr) {
		data->domain_violations+=numerr;
		string message("Error ");
		char buf[5];
		snprintf(buf,5,"%d",numerr);
		message+=buf;
		message+=" computing Hessian-Vector product of constraint ";
		addRowName(message);
		throw FunctionEvaluationError(message, "GamsFunction", "hessianVectorProduct");
	}
}

void GamsFunction::print(ostream& out) const {
	char namebuf[50];
	if (data->dict && GamsReader::getRowName(data->dict, connr, namebuf, 50))
		out << namebuf;
	else
		out << "GamsFunction " << connr;	
}	

#ifdef COIN_HAS_FILIB
interval<double> GamsFunction::eval(const IntervalVector& x) const {
	double* xmin=new double[x.getNumElements()];
	double* xmax=new double[x.getNumElements()];
	double* xmin_=xmin;
	double* xmax_=xmax;
	const interval<double>* x_=x.getElements();
	for (int i=x.getNumElements(); i>0; --i) {
		*xmin_=x_->inf();
		*xmax_=x_->sup();
		++xmin_;
		++xmax_;
		++x_;
	}

	double valmin, valmax;
	G2DINTERVAL0X(xmin, xmax, &valmin, &valmax, data->s, data->sbar, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons);
	if (valmin==-1E+20) valmin=filib::fp_traits<double>::ninfinity();
	if (valmax==1E+20)  valmax=filib::fp_traits<double>::infinity();

	delete[] xmin;
	delete[] xmax;
	
	return interval<double>(valmin, valmax);
}

void GamsFunction::evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const {
	assert(x.getNumElements()==grad.getNumElements());
	double* xmin=new double[x.getNumElements()];
	double* xmax=new double[x.getNumElements()];
	double* ymin=new double[x.getNumElements()];
	double* ymax=new double[x.getNumElements()];
	double* xmin_=xmin;
	double* xmax_=xmax;
	double* ymin_=ymin;
	double* ymax_=ymax;
	const interval<double>* x_=x.getElements();
	for (int i=x.getNumElements(); i>0; --i) {
		*xmin_=x_->inf();
		*xmax_=x_->sup();
		*ymin_=0.;
		*ymax_=0.;
		++xmin_;
		++xmax_;
		++ymin_;
		++ymax_;
		++x_;
	}

	double valmin, valmax;
	G2DINTERVAL0X(xmin, xmax, &valmin, &valmax, data->s, data->sbar, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons);
	if (valmin==-1E+20) valmin=filib::fp_traits<double>::ninfinity();
	if (valmax==1E+20)  valmax=filib::fp_traits<double>::infinity();
	value=interval<double>(valmin, valmax);

	double* s2=new double[data->lenins];
	double* sbar2=new double[data->lenins];

	G2DINTERVAL1X(xmin, xmax, ymin, ymax, data->s, data->sbar, s2, sbar2, data->instr+data->startIdx[connr]-1, &(data->numInstr[connr]), data->nlCons);

	ymin_=ymin;
	ymax_=ymax;
	interval<double>* grad_=grad.getElements();
	for (int i=x.getNumElements(); i>0; --i, ++ymin_, ++ymax_, ++grad_) {
		if (*ymin_==-1E+20) *ymin_=filib::fp_traits<double>::ninfinity();
		if (*ymax_==1E+20)  *ymax_=filib::fp_traits<double>::infinity();
		*grad_=interval<double>(*ymin_, *ymax_);
	}

	delete[] xmin;
	delete[] xmax;
	delete[] ymin;
	delete[] ymax;
	delete[] s2;
	delete[] sbar2;
}
#endif
	
} // namespace LaGO
