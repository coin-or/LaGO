// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOQuadraticEstimation.hpp"
#include "LaGOSampling.hpp"
#include "LaGOBoxMinimizationProblem.hpp"
#include "LaGOSumFunction.hpp"
#include "LaGOScaledFunction.hpp"

namespace LaGO {
	
QuadraticEstimation::QuadraticEstimation()
: eps(1E-5), iter_max(100)
{
//	ipopt.Options()->SetStringValue("hessian_approximation", "limited-memory");
	ipopt.Options()->SetNumericValue("tol", eps);
	ipopt.Options()->SetNumericValue("dual_inf_tol", eps);
	ipopt.Options()->SetIntegerValue("print_level", 0);
	//	ipopt.Initialize(); // this reads ipopt.opt
	ipopt.Initialize("");
}

int QuadraticEstimation::computeEstimators(MINLPData& data) {
	int nr=computeEstimators(data, data.obj, true, false); // compute underestimators of objective
	for (int c=0; c<data.numConstraints(); ++c)
		nr+=computeEstimators(data, data.con[c], data.con[c].upper<getInfinity(), data.con[c].lower>-getInfinity());
	return nr;
}

int QuadraticEstimation::computeImprovingEstimators(MINLPData& data, const DenseVector& x) {
	int nr=computeImprovingEstimators(data, data.obj, x, true, false); // compute underestimators of objective
	for (int c=0; c<data.numConstraints(); ++c)
		nr+=computeImprovingEstimators(data, data.con[c], x, data.con[c].upper<getInfinity(), data.con[c].lower>-getInfinity());
	return nr;
}

int QuadraticEstimation::computeEstimators(MINLPData& data, MINLPData::ObjCon& objcon, bool need_lower, bool need_upper) {
	int count=0;
	for (unsigned int k=0; k<objcon.decompfuncNL.size(); ++k) {
		BlockFunction& func(*objcon.decompfuncNL[k]);
		if (IsNull(func.nonquad)) continue;

		// check which under/overestimator we really need to compute
		bool do_lower=need_lower && !(func.curvature&CONVEX);
		bool do_upper=need_upper && !(func.curvature&CONCAVE);
		if ((!do_lower) && (!do_upper)) continue;

		clog << "Quadratic Estimators " << objcon.name << " block " << k << ':' << ' ';

		lp.reset();

		DenseVector lower, upper;
		data.getBox(lower, upper, func.indices);

		BlockFunctionWrapper funcwrap(func);

		pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > estimators=computeFirstEstimator(funcwrap, lower, upper, do_lower, do_upper);
		if (do_lower && IsValid(estimators.first)) {
			func.underestimators.push_back(new QuadraticEstimator(estimators.first));
			++count;
		}
		if (do_upper && IsValid(estimators.second)) {
			func.overestimators.push_back(new QuadraticEstimator(estimators.second));
			++count;
		}
		
		clog << endl;
	}
	return count;
}

int QuadraticEstimation::computeImprovingEstimators(MINLPData& data, MINLPData::ObjCon& objcon, const DenseVector& refpoint, bool need_lower, bool need_upper) {
	int count=0;
	for (unsigned int k=0; k<objcon.decompfuncNL.size(); ++k) {
		BlockFunction& func(*objcon.decompfuncNL[k]);
		if (IsNull(func.nonquad)) continue;
		
		DenseVector refpoint_block(refpoint, func.indices);
		double origval=func.nonquad->eval(refpoint_block);
		
		bool do_lower=need_lower && !(func.curvature&CONVEX);
		bool do_upper=need_upper && !(func.curvature&CONCAVE);
		if ((!do_lower) && (!do_upper)) continue;

		clog << "Quadratic Estimators " << objcon.name << " block " << k << ": orig. val: " << origval << '\t';
		
		// check which under/overestimator we want to compute
		if (do_lower) {
			double lowerest=func.evalUnderEstimator(refpoint_block);
			clog << "lower est.: " << lowerest << '\t';
			// skip new underestimator if gap is small
			if ((origval-lowerest)<.1*CoinAbs(origval)) do_lower=false;
		}
		
		if (do_upper) {
			double upperest=func.evalOverEstimator(refpoint_block);			
			clog << "upper est.: " << upperest << '\t';
			// skip new overestimator if gap is small
			if ((upperest-origval)<.1*CoinAbs(origval)) do_upper=false;
		}
		
		if ((!do_lower) && (!do_upper)) {
			clog << endl;
			continue;
		}

		// add reference point to sample set
		SampleSet::iterator refpoint_it=func.samplepoints.insert(refpoint_block).first;
		refpoint_it->funcvalue=origval;
		
		lp.reset();

		DenseVector lower, upper;
		data.getBox(lower, upper, func.indices);

		BlockFunctionWrapper funcwrap(func);

		pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > estimators=computeAdditionalEstimator(funcwrap, lower, upper, refpoint_it, do_lower, do_upper);
		if (do_lower && IsValid(estimators.first)) {
			func.underestimators.push_back(new QuadraticEstimator(estimators.first));
			double lowerest=func.evalUnderEstimator(refpoint_block);
			clog << "new lower est.: " << lowerest << '\t';
			++count;
		}
		if (do_upper && IsValid(estimators.second)) {
			func.overestimators.push_back(new QuadraticEstimator(estimators.second));
//			double upperest=func.evalOverEstimator(refpoint_block);
			double upperest=estimators.second->eval(refpoint_block);
			clog << "new upper est.: " << upperest << '\t';
			++count;
		}
		
		clog << endl;
	}
	
	return count;
}


pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > QuadraticEstimation::computeFirstEstimator(NonconvexFunction& func, const DenseVector& lower, const DenseVector& upper, bool do_lower, bool do_upper) {
	pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > estimators;
	if ((!do_lower) && (!do_upper)) return estimators; // nothing todo

	// set sample point function values
	SampleSet::iterator it_sp(func.getSamplePoints().begin());
	while (it_sp!=func.getSamplePoints().end()) {
		try {
			if (it_sp->funcvalue==getInfinity())
				it_sp->funcvalue=func.getFunction()->eval(*it_sp);
		} catch (FunctionEvaluationError error) {
			cerr << "QuadraticEstimation::computeEstimator: Error evaluation function in sample point, skipping point. " << error << endl;
			it_sp=func.getSamplePoints().eraseAndGetNext(it_sp);
			continue;
		}
		++it_sp;
	}
	
	Sampling sampling;

	if (do_lower) {
		clog << "lower: ";
		// compute a local minimizer of our function and add it to the sample set
		SampleSet::iterator enforce_tightness=sampling.addMinimizer(func.getSamplePoints(), *func.getFunction(), lower, upper, func.getSparsityGraph(), &func.getSamplePoints().begin()->getPoint());
		// if we couldn't find one, search for a sample point that is a startingpoint of the MINLP
		if (enforce_tightness==func.getSamplePoints().end()) {
			enforce_tightness=func.getSamplePoints().begin();
			while (enforce_tightness!=func.getSamplePoints().end() && !enforce_tightness->is_startpoint)
				++enforce_tightness;
			// if we also don't have such one, pick just the first sample point
			if (enforce_tightness==func.getSamplePoints().end())
				enforce_tightness=func.getSamplePoints().begin();
		}

		// initialize the LP which we use to generate our quad. underestimator
		int enforce_tightness_index=initLP(func, enforce_tightness);

		// and set the right-hand-side such that we get an underestimator 
		updateLP(true);

		estimators.first=getEstimator(func, true, lower, upper);

		if (do_upper) { // release tightness in reference sample point
			lp.setColBounds(enforce_tightness_index, 0, getInfinity());
		}
	}

	if (do_upper) {
		clog << "upper: ";
		// generate maximizer
		// compute a local maximizer of our function and add it to the sample set
		ScaledFunction minusf(func.getFunction(), -1.);
		SampleSet::iterator enforce_tightness=sampling.addMinimizer(func.getSamplePoints(), minusf, lower, upper, func.getSparsityGraph(), &func.getSamplePoints().begin()->getPoint());
		if (enforce_tightness==func.getSamplePoints().end()) {
			enforce_tightness=func.getSamplePoints().begin();
			while (enforce_tightness!=func.getSamplePoints().end() && !enforce_tightness->is_startpoint)
				++enforce_tightness;
			// if we also don't have such one, pick just the first sample point
			if (enforce_tightness==func.getSamplePoints().end())
				enforce_tightness=func.getSamplePoints().begin();
		} else { // the maximizer of f was computed by minimizing -f
			enforce_tightness->funcvalue*=-1;
		}
		if (!do_lower) // if we didn't compute a lower underestimator, we need to initialize the LP
			initLP(func, enforce_tightness);
		else { // otherwise we add only the row corresponding to the new sample point 
			double scale;
			SparseVector* row=constructRow(func, *enforce_tightness, -1, scale);

			double rhs=enforce_tightness->funcvalue/scale;
			double lhs=rhs-eps*fabs(rhs);

			lp.addRow(*row, lhs, rhs);
			delete row;

			sampleset_newsort.push_front(SampleSetItem(*enforce_tightness, -1, lp.getNumRows()-1, scale));
		}
		// set the right-hand-side such that we get an overestimator
		updateLP(false);

		estimators.second=getEstimator(func, false, lower, upper);
	}

	return estimators;
}


pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > QuadraticEstimation::computeAdditionalEstimator(NonconvexFunction& func, const DenseVector& lower, const DenseVector& upper, SampleSet::iterator enforce_tightness, bool do_lower, bool do_upper) {
	pair<SmartPtr<QuadraticFunction>, SmartPtr<QuadraticFunction> > estimators;
	if ((!do_lower) && (!do_upper)) return estimators; // nothing todo

	if (do_lower) {
		clog << "lower: ";
		// compute a local minimizer of our function and add it to the sample set

		// initialize the LP which we use to generate our quad. underestimator
		int enforce_tightness_index=initLP(func, enforce_tightness);

		// and set the right-hand-side such that we get an underestimator 
		updateLP(true);

		estimators.first=getEstimator(func, true, lower, upper);

		if (do_upper) { // release tightness in reference sample point
			lp.setColBounds(enforce_tightness_index, 0, getInfinity());
		}
	}

	if (do_upper) {
		clog << "upper: ";
		// generate maximizer
		if (!do_lower) // if we didn't compute a lower underestimator, we need to initialize the LP
			initLP(func, enforce_tightness);
		else { // otherwise we add only the row corresponding to the new sample point 
//			double scale;
//			SparseVector* row=constructRow(func, *enforce_tightness, -1, scale);
//
//			double rhs=enforce_tightness->funcvalue/scale;
//			double lhs=rhs-eps*fabs(rhs);
//
//			lp.addRow(*row, lhs, rhs);
//			delete row;
//
//			sampleset_newsort.push_front(SampleSetItem(*enforce_tightness, -1, lp.getNumRows()-1, scale));
		}
		// set the right-hand-side such that we get an overestimator
		updateLP(false);

		estimators.second=getEstimator(func, false, lower, upper);
	}

	return estimators;
}



int QuadraticEstimation::initLP(NonconvexFunction& func, const SampleSet::iterator& enforce_tightness) {
// create auxiliary LP with initial sample set
	int hess_entries=0;
	for (SparsityGraph::arc_iterator it_arc(func.getSparsityGraph()->arc_begin()); it_arc!=func.getSparsityGraph()->arc_end(); ++it_arc) {
		const SparsityGraph::Arc& arc(*it_arc);
		int var1=(**arc.head()).varindex;
		int var2=(**arc.tail()).varindex;
		if (var1<=var2) ++hess_entries;
	}
	nr_coeff=1+func.dim()+hess_entries;
	nr_auxvars=func.getSamplePoints().size();
	int nr_cols=nr_coeff+nr_auxvars;

	SparseVector zero;
	CoinPackedVectorBase** cols=new CoinPackedVectorBase*[nr_cols];
	for (int i=0; i<nr_cols; ++i) cols[i]=&zero;

	DenseVector collb(nr_cols, 0.);
	DenseVector colub(nr_cols, getInfinity());
//	collb[0]=-getInfinity();
	for (int i=0; i<nr_coeff; ++i)
		collb[i]=-getInfinity();

//	// setting up some bounds on the coefficients:
//	// let coefficients in linear and quadratic part not exceed 10 times the absolute value of corresponding coefficient in reference point

//	DenseVector grad_in_refpoint(func.indices.size());
//	func.nonquad->gradient(grad_in_refpoint, *enforce_tightness);
//	int index=1;
//	for (; index<1+(int)func.indices.size(); ++index) {
//		double bound=CoinMax(100., 10*CoinAbs(grad_in_refpoint(index-1)));
//		collb[index]=-bound;
//		colub[index]=bound;
//	}
//
//	SymSparseMatrixCreator hessian_in_refpoint(func.indices.size());
//	func.nonquad->fullHessian(hessian_in_refpoint, *enforce_tightness);
//	for (SparsityGraph::arc_iterator it_arc(func.sparsitygraph->arc_begin()); it_arc!=func.sparsitygraph->arc_end(); ++it_arc) {
//		const SparsityGraph::Arc& arc(*it_arc);
//		int var1=(**arc.head()).varindex;
//		int var2=(**arc.tail()).varindex;
//		if (var1<=var2) {
//			double hessval=hessian_in_refpoint[pair<int,int>(var2,var1)];
//			double bound=CoinMax(100., 10*CoinAbs(hessval));
//			collb[index]=-bound;
//			colub[index]=bound;
//			++index;
//		}
//	}


	DenseVector objcoeff(nr_cols);
	for (int i=nr_coeff; i<nr_cols; ++i) objcoeff[i]=1.;

	CoinPackedVectorBase** rows=new CoinPackedVectorBase*[nr_auxvars];
	DenseVector rowlb(nr_auxvars, -getInfinity());
	DenseVector rowub(nr_auxvars);

	int samplepoint_index=0;
	int enforce_tightness_index=-1;
	sampleset_newsort.clear();
	for (SampleSet::iterator it(func.getSamplePoints().begin()); it!=func.getSamplePoints().end(); ++it, ++samplepoint_index) {
		double scale;
		SparseVector* row=constructRow(func, *it, nr_coeff+samplepoint_index, scale);
		rows[samplepoint_index]=row;

		if (it==enforce_tightness) {
			colub[nr_coeff+samplepoint_index]=eps*CoinAbs(it->funcvalue/scale);
			sampleset_newsort.push_front(SampleSetItem(*it, nr_coeff+samplepoint_index, samplepoint_index, scale));
			enforce_tightness_index=nr_coeff+samplepoint_index;
		} else {
			sampleset_newsort.push_back(SampleSetItem(*it, nr_coeff+samplepoint_index, samplepoint_index, scale));
		}
	}

	lp.addCols(nr_cols, cols, collb.getElements(), colub.getElements(), objcoeff.getElements());
	delete[] cols;

	lp.addRows(nr_auxvars, rows, rowlb.getElements(), rowub.getElements());
	for (int i=0; i<nr_auxvars; ++i)
		delete rows[i];
	delete[] rows;

	lp.messageHandler()->setLogLevel(0);

	return enforce_tightness_index;
}


SparseVector* QuadraticEstimation::constructRow(NonconvexFunction& func, const DenseVector& point, int samplepoint_col, double& scale) {
	int rownz=nr_coeff;
	if (samplepoint_col>=0) ++rownz;
	int* indices=new int[rownz];
	CoinIotaN(indices, nr_coeff, 0); // put indices 0..nr_coeff-1
	double* elements=new double[rownz];

	double* el=elements;
	*el++=1.; // corresponding to constant term of quad. function

	// corresponding to linear term of poly of quad. function == coefficients of samplepoint
	assert(point.getNumElements()==func.dim());
	CoinMemcpyN(point.getElements(), func.dim(), el);
	el+=func.dim();

	for (SparsityGraph::arc_iterator it_arc(func.getSparsityGraph()->arc_begin()); it_arc!=func.getSparsityGraph()->arc_end(); ++it_arc) {
		const SparsityGraph::Arc& arc(*it_arc);
		int var1=(**arc.head()).varindex;
		int var2=(**arc.tail()).varindex;
		if (var1==var2) *el++=point[var1]*point[var1];
		else if (var1<var2) *el++=2*point[var1]*point[var2];
	}

	if (samplepoint_col>=0) {
		indices[nr_coeff]=samplepoint_col;
		*el=1.;
	}

	SparseVector* row=new SparseVector;
	row->assignVector(rownz, indices, elements, false);
	
//	scale = 1.;
	scale=row->infNorm();
//	printf("scale: %g\n", scale);
//	cout << *row << endl;
	if (scale<1.) scale=1.;
	if (scale>1.) *row/=scale; // scale row

	return row;
}

void QuadraticEstimation::updateLP(bool as_underestimator) {
	for (list<SampleSetItem>::iterator it(sampleset_newsort.begin()); it!=sampleset_newsort.end(); ++it) {
		const SamplePoint& point(it->sample_point);
		assert(point.funcvalue!=getInfinity());

		double rhs=point.funcvalue/it->maxcoeff;
		if (!as_underestimator) rhs=-rhs;
		// a sample point which has an own column need an equality constraint
		double lhs=(it->colnr>=0) ? rhs : -getInfinity();
		lp.setRowBounds(it->rownr, lhs, rhs);
	}
}

SmartPtr<QuadraticFunction> QuadraticEstimation::getEstimator(NonconvexFunction& func, bool as_underestimator, const DenseVector& lower, const DenseVector& upper) {
	SmartPtr<SymSparseMatrix> best_A;
	SmartPtr<SparseVector> best_b;
	double best_constant=0.;
	double best_violation=-getInfinity();
	double best_U3_val=getInfinity();

	SmartPtr<SparsityGraph> sparsitygraph(func.getSparsityGraph());
//	Timer timer, timer2;

	SmartPtr<QuadraticFunction> quadfunc=new QuadraticFunction(NULL, NULL, 0);

	// f - p if underestimator; -f - p if overestimator
	SumFunction fpdiff(as_underestimator ? 1 : -1, func.getFunction(), -1, GetRawPtr(quadfunc));

	// maximization of estimation error
	SmartPtr<BoxMinimizationProblem> errormaxprob=new BoxMinimizationProblem(fpdiff, lower, upper, sparsitygraph);

	int nr_locmin=0;
	bool finished;
	int iter=0;
	do {
		// remind Clp about the loglevel, maybe it forgot it already
		lp.messageHandler()->setLogLevel(0);

//		timer.start();
		if (iter==0) lp.initialSolve();
		else lp.resolve();
//		U3_time+=timer.stop();

		if (!lp.isProvenOptimal()) {
			cerr << "U3 not solved to optimality.";
			if (best_violation>-getInfinity()) {
				clog << " Keeping prior found solution with violation " << best_violation << ". ";
				break;
			} else {
				lp.writeMps("U3_fail");
				lp.writeLp("U3_fail");
				cerr << " Computation of quad. estimator failed!" << endl;
				return NULL;
//				exit(EXIT_FAILURE);
			}
		}
		double U3_val=lp.getObjValue();

		DenseVector U3sol(lp.getNumCols(), lp.getColSolution());
		DenseVector coeff(nr_coeff, U3sol.getElements());

		DenseVector rowact(lp.getNumRows(), lp.getRowActivity());

		// setup p(x)
		SymSparseMatrixCreator A_(func.dim());
		SparseVectorCreator b_;
		double constant=coeff[0];

		int coeff_index=1;
		for (int i=0; i<func.dim(); ++i, ++coeff_index)
			b_.insert(i, coeff[coeff_index]);

		for (SparsityGraph::arc_iterator it_arc(sparsitygraph->arc_begin()); it_arc!=sparsitygraph->arc_end(); ++it_arc) {
			const SparsityGraph::Arc& arc(*it_arc);
			int var1=(**arc.head()).varindex;
			int var2=(**arc.tail()).varindex;
			if (var1<=var2) A_.insert(var1, var2, coeff[coeff_index++]);
		}

		quadfunc->A=new SymSparseMatrix(A_);
		quadfunc->b=new SparseVector(b_);
		quadfunc->constant=constant;

//		clog << "quadfunc: " << *quadfunc << endl;

		finished=true;
		double maxviol=0;
		double maxviol_unscaled=0;
		// check active constraints to determine ''active'' sample points
		// and start local minimization of f-p from active sample points
		for (list<SampleSetItem>::iterator it(sampleset_newsort.begin()); it!=sampleset_newsort.end() && finished; ++it) {
			int j=it->rownr;
			// for original sample points we have a column that reports the distance between f and p
			if (j<nr_auxvars && U3sol(nr_coeff+j)>1E-4) continue;
			// for newer sample points the distance is given by the row activity
			double rhs=it->sample_point.funcvalue/it->maxcoeff;
			if (!as_underestimator) rhs=-rhs;
			if (j>=nr_auxvars && rowact[j]-rhs<-1E-4) continue;
// 			out_log << "U3 constraint active " << rowact(j)-*it_rhs << " for sample point " << it_sample_point->second;
//			if (!finished) continue;

			double viol1, viol2, scale2, f_val;
			SparseVector* new_row=NULL;
			++nr_locmin;
			if (!doLocMin(errormaxprob, func, it->sample_point, f_val, viol1, viol2, scale2, new_row, nr_locmin>1)) continue;

			if (maxviol>viol1) { maxviol=viol1; maxviol_unscaled=errormaxprob->getOptimalValue(); }

			if (viol1<-1E-2 && viol2<-1E-2) { // too large, add point to sample set and LP and restart
				pair<set<SamplePoint>::iterator, bool> newpoint(func.getSamplePoints().insert(errormaxprob->getSolution()));
				if (newpoint.second) {
					newpoint.first->funcvalue=f_val;
					double buffer=eps*CoinMin(1., sqrt(newpoint.first->getPoint().euclidianDistance(it->sample_point)/CoinMax(1., it->sample_point.twoNorm())))*CoinAbs(f_val/scale2);
					rhs=f_val/scale2-buffer;
					if (!as_underestimator) rhs=-rhs;
					lp.addRow(*new_row, -getInfinity(), rhs);
					finished=false;
// 	 				out_log << "\tadded row with rhs " << rhs.back() << " and coeff " << b1;
					sampleset_newsort.push_front(SampleSetItem(*newpoint.first, -1, lp.getNumRows()-1, scale2));
				} else {
					clog << "Local minimizer already in sample set, not adding again." << endl;
				}
			}

			delete new_row;

//			if (time_max && timer.stop()>time_max && locminsolver.iter_max>10*f->dim()) {
//				if (best_violation==-INFINITY) { // first iteration -> give him some chance to check the other points
//					out_log << "Time limit (" << time_max << ") exceeded. Reduce locmin iteration limit to " << 10*f->dim();
//					locminsolver.iter_max=10*f->dim();
//				} else { // some U3 solution already checked -> abort
//					out_log << "Time limit (" << time_max << ") exceeded. ";
//					maxviol=-INFINITY;
//					break;
//				}
//			}
		}
		clog << 'v' << -maxviol_unscaled << ' ';

		if (maxviol_unscaled>best_violation || finished) {
			best_violation=maxviol_unscaled;
			best_A=quadfunc->A;
			best_b=quadfunc->b;
			best_constant=quadfunc->constant;
			best_U3_val=U3_val;
		}

		if (iter>=iter_max) {
			clog << "Iteration limit (" << iter_max << ") exceeded. Keeping solution with violation " << best_violation << ' ';
			finished=true;
		}
//		if (time_max && timer.stop()>time_max) {
//			out_log << "Keeping solution with violation " << best_violation << ' ';
//			finished=true;
//		}

		++iter;
	} while (!finished);

	best_constant+=best_violation; // lower underestimator by maximum known violation
	clog << "U3:" << best_U3_val << ' ';

//	// set A, b, and c.
//	int i=0;
//	list<MultiIndex>::iterator it_mind(multiindices.begin());
//	if (it_mind!=multiindices.end() && it_mind->size()==0) {
//		if (fabs(best_coeff[i])>max_abscoeff) max_abscoeff=fabs(best_coeff[i]);
//		c+=best_coeff[i];
//		it_mind++;
//		i++;
//	}
//
//	while (it_mind!=multiindices.end() && it_mind->size()==1) {
//		if (fabs(best_coeff[i])>max_abscoeff) max_abscoeff=fabs(best_coeff[i]);
//		b[indices(*it_mind->begin())]+=.5*best_coeff[i];
//		it_mind++;
//		i++;
//	}
//
//	while (it_mind!=multiindices.end()) {
//		if (fabs(best_coeff[i])>max_abscoeff) max_abscoeff=fabs(best_coeff[i]);
//		int second=*(++it_mind->begin());
//		if (*it_mind->begin()!=second) {
//		  A.AddToElement(indices(*it_mind->begin()), indices(second), .5*best_coeff[i]);
//		  A.AddToElement(indices(second), indices(*it_mind->begin()), .5*best_coeff[i]);
//		} else {
//		  A.AddToElement(indices(second), indices(second), best_coeff[i]);
//		}
//		it_mind++;
//		i++;
//	}

//	if (nr_locmin>max_locmin) max_locmin=nr_locmin;

	quadfunc->A=best_A;
	quadfunc->b=best_b;
	quadfunc->constant=best_constant;
	
	if (!as_underestimator) {
		quadfunc->A->scale(-1);
		quadfunc->b->scale(-1);
		quadfunc->constant*=-1;
	}

	return quadfunc;
}

bool QuadraticEstimation::doLocMin(SmartPtr<BoxMinimizationProblem>& prob, NonconvexFunction& func, const SamplePoint& sample_point, double& f_val, double& viol1, double& viol2, double& scale2, SparseVector*& row, bool do_resolve) {
//	Timer timer;
	prob->startpoint=&sample_point.getPoint();

	SmartPtr<Ipopt::TNLP> tnlp(GetRawPtr(prob));
	Ipopt::ApplicationReturnStatus ret;
	if (do_resolve)
		ret=ipopt.ReOptimizeTNLP(tnlp);
	else
		ret=ipopt.OptimizeTNLP(tnlp);

	if (ret!=Ipopt::Solve_Succeeded)
		clog << 'm' << ret;
	else
		clog << 's';

	switch (ret) {
		case Ipopt::Solve_Succeeded:
		case Ipopt::Solved_To_Acceptable_Level:
		case Ipopt::Maximum_Iterations_Exceeded:
			break;
		default:
			return false;

	}

	f_val=func.getFunction()->eval(prob->getSolution());
	viol1=prob->getOptimalValue()/CoinMax(1., CoinAbs(f_val));

	viol2=0.;
	scale2=1.;
	if (viol1<-1E-4) { // compute violation according to second (LP) scaling
		row=constructRow(func, prob->getSolution(), -1, scale2);
//		if (scale2>1.) *new_row/=scale2; // scale row
		viol2=prob->getOptimalValue()/scale2;
	}

//	locopt_time+=timer.stop();

	return true;
}

void QuadraticEstimation::testEstimators(const MINLPData& data) const {
	for (int c=0; c<=data.numConstraints(); ++c) {
		const MINLPData::ObjCon& objcon(c<data.numConstraints() ? (MINLPData::ObjCon)data.con[c] : (MINLPData::ObjCon)data.obj);
		for (int k=0; k<(int)objcon.decompfuncNL.size(); ++k) {
			const BlockFunction& blockfunc(*objcon.decompfuncNL[k]);
			if (IsNull(blockfunc.nonquad)) continue;
			DenseVector low, up;
			data.getBox(low, up, blockfunc.indices);
			for (list<SmartPtr<QuadraticEstimator> >::const_iterator it(blockfunc.underestimators.begin()); it!=blockfunc.underestimators.end(); ++it)
				testEstimator(*blockfunc.nonquad, *(*it)->func, true, low, up);
			for (list<SmartPtr<QuadraticEstimator> >::const_iterator it(blockfunc.overestimators.begin()); it!=blockfunc.overestimators.end(); ++it)
				testEstimator(*blockfunc.nonquad, *(*it)->func, false, low, up);
		}
	}
}
	
void QuadraticEstimation::testEstimator(const Function& orig, const QuadraticFunction& estimator, bool is_underestimator, const DenseVector& lower, const DenseVector& upper) const {
	DenseVector x;
	for (int i=0; i<20; ++i) {
		x.setRandom(lower, upper);
		double origval=orig.eval(x);
		double estval=estimator.eval(x);
		double diff=(origval-estval)/(1+CoinAbs(origval));
		if (is_underestimator) diff*=-1;
		if (diff>1e-4) {
			cerr << "QuadraticEstimation test: " << (is_underestimator ? "underestimator" : "overestimator") << " of function " << orig << endl;
//			cerr << "\t Estimator: " << estimator << endl;
			cerr << "\t Violation by rel. dist. " << diff << "\t origval: " << origval << "\t estval: " << estval << endl;
		}
	}
	
}   


} // namespace LaGO
