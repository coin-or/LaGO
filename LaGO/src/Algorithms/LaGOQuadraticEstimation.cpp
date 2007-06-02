// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOQuadraticEstimation.hpp"
#include "LaGOSampling.hpp"

namespace LaGO {
	
void QuadraticEstimation::computeEstimators() {
	computeEstimators(data.obj, true, false); // compute underestimators of objective
	for (int c=0; c<data.numConstraints(); ++c)
		computeEstimators(data.con[c], data.con[c].upper<getInfinity(), data.con[c].lower>-getInfinity());
}
	
void QuadraticEstimation::computeEstimators(MINLPData::ObjCon& objcon, bool need_lower, bool need_upper) {
	for (unsigned int k=0; k<objcon.decompfuncNL.size(); ++k) {
		BlockFunction& func(*objcon.decompfuncNL[k]);
		if (IsNull(func.nonquad)) continue;
		lp.reset();
		clog << "Quadratic Estimators " << objcon.name << " block " << k << '(' << need_lower << need_upper << "): ";
		computeEstimator(func, need_lower, need_upper);
		clog << endl;
	}	
}

void QuadraticEstimation::computeEstimator(BlockFunction& func, bool need_lower, bool need_upper) {
	bool do_lower=need_lower && !(func.curvature&CONVEX);
	bool do_upper=need_upper && !(func.curvature&CONCAVE);
	clog << do_lower << do_upper;
	if ((!do_lower) && (!do_upper)) return;
	
	DenseVector lower, upper;
	data.getBox(lower, upper, func.indices); 

	Sampling sampling;
	
	// set sample point function values 
	SampleSet::iterator it_sp(func.samplepoints.begin());
	while (it_sp!=func.samplepoints.end()) {
		try {
			if (it_sp->funcvalue==getInfinity())
				it_sp->funcvalue=func.nonquad->eval(*it_sp);
		} catch (FunctionEvaluationError error) {
			cerr << "QuadraticEstimation::computeEstimator: " << error << endl;
			it_sp=func.samplepoints.eraseAndGetNext(it_sp);
			continue; 
		}
		++it_sp;
	}
	
	if (do_lower) {
		SampleSet::iterator enforce_tightness=sampling.addMinimizer(func.samplepoints, func, lower, upper, func.sparsitygraph, &func.samplepoints.begin()->getPoint());
		if (enforce_tightness==func.samplepoints.end()) {
			enforce_tightness=func.samplepoints.begin();
			while (enforce_tightness!=func.samplepoints.end() && !enforce_tightness->is_startpoint)
				++enforce_tightness;
			if (enforce_tightness==func.samplepoints.end())
				enforce_tightness=func.samplepoints.begin();
		}

		initLP(func, enforce_tightness);
		
		updateLP(func, true);
		
//		getEstimator(func);
	}
	
	if (do_upper) {
		//TODO generate maximizer 
		SampleSet::iterator enforce_tightness;
		if (!do_lower)
			initLP(func, enforce_tightness);
		else {
			// add enforce_tightness at front of sampleset_newsort and as row
		}			
			
		updateLP(func, false);
		
	}
	
}

void QuadraticEstimation::initLP(BlockFunction& func, const SampleSet::iterator& enforce_tightness) {
// create auxiliary LP with initial sample set
	int nr_coeff=1+func.indices.size()+func.sparsitygraph->arc_size();
//	int multiindices_size=multiindices.size(); // sparsity of p(x)
	int nr_auxvars=func.samplepoints.size();
	int nr_cols=nr_coeff+nr_auxvars;

	SparseVector zero;
	CoinPackedVectorBase** cols=new CoinPackedVectorBase*[nr_cols];
	for (int i=0; i<nr_cols; ++i) cols[i]=&zero;

	DenseVector collb(nr_cols, 0.);
	DenseVector colub(nr_cols, getInfinity());
	for (int i=0; i<nr_coeff; ++i) collb[i]=-getInfinity();
	
	DenseVector objcoeff(nr_cols);
	for (int i=nr_coeff; i<nr_cols; ++i) objcoeff[i]=1.;

	lp.addCols(nr_cols, cols, collb.getElements(), colub.getElements(), objcoeff.getElements());
	delete[] cols;
	
	CoinPackedVectorBase** rows=new CoinPackedVectorBase*[nr_auxvars];
	DenseVector rowlb(nr_auxvars, -getInfinity());
	DenseVector rowub(nr_auxvars);	
	
	int samplepoint_index=0;
	sampleset_newsort.clear();
	for (SampleSet::iterator it(func.samplepoints.begin()); it!=func.samplepoints.end(); ++it, ++samplepoint_index) {
		SparseVector* row=constructRow(func, *it, nr_coeff+samplepoint_index);
		rows[samplepoint_index]=row;
		
		double scale=row->infNorm();
		if (scale<1.) scale=1.;
		if (scale>1.) { // scale row, except for last element (samplepoint_index)
			for (int i=0; i<nr_coeff; ++i)
				row->getElements()[i]/=scale;
		}
		
		if (it==enforce_tightness)
			sampleset_newsort.push_front(SampleSetItem(*it, nr_coeff+samplepoint_index, samplepoint_index, scale));
		else
			sampleset_newsort.push_back(SampleSetItem(*it, nr_coeff+samplepoint_index, samplepoint_index, scale));
	}
	
	lp.addRows(nr_auxvars, rows, rowlb.getElements(), rowub.getElements());
	
	for (int i=0; i<nr_auxvars; ++i)
		delete rows[i];
	delete[] rows;
}


SparseVector* QuadraticEstimation::constructRow(BlockFunction& func, const SamplePoint& point, int samplepoint_col) {
	int nr_coeff=1+func.indices.size()+func.sparsitygraph->arc_size();
	
	DenseVector point_(point.getPoint());

	int* indices=new int[nr_coeff+1];
	CoinIotaN(indices, nr_coeff, 0); // put indices 0..nr_coeff-1
	double* elements=new double[nr_coeff+1];
	
	double* el=elements;	
	*el++=1.; // corresponding to constant term of quad. function
	
	// corresponding to linear term of poly of quad. function == coefficients of samplepoint
	assert(point_.getNumElements()==(int)func.indices.size());
	CoinMemcpyN(point_.getElements(), func.indices.size(), el);
	el+=func.indices.size();
	
	for (SparsityGraph::arc_iterator it_arc(func.sparsitygraph->arc_begin()); it_arc!=func.sparsitygraph->arc_end(); ++it_arc, ++el) {
		const SparsityGraph::Arc& arc(*it_arc);
		int var1=(**arc.head()).varindex;
		int var2=(**arc.tail()).varindex;
		if (var1==var2) *el=point_[var1]*point_[var1];
		else *el=2*point_[var1]*point_[var2];
	}
	
	indices[nr_coeff]=samplepoint_col;
	*el=1.;
	
	SparseVector* row=new SparseVector;
	row->assignVector(nr_coeff+1, indices, elements, false);	
	
	return row; 
}

void QuadraticEstimation::updateLP(BlockFunction& func, bool as_underestimator) {
	for (list<SampleSetItem>::iterator it(sampleset_newsort.begin()); it!=sampleset_newsort.end(); ++it) {
		const SamplePoint& point(it->sample_point);
		assert(point.funcvalue!=getInfinity());
		
		double rhs=point.funcvalue/it->maxcoeff;
		if (!as_underestimator) rhs=-rhs;
		if (it==sampleset_newsort.begin()) { // we wanna enforce tightness in this sample point, at least almost tightness
			lp.setColBounds(it->colnr, 0, eps*CoinAbs(rhs));
		} else {
			lp.setColBounds(it->colnr, 0, getInfinity());
		}		
		lp.setRowBounds(it->rownr, rhs, rhs);
	}
}
	
	
} // namespace LaGO
