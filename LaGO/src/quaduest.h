// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef QUADUEST_H
#define QUADUEST_H

#include "standard.h"
#include "func.h"
#include "sampling.h"
#include "decomp.h"
#include "MINLPData.h"
#include "polynom.h"

/** Class to compute a quadratic underestimators of nonconvex functions.
		@class QuadraticUnderestimator
		@param Quadratic Underestimator sample set ...
		%level 1
		All sample set parameters starts with "Quadratic Underestimator" here. Check Sampling section for additional parameters.
		@param Quadratic Underestimator sample set initial
		%options 0 or 1
		%default 0
		%level 1
		@param Quadratic Underestimator sample set initial
		%options 0 or 1
		%default 0
		%level 1
		@param Quadratic Underestimator iteration limit
		%options >0
		%default 100
		%level 1
		How often we should refit the underestimator of a function.
		@param Quadratic Underestimator epsilon
		%options >=0
		%default 0.001
		%level 0
		How big should the buffer in the underestimator for newly computed sample points be.
*/
class QuadraticUnderestimator {
	private:
		Pointer<Param> param;
		
		double eps;
		int iter_max;

		Sampling sampling;
		Sampling_Vertices sampling_vertices;
		Sampling_Minimizer sampling_minimizer;
		int sampling_initial; // should we add the initial point
		
		multimap<double, dvector> sample_set;
		multimap<double, dvector>::iterator enforce_tightness;

		list<MultiIndex> multiindices;
		list<Monom> monoms;

		Decomposition decomp;

		void new_multiindices(const SparsityInfo& si, int n);
		
		void new_sampleset(const dvector& lower, const dvector& upper);
		void check_for_nan(const Func& f);
		multimap<double, dvector>::iterator add_point_to_sampleset(const dvector& point);
		multimap<double, dvector>::iterator add_minimizer_to_sample(Pointer<Func> f, const dvector& lower, const dvector& upper);
		
		/** new quadratic underestimator, based on A. Neumaiers idea
		*/
		void quadratic_underestimator(SparseMatrix2& A, SparseVector<double>& b, double& c, const Pointer<Func>& f, ivector& indices, const Pointer<dvector>& lower, const Pointer<dvector>& upper);
		
	public:
		SparseVector<double> c_add1, c_add2;

		QuadraticUnderestimator(Pointer<Param> param_=NULL);

		void quadratic_underestimator(MinlpProblem& prob, MINLPData& minlpdata, ivector& ineq_index, SparseVector<double>& obj_c_add, vector<SparseVector<double> >& con_c_add);
		Pointer<SepQcFunc> quadratic_underestimator(Pointer<SepQcFunc> f, bool eq, dvector& lower, dvector& upper, dvector& primal);
};

#endif
