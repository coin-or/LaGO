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
		@param Quadratic Underestimator iteration limit
		%options >0
		%default 10
		%level 1
		How often (at most) we should refit the underestimator of a function.
		@param Quadratic Underestimator time limit
		%options >=0
		%default 0
		%level 1
		How long (at most) we should refit the underestimator of a function. 0 means no timelimit.
		@param Quadratic Underestimator epsilon
		%options >=0
		%default 0.0001
		%level 0
		How big should the buffer in the underestimator for newly computed sample points be.
		@param Quadratic Underestimator sample set Monte Carlo
		%options integer >= 0
		%default 20
		%level 1
		@param Quadratic Underestimator sample set mid point
		%options 0, 1
		%default 0
		%level 1
		@param Quadratic Underestimator sample set box ends
		%options 0, 1
		%default 0
		%level 1
		@param Quadratic Underestimator sample set vertices2
		%options integer >= 0
		%default 200
		%level 1
		@param Quadratic Underestimator sample set minimizer
		%options 0, 1
		%default 1
		%level 1
		@param Quadratic Underestimator sample set initial
		%options 0, 1
		%default 1
		%level 1
*/
class QuadraticUnderestimator {
	private:
		Pointer<Param> param;
		
		double eps;
		int iter_max;
		double time_max;

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
		multimap<double, dvector>::iterator add_minimizer_to_sample(Pointer<Func> f, const dvector& lower, const dvector& upper, dvector& start);
		
		/** new quadratic underestimator, based on A. Neumaiers idea
		*/
		void quadratic_underestimator(SparseMatrix2& A, SparseVector<double>& b, double& c, const Pointer<Func>& f, ivector& indices, const Pointer<dvector>& lower, const Pointer<dvector>& upper);
		
		bool do_locmin(BoxLocOpt& locminsolver, const Pointer<Func>& f, dvector& start, double& f_val, double& viol1, double& viol2, double& scale2, dvector& b1);
		
	public:
		SparseVector<double> c_add1, c_add2;
		
		double U3_time, locopt_time, max_abscoeff;
		int max_locmin, nr_estimators;

		QuadraticUnderestimator(Pointer<Param> param_=NULL);

		void quadratic_underestimator(MinlpProblem& prob, MINLPData& minlpdata, ivector& ineq_index, SparseVector<double>& obj_c_add, vector<SparseVector<double> >& con_c_add);
		Pointer<SepQcFunc> quadratic_underestimator(Pointer<SepQcFunc> f, bool eq, dvector& lower, dvector& upper, dvector& primal);
};

#endif
