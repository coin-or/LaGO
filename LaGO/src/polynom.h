// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef POLYNOM_H
#define POLYNOM_H

#include "standard.h"
#include "func.h"
#include "sampling.h"
#include "decomp.h"
#include "MINLPData.h"

/** Class to represent a MultiIndex.
*/
class MultiIndex : public multiset<int> {
	friend ostream& operator<<(ostream& out, const MultiIndex& x) {
		for (MultiIndex::const_iterator it(x.begin()); it!=x.end(); it++) out << (it==x.begin() ? "" : " ") << *it;
		return out;
	}
	
	public:
		/** Default-Constructor.
		*/
		MultiIndex()
		: multiset<int>()
		{ }
		
		/** Constructor for one element.
		    Constructs a multiset with given element.
		    @param i Element to put in multiset.
		*/
		MultiIndex(int i)
		: multiset<int>()
		{ insert(i);
		}
		
		/** Constructor for two elements.
		    Constructs a multiset with two given element.
		    @param i First element to put in multiset.
		    @param j Second element to put in multiset.
		*/
		MultiIndex(int i, int j)
		: multiset<int>()
		{ insert(i);
		  insert(j);
		}
		
		/** Comparision operator.
		    @param x The MultiIndex to compare this MultiIndex with.
		    @return If this->size() < x.size(), returns true. If this->size()>x.size(), returns false. If both have same size, calls operator< from multiset-class, which does a lexicographic comparision.
		*/
		bool operator<(const MultiIndex& x) {
//			if (size()==x.size()) return ::operator<(*(const multiset<int>*)this, *(const multiset<int>*)&x);
			if (size()==x.size()) return std::operator<(*this, x);
			if (size()<x.size()) return true;
			return false;
		}
	
};

/** Class to represent a Monom.
*/
class Monom : public Func {
	private:
		MultiIndex indices;

		/** Recursive algorithm to compute a partial derivate of a Monom.
				Computes @f$\frac{\partial x^\alpha}{\partial x^{\beta_1}}(x)@f$
		    @param x The UserVector, for which we want to compute the derivative.
		    @param alpha The indices, which are left in our Monom.
		    @param beta The indices, we need to derive for.
		    @return The derivative @f$\frac{\partial^{|\beta|} x^\alpha}{\partial x^\beta}(x)@f$
		*/
		double part_derivate_rek(const UserVector<double>& x, MultiIndex& alpha, MultiIndex& beta) const;
		
	public:
		/** Constructor for a MultiIndex of indices.
		    @param dim_ The dimension of this function.
		    @param indices_ The indices, for which this class represents the monom. If an empty set is given, this monom will be the constant function 1.
		*/
		Monom(int dim_, MultiIndex& indices_)
		: Func(dim_), indices(indices_)
		{ }
		
		
//		Monom(int dim_, MultInd& indices_);	

		double eval(const UserVector<double>& x) const {
			double val=1.;
			for (MultiIndex::const_iterator it(indices.begin()); it!=indices.end(); it++) val*=x(*it);
			return val;
		}
		
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Func::eval;
#endif

		void grad(UserVector<double>& grad, const UserVector<double>& x) const;
		
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Func::grad;
#endif

		/** Computes the multiplication of the hessian of this Monom with a UserVector.
				Sets @f$y_i = \sum_{j,z_j\neq 0} \frac{\partial *this}{\partial x_i x_j}(x) z_j@f$
		    @param y The UserVector to store the result in.
		    @param x The UserVector to compute the hessian for.
		    @param z The UserVector to multiply with.
		*/
		void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const;
		
#if (!defined(__GNUC__)) || (GCC_VERSION>=3000)
		using Func::HessMult;
#endif

		/** Gives a partial derivate in a point.
		    @param x The point to compute the derivative for.
		    @param ind The index-set of variables, which we want to derive.
		    @return @f$\frac{\partial^{|\alpha|} *this}{\partial \alpha_1\cdots \partial \alpha_{|\alpha|}}@f$
		    @see part_derivate_rek(const UserVector<double>&, MultiIndex&, MultiIndex&) const;
		*/
		double part_derivate(const UserVector<double>& x, const MultiIndex& ind) const {
			MultiIndex alpha(indices);
			MultiIndex beta(ind);
			return part_derivate_rek(x, alpha, beta);
		}

		virtual void set_curvature(CurvatureType ct) { out_err << "Monom::set_curvature() not implemented. Aborting." << endl; exit (-1); };
		virtual CurvatureType get_curvature() const  { out_err << "Monom::get_curvature() not implemented. Aborting." << endl; exit (-1); return Func::UNKNOWN; };

		void print(ostream& out) const {
			out << "[" << indices << "]";
		}
};


/** Class to compute Polynomial underestimators of nonconvex functions.
    @class PolynomialUnderestimator2
    @param Polynomial Underestimator K? sample set ...
		%level 1
    All sample set parameters starts with "Polynomial Underestimator Kx" here, with x=0,1,2. Check Sampling section for additional parameters.
    @param Polynomial Underestimator K? sample set initial
    %options 0 or 1
    %default 0
		%level 1
    If 1, adds the point given in add\_point\_to\_sample to the sample set K?.
    @param Polynomial Underestimator max polynom degree
    %options integer $\geq 0$
    %default 2
    The maximum degree, a polynomial understimator is allowed to have.
*/
class PolynomialUnderestimator2 {
	private:
		Pointer<Param> param;

		int max_degree;
		int maxdegree1_size, maxdegree2_size;

		Sampling sampling0; // for K0
		Sampling sampling1; // for K1
		Sampling sampling2; // for K2
		Sampling_Vertices sampling_vertices;
		Sampling_Minimizer sampling_minimizer;
		int sampling_initial; // should we add the initial point
		vector<dvector> sample_set;
		ivector ss_size;

		list<MultiIndex> multiindices;
		list<Monom> monoms;

		Decomposition decomp;

		void new_multiindices(const SparsityInfo& si, int n);
		void polynomial_underestimator(SparseMatrix2& A, SparseVector<double>& b, double& c, Func& f, ivector& indices);

	public:
		SparseVector<double> c_add1, c_add2;

		PolynomialUnderestimator2(Pointer<Param> param_=NULL);

		void polynomial_underestimator(MinlpProblem& prob, MINLPData& minlpdata, ivector& ineq_index, SparseVector<double>& obj_c_add, vector<SparseVector<double> >& con_c_add);
		Pointer<SepQcFunc> polynomial_underestimator(Pointer<SepQcFunc> f, bool eq, dvector& lower, dvector& upper, dvector& primal);
		void new_sampleset(const dvector& lower, const dvector& upper);
		void check_for_nan(const Func& f);
		bool add_point_to_sampleset(const dvector& point);
		bool add_minimizer_to_sample(Pointer<Func> f, const dvector& lower, const dvector& upper, dvector& start);
		void remove_last_point_from_sample();

		void check(MinlpProblem& prob, MinlpProblem& quad, ivector& ineq_index);
};

#endif
