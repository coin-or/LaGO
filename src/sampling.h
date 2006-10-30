// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef SAMPLING_H
#define SAMPLING_H

#include "standard.h"
#include "problem.h"
#include "func.h"
#include "opt.h"

/** Class to generate a sample set.
    @class Sampling
    @param sample set Monte Carlo
    %options integer $\geq0$
    %default 5
		%level 1
    The number of random points, which are added to the sample set. Works also for unbounded variables.
    @param sample set mid point
    %options 0 or 1
    %default 1
		%level 1
    If 1, the middle point of the box is added to the sample set. If the box is incomplete, this is not done.
    @param sample set box ends
    %options 0 or 1
    %default 1
		%level 1
    If 1, the lower and upper bounds of the box are added to the sample set. If the box is incomplete, this is not done.
    @param sample set vertices
    %options integer $\geq0$
    %default 5
		%level 1
    The number of random points, which lie on the vertices. If the box is incomplete, this is not done.

*/
class Sampling {
	/** Output operator.
	    Calls cr.print(out).
	    @param out The ostream to print to.
	    @param cr The Relax to print.
	    @return out
	*/
	friend ostream& operator<<(ostream& out, const Sampling& s) {
		s.print(out);
		return out;
	}

	protected:
		/** Number of random sample points to generate.
		*/
		int montecarlo;
		/** Number of random points at the vertices to generate.
		*/
		int vertices;
		/** Indicates, whether the midpoint should be added.
		*/
		bool midpoint;
		/** Indicates, whether the lower and upper bounds should be added.
		*/
		bool bounds;

	public:
		/** (Default-)Constructor.
		    @param param Parameters.
		    @param param_prefix A prefix for the parameter names.
		*/
		Sampling(Pointer<Param> param=NULL, char* param_prefix=NULL);
		
		/** Destructor.
		    Does nothing.
		*/
		virtual ~Sampling() { }
		
		/** Adds sample points to a sample set.
		    @param sample_set The sample set to add the points to.
		    @param lower Lower bounds of variables.
		    @param upper Upper bounds of variables.
		    @return The number of added sample points.
		*/
		virtual int get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper);
		
		/** Adds sample points to a sample set, given a block-structure.
		    @param sample_set The sample set to add the points to. sample_set.size() must be the number of blocks.
		    @param lower Lower bounds of variables.
		    @param upper Upper bounds of variables.
		    @param block The block structure.
		*/
		void get_points(vector<vector<dvector> >& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, const vector<ivector>& block) {
			assert(sample_set.size()==block.size());
			for (int k=0; k<block.size(); k++) get_points(sample_set[k], lower(block[k]), upper(block[k]));
		}

		/** Gives a sample set.
		    @param lower Lower bounds of variables.
		    @param upper Upper bounds of variables.
		    @param block A block structure.
		    @return The computed sample set.
		*/
		Pointer<vector<vector<dvector> > > get_points(const UserVector<double>& lower, const UserVector<double>& upper, const vector<ivector>& block) {
			Pointer<vector<vector<dvector> > > ss(new vector<vector<dvector> >(block.size()));
			get_points(*ss, lower, upper, block);
			return ss;
		}

		/** Prints parameters.
		    @param out The ostream to print to.
		*/
		virtual void print(ostream& out) const {
			out << "Montecarlo points: " << montecarlo << endl
			    << "Vertices points  : " << vertices << endl
			    << "Midpoint         : " << (midpoint ? "yes" : "no") << endl
		  	  << "Bounds           : " << (bounds ? "yes" : "no") << endl
			;
		}

};

/** Sampling technique to generate points at vertices.
    @class Sampling_Vertices
    @param sample set vertices2
    %options integer $\geq 0$
    %default 0
		%level 1
    The number of vertices-points to generate.
*/
class Sampling_Vertices {
	public:
		/** The number of vertices to generate.
	      Set in constructor according to parameters.
		*/
		int vertices;

		/** (Default-)Constructor.
		    @param param Parameters.
		    @param param_prefix A prefix for the parameter names.
		*/
		Sampling_Vertices(Pointer<Param> param=NULL, char* param_prefix=NULL);

		void get_points(vector<vector<dvector> >& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, const vector<ivector>& block, const vector<Pointer<set<int> > >& i_quad, const vector<Pointer<set<int> > > &i_nonquadlin);

		void get_points(vector<vector<dvector> >& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, const SepQcFunc& f);

		/** Generates points, which lie on the vertices with respect to a set of variable indices.
		    Generates points, which values lie on the bounds for the indices, given in i_quad and i_nonquadlin.
		    The other variables are set to the midpoint.
		    @param sample_set The sample set to add the points to.
		    @param lower The lower bounds of the variables.
		    @param upper The upper bounds of the variables.
		    @param i_quad The indices of the quadratic variables.
				@param i_nonquadlin The indices of the nonquadratic/nonlinear variables.
				@return The number of added points.
		*/
		int get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, set<int>& i_quad, set<int>& i_nonquadlin);

		int get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper, const SparsityInfo& si);

		int get_points(vector<dvector>& sample_set, const UserVector<double>& lower, const UserVector<double>& upper);
};

/** Simple algorithm to solve the equation f==0 over a box.
    Starting from a point, finds the other extreme of the box towards gradient direction, and do bisection to find the zero inside the box.
    @class SimpleEquationSolver
    @param SimpleEquationSolver max iter
    %options integer $\geq 0$
    %default 4
    Maximum number of gradient steps.
    @param SimpleEquationSolver tolerance
    %options double $\geq 0$
    %default 1E-4
    Function value tolerance. If absolut value of point is <1E-4, we accept the point.
    @param SimpleEquationSolver bisection max iter
    %options  integer $\geq 0$
    %default 1000
    Maximum number of iterations in bisection algorithm.
*/
class SimpleEquationSolver : public Solver {
	private:
		const Func& f;
		const UserVector<double>& lower;
		const UserVector<double>& upper;

		int bisection_iter_max;

		/** Determines point at other extreme of box and tries bisections-algorithm, to find zero between them.
	      @param start Starting point.
	      @param val Value at starting point.
	      @param dir The direction.
	      @return 0, if we found a zero; 1, if functions values of current point and point at other side of the box have the same sign; 2, if we got a nan in the gradient or as function value; 3, if the iteration limit of the bisection method was reached; 4, if we couldn't compute an initial steplength
		*/
		int run(const dvector& start, double val, dvector& dir);

	public:
		SimpleEquationSolver(const Func& f_, const UserVector<double>& lower_, const UserVector<double>& upper_, Pointer<Param> param=NULL)
		: Solver(f_.dim()), f(f_), lower(lower_), upper(upper_)
		{ iter_max=param ? param->get_i("SimpleEquationSolver max iter", 4) : 4;
		  tol=param ? param->get_d("SimpleEquationSolver tolerance", 1e-4) : 1e-4;
			bisection_iter_max=param ? param->get_i("SimpleEquationSolver bisection max iter", 1000) : 1000;
		}

		int solve(dvector& start);

		int solve() {
			dvector start(lower); start+=upper; start*=.5;
			return solve(start);
		}


};

/** Class to check an existing sample set, if a constraint is fullfilled.
    If the constraint is not fullfilled in a point, the point is tried to make feasible by adjusting the linear variables.
*/
class Sampling_check {
	private:
		Pointer<Param> param;

	public:
		/** (Default-)Constructor.
		*/
		Sampling_check(Pointer<Param> param_)
		: param(param_)
		{ }

		/** Computes sample points by sampling in one dimension less and adjusting the linear variable, such that a constraint is fullfilled.
		    If j is a linear variable, samples in {1,..,j-1,j+1,..}. If x_j can be set, such that f(x)==0 or f(x)<=0 and x_j is in its bounds, the point is added to the sample set.
		    If x_j cannot set this way, it's set to on of its bounds and the next variable is tried.
		    This method prints a warning, if no points were added.
				If no linear variables exists, it generates a normal sample set.
		    @param sample_set The sample sets for each block of f to add the points to.
		    @param f The function for the constraint f(x) == 0 or f(x) <= 0.
		    @param eq Indicates, whether the function is an equation or inequation.
		    @param i_lin The indices of the linear variables for each block.
		    @param lower The lower bounds of the variables.
		    @param upper The upper bounds of the variables
		    @param start The first sample point in the sample set, we should check.
				@param min_keep The minimum number of sample points to keep in the sample set, if somehow possible.
		*/
		void check(vector<vector<dvector> >& sample_set, const SepQcFunc& f, bool eq, const vector<Pointer<set<int> > >& i_lin, const UserVector<double>& lower, const UserVector<double>& upper, const vector<int>& start, int min_keep=0);

		void check(vector<vector<dvector> >& sample_set, const SepQcFunc& f, bool eq, const UserVector<double>& lower, const UserVector<double>& upper, const vector<int>& start, int min_keep);
};

/** Computing sample point by minimizing a given function, starting from the point from a sample set with the lowest sample function value.
    @class Sampling_Minimizer
    @param sample set minimizer
    %options 0 or 1
    %default 0
		%level 1
    If 1, tries to find a minimizer of the function and adds this one to the sample set.
    @param Sampling Minimizer snopt nospecs
    %options 0 or 1
    %default 0
    @param Sampling Minimizer snopt specs
    %options filename
    %default resource/superbasic.snopt
*/
class Sampling_Minimizer {
	private:
		bool minimizer;
		Pointer<Param> param;

	public:
		/** (Default-)Constructor.
		    @param param Parameters.
		    @param param_prefix A prefix for the parameter names.
		*/
		Sampling_Minimizer(Pointer<Param> param_=NULL, char* param_prefix=NULL);

		/** Uses the best sample point to start the search of a minimizer of a function to add this to the sample set.
		    @param sample_set A sample set to pick a starting point from.
		    @param f The function to minimize.
		    @param lower Lower bounds of variables.
		    @param upper Upper bounds of variables.
		    @return True, if a minimizer was found and added. False, else.
		*/
		bool add_minimizer(vector<vector<dvector> >& sample_set, Pointer<SepQcFunc> f, const UserVector<double>& lower, const UserVector<double>& upper);

		bool add_minimizer(vector<dvector>& sample_set, Pointer<Func> f, const UserVector<double>& lower, const UserVector<double>& upper);
};

#endif
