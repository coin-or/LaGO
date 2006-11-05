// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef NODE_H
#define NODE_H

#include "standard.h"
#include "problem.h"
#include "minlpopt.h"
#include "bcp.h"

class ColumnGenerator;
class RMPManager;
class LinearRelax;
class LinearRelaxSolverGeneral;
class LinearRelaxSolverMIP;
class LagHeu;
class LagHeu1;
class LagHeu_SimAnnealing;
class LagHeu2;
class LagHeu2b;

/** A class containing methods needed for partitioning algorithms.
    @class MinlpNode
*/
class MinlpNode {
	friend class MinlpBCP;
	friend class ColumnGenerator;
	friend class RMPManager;
	friend class LinearRelax;
	friend class LinearRelaxSolverGeneral;
	friend class LinearRelaxSolverMIP;
	friend class RelaxationSolver;
	friend class LagHeu;
	friend class LagHeu1;
	friend class LagHeu_SimAnnealing;
	friend class LagHeu2;
	friend class LagHeu2b;

  private:
		/** A lower bound of the optimal value.
		*/
		double low_bound;

		bool update_subdiv_bound_called;

		/** The indices of the fixed binary variables for MinlpBCP, seperated by block.
		*/
		map<int, set<int> > bcp_fixed_var;

		/** The partition cuts from MinlpBCP, seperated by blocks.
		*/
		vector<list<Pointer<SepQcFunc> > > part_con;

		/** Box for this node.
		*/
		dvector lower, upper;

		/** List of fixed branching variables.
		    Variables, which were just subdivided and shouln't be subdivided in the near future again.
		*/
		set<int> fix_branch_var;

		/** Solution points from the lagrangian problems.
		    To use as starting points when solving the lag problems again.
		*/
		vector<dvector> lagprob_solutions;

		/** Point, used by feasible_linrelax to store the solution of (R[U]).
		    Especially, it stores the values of the fixed binaries from fixed_var for RoundPartHeu.
		*/
		dvector ref_point;

		/** The dual point to the refpoint.
		*/
		dvector dual_point;

		/** Solution point of the RMP */
		dvector yz_RMP;

		/** The points and columns of the RMP for each block */
		vector<list<list<ExtremePoint>::iterator> > i_ExtremePoints;
		/** Points to the last RMP_point (in MinlpBCP), we checked.
		*/
		vector<list<ExtremePoint>::iterator> i_ExtremePoints_limit;

	public:
		/** Constructor for a MinlpNode.
		    @param lower_ Lower bounds on the variables.
				@param upper_ Upper bounds on the variables.
	      @param already Variables, which are already fixed. Can be left out.
		*/
		MinlpNode(const dvector& lower_, const dvector& upper_/*, const map<int, set<int> >& already_fixed/*=map<int, set<int> >*/);

		/** Copy-Constructor.
		    @param node The MinlpNode to copy.
		*/
		MinlpNode(MinlpNode& node);

		/** Evaluates the key and the related block and index.
		    Computes the minimum distance of the variables in x to the middle of its box regarding the discrete unfixed variables.
		    And gives the block and index of a variable, where this minimum is attached.
		    The indices of fixed variables are taken from bcp_fixed_var.
		*/
		pair<double, pair<int,int> > bcp_rho(const dvector& x, const vector<ivector>& block, const vector<bool>& discr);

		/** Checks if a sub-vector is inside the partition set
		    @param point A point.
		    @param k The block number.
				@param block The block structure.
		    @return True, if the point is inside.
		*/
		bool inside_part_set(dvector &point, int k, const vector<ivector>& block);

		/** Checks if a set of points is outside the partition set.
		    @param points A set of points.
		    @return A pointer to a point, which is not outside, or points.end(), if the set is outside.
		*/
		set<SolCandidate>::const_iterator outside_part_set(const set<SolCandidate>& points);

		/** Evaluates the key of a node.
		*/
		double key() { return low_bound; }
		
		/** Evaluates the key of the node which is currently the lower bound.
		 * If all discrete variables are fixed, the key in increased by rtol.
		 */
		double key(int nr_discr_var);

};

//------------------------------------------

#endif
