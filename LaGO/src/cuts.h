// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef CUTS_H
#define CUTS_H

#include "standard.h"
#include "opt.h"
#include "MINLPData.h"

/** A class for storing the information that is neccessary to generate an interval gradient cut.
    The cut has the form b_x * x + b_w * w + b_z * z + c <= 0, w-z=x-ref_x, w,z>=0, w<=ref_x-\lb x, z<=\ub x-ref_x, (w+z)_i<=max(\ub x_i-ref_x_i, ref_x_i-\lb x_i)
*/
class IntervalGradientCut {
	public:
		class ConstraintInfo {
			public:
				/** A constant */
				double c;
				/** Linear parts, not scaled by .5 !*/
				Pointer<UserVector<double> > b_x;
				Pointer<UserVector<double> > b_w;
				Pointer<UserVector<double> > b_z;

				ConstraintInfo(const Pointer<UserVector<double> >& b_x_, const Pointer<UserVector<double> >& b_w_, const Pointer<UserVector<double> >& b_z_, double c_)
				: b_x(b_x_), b_w(b_w_), b_z(b_z_), c(c_)
				{ }
		};

		/** Information for the cuts from the several constraints.
		*/
		list<ConstraintInfo> coninfos;
		/** Length of coninfos list.
		*/
		int coninfos_size;

		/** The reference point.
		*/
		dvector ref_x;
		/** The upper bounds of w_i and z_i.
		*/
		dvector w_z_up;
		/** The indices in the block of the variables in the (condensed) w and z.
		*/
		ivector indices;

		IntervalGradientCut(const dvector& ref_x_, const dvector& w_z_up_, const ivector& indices_)
		: ref_x(ref_x_), w_z_up(w_z_up_), indices(indices_), coninfos_size(0)
		{ }
};


class SimpleCut {
	friend ostream& operator<<(ostream& out, const SimpleCut& cut);
	public:
		Pointer<UserVector<double> > coeff;
		double constant;
		
		SimpleCut()
		: constant(0.)
		{ }
	
		SimpleCut(const Pointer<UserVector<double> >& coeff_, double constant_)
		: coeff(coeff_), constant(constant_)
		{ }
		
		bool feasible(const dvector& x, double tol) const {
			return 2*(*coeff*x)+constant<tol;
		}
		
		/** Scales the coefficients (and constant) of the generates cut such that the inf norm of the coefficients does not exceed max_coeff
		 */
		void scale(double max_coeff=100.);
};

class LinearizationCut : public SimpleCut {
	public:
		LinearizationCut()
		: objcon_nr(-2), derived_from_lower_part(true)
		{ }

		LinearizationCut(const Pointer<UserVector<double> >& coeff_, double constant_)
		: SimpleCut(coeff_, constant_), objcon_nr(-2), derived_from_lower_part(true)
		{ }

		LinearizationCut(const Pointer<UserVector<double> >& coeff_, double constant_, int objcon_nr_, bool derived_from_lower_part_, const dvector& x_, double t_)
		: SimpleCut(coeff_, constant_), objcon_nr(objcon_nr_), derived_from_lower_part(derived_from_lower_part_), x(x_), t_value(t_)
		{ }
		
		/** The number of the constraint from which the cut was derived. -1 if it was derived from the objective function. -2 if not set.
		 */
		int objcon_nr;
		/** Indicates whether the cut was derived from the lower or upper part of the constraint.
		 */
		bool derived_from_lower_part;
		/** The point where the cut was generated (as a block in the non reformulated problem).
		 */
		dvector x;
		/** Value of t-variable for this block, if reformulated.
		 */
		double t_value;
		
};


/** Constructs a interval gradient cut.
*/
class IntervalGradientCutGenerator {
	private:
		/** The problem, we generate the cuts for.
		*/
		Pointer<MinlpProblem> prob;

		vector<Pointer<SparsityInfo> > sparsity;

	public:
		IntervalGradientCutGenerator(Pointer<MinlpProblem> prob_=NULL)
		: prob(prob_)
		{ prob->get_sparsity(sparsity); }

		void set_problem(Pointer<MinlpProblem> prob_) { prob=prob_; prob->get_sparsity(sparsity); }

		/** Generate interval gradient cut information according to a solution estimate and a block.
		    @param x The solution estimate. (reference point)
		    @param k The block-number.
				@param low Current lower bounds.
				@param up Current upper bounds.
		*/
		Pointer<IntervalGradientCut> get_cuts(const dvector& x, int k, const dvector& low, const dvector& up);

		Pointer<IntervalGradientCut> update_cuts(Pointer<IntervalGradientCut> cut, int k, const dvector& low, const dvector& up) {
			return get_cuts(cut->ref_x, k, low, up);
		}
};

class CutPool;
class MinlpNode;

/** Represents a Cut in the CutPool.
*/
template <class CutType> class Cut {
	friend class CutPool;
//	friend ostream& operator<<(ostream& out, Cut<CutType>& cut);

	class NodeInfo {
		public:
			/** The time (as number or (R[U]) solves), the cut was inactive for this node.
			*/
			int inactive_time;

			NodeInfo()
			: inactive_time(0)
			{ }

			NodeInfo(const NodeInfo& cni)
			: inactive_time(cni.inactive_time)
			{ }
	};

	private:
		/** The cut.
		*/
		Pointer<CutType> cut;

		/** The nodes, which use this cut.
		*/
		map<Pointer<MinlpNode>, NodeInfo> nodes;

		/** Indicates, whether this is a global cut (valid for all nodes).
		*/
		bool global;
		/** If global, the time, this cut was inactive now.
		*/
		int inactive_time_global;

		/** Used to distinguish between cuts, which were added recently, and cuts, which existed "from the beginning".
		*/
		bool tagged;

	public:
		Cut(Pointer<CutType> cut_, Pointer<MinlpNode> node=NULL)
		: cut(cut_), global(false), tagged(false), inactive_time_global(0)
		{ if (node) add_node(node);
			else global=true;
		}

		Cut(const Cut<CutType>& cut_)
		: cut(cut_.cut), nodes(cut_.nodes), global(cut_.global), tagged(cut_.tagged), inactive_time_global(cut_.inactive_time_global)
		{ }

		virtual ~Cut();

		bool valid(const Pointer<MinlpNode>& node=NULL) const {
			return global || (node && nodes.count(node));
		}

		void add_node(Pointer<MinlpNode> node) {
			nodes.insert(pair<Pointer<MinlpNode>, NodeInfo>(node, NodeInfo()));
		}

		/** Duplicates a CutNodeInfo and uses it for a new node.
		    Needed, when a node is splitted into new nodes.
		*/
		void duplicate_nodeinfo(Pointer<MinlpNode> oldnode, Pointer<MinlpNode> newnode);

		/** Removed a node from the nodes-list, if it exists there.
		    @return True, if the cut was local and this was the last node, which used this cut. False, else.
		*/
		bool remove_node(Pointer<MinlpNode> node);

		/** Increases or resets the inactivity time for a node or global.
		    If this time exceeds a limit, the node is removed from the nodes set.
				For a global cut, the value of node is nonrelevant.
				@param increase Whether to increase the inactive timer or to set it to 0.
				@param limit The limit for inactivity.
				@return If the cut isn't used by any node anymore or was global and the time limit is exceeded, the first return value is true. False otherwise. If the cut is invalid for the node now, true is returned in second return value. False otherwise.
		*/
		pair<bool, bool> set_inactivetime(Pointer<MinlpNode> node, bool increase, int limit);

		const CutType& get_cut() const { return *cut; }
};

/** A Cut management system.
    @class CutPool
		@param Cut inactive time limit
		%options integer
		%defaults 10
		%level 2
		Determines, after how many consecutive (R)-solves, in which a cut was inactive, will be removed. 0 means no removals.
		Affects only cuts which are not specific for a node (global). So called local cuts are removed (from their node) as soon as they become inactive.
*/
class CutPool {
	private:
		/** Cuts which are defined by a seperating hyperplane.
		    Given a vector b and a double c, the cut has the form 2*b^T x + c <= 0.
		*/
		list<Cut<SimpleCut> > simplecuts;

		/** Cuts also defined by a seperating hyperplane. The hyperplane was derived by linearization of a (convexified) constraint or objective function.
		 */
		list<Cut<LinearizationCut> > linearizationcuts;

		/** IntervalGradientCuts.
		*/
		list<Cut<IntervalGradientCut> > intgradcuts;

		/** The number of cuts in this CutPool.
		*/
		int cuts_size;

		int inactivetime_limit_global;
		int inactivetime_limit_local;

	public:
		typedef enum { SIMPLE, LINEARIZATION, INTERVALGRADIENT } CutType;

		/** Class for communication between CutPool and LinearRelax.
		    Adding a cut, returnes a CutInfo, so the linear relaxation can remember a cut, which is added to it.
				Solving the linear relaxation sets the inactive flag. This is checked by update_nodeinfo, which sets the removed flag.
				This flag can be checked by LinearRelax again.
		*/
		class CutInfo {
			friend class CutPool;
			public:
				/** A pointer to the cut, if it's a simple cut.
				*/
				list<Cut<SimpleCut> >::iterator it_simplecuts;

				/** A pointer to the cut, if it's a IntervalGradientCut.
				*/
				list<Cut<IntervalGradientCut> >::iterator it_intgradcuts;
				
				/** A pointer to the cut, if it's a LinearizationCut.
				*/
				list<Cut<LinearizationCut> >::iterator it_linearizationcuts;

				/** The kind of cut.
				*/
				CutType type;

				int block_nr;
				/** Set by LinearRelax, if the cut was inactive.
				*/
				bool inactive;
				/** Set by CutPool, if the cut was removed or the node was removed for the cut.
				*/
				bool removed;

				/** Set by LinearRelaxSolver, to remember the constraints in the MIPSolver, which were added for this cut.
				*/
				list<const MIPSolver::RowItem*> rowitems;

				/** Set by LinearRelaxSolver, to remember the variables in the MIPSolver, which were added for this cut.
				*/
				list<const MIPSolver::ColItem*> colitems;

				CutInfo(list<Cut<SimpleCut> >::iterator& it, int bnr)
				: it_simplecuts(it), type(SIMPLE), block_nr(bnr), inactive(false), removed(false)
				{ }

				CutInfo(list<Cut<LinearizationCut> >::iterator& it, int bnr)
				: it_linearizationcuts(it), type(LINEARIZATION), block_nr(bnr), inactive(false), removed(false)
				{ }

				CutInfo(list<Cut<IntervalGradientCut> >::iterator& it_intgradcuts_, int bnr)
				: it_intgradcuts(it_intgradcuts_), type(INTERVALGRADIENT), block_nr(bnr), inactive(false), removed(false)
				{ }

				pair<bool, bool> set_inactivetime(Pointer<MinlpNode> node, int limit);
		};

		CutPool(int inactivetime_limit_global_=10, int inactivetime_limit_local_=3)
		: inactivetime_limit_global(inactivetime_limit_global_), inactivetime_limit_local(inactivetime_limit_local_), cuts_size(0)
		{ }

		/** A copy of a given cutpool, containing only these cuts, which are defined for a specific node, or are global.
		*/
		CutPool(const CutPool& cutpool, Pointer<MinlpNode> node=NULL);
		
		/** Adds a cut to the pool.
		    @param cut The cut to add.
				@param node For local cuts, put the node, the cut belongs to, here. For global cuts, put NULL here.
				@param blocknr The block number to put into the CutInfo. -1 for cuts that couple several blocks.
		    @return cutinfo A CutInfo for the new cut.
		*/
		CutInfo add_cut(Pointer<SimpleCut>, Pointer<MinlpNode> node, int blocknr=-1);

		/** Adds a IntervalGradientCut to the pool.
		    @param cut The cut to add.
				@param node For local cuts, put the node, the cut belongs to, here. For global cuts, put NULL here.
				@param blocknr The block number to put in the CutInfo.
		    @return cutinfo A CutInfo for the new cut.
		*/
		CutInfo add_cut(Pointer<IntervalGradientCut> intgradcut, Pointer<MinlpNode> node, int blocknr);

		/** Adds a LinearizationCut to the pool.
		 */
		CutInfo add_cut(Pointer<LinearizationCut> linearizationcut, Pointer<MinlpNode> node, int blocknr);

		/** Adds the untagged cuts of a cutpool to this one.
		    @param cutpool The cuts to add.
				@param node The node, which should use these cuts. Or NULL, if intended as global cuts.
		*/
		void integrate(const CutPool& cutpool, Pointer<MinlpNode> node=NULL);

		void duplicate_nodeinfo(Pointer<MinlpNode> oldnode, Pointer<MinlpNode> newnode);

		void remove_node(Pointer<MinlpNode> node);

		/** Updates the information, which cuts were inactive.
		    If a cut is not needed anymore, it's removed from the cutlist.
				@param cutinfos On input, info about the inactivity of the cuts. On output, info about which cuts are invalid for the node now.
				@param block_nr Only CutInfo's with this block number are checked.
		*/
		void update_nodeinfo(list<CutInfo>& cutinfos, int block_nr, Pointer<MinlpNode> node);

		/** Builds a list of CutInfo's to point to these cuts, which belong to a given node.
		    This method is used for the communication between the CutPool and the LinearRelaxSolver about inactivity of cuts.
				@cutinfos A list, where we can push back the CutInfos.
				@blocknr The number of the block to put into the CutInfos. Give -1 for coupling cuts.
		*/
		void get_cuts(list<CutInfo>& cutinfos, int blocknr, Pointer<MinlpNode> node=NULL);

		/** Updates the IntervalGradientCuts of a specific node after the bounds changed.
		*/
		void update_cuts(Pointer<MinlpNode> node, int blocknr, const dvector& low, const dvector& up, IntervalGradientCutGenerator& generator, LinearizedConCutGenerator& linconcutgen);

		/** Checks, whether a point is feasible for the cuts.
		    Checks each global cut and the local cuts, which belong to the given node (if not NULL), if the cut is valid for the constraints.
		*/
		bool feasible(Pointer<MinlpNode> node, const dvector& x, double tol=1E-4) const;

		/** The number of local cuts for a specific node.
		    Does not count global cuts.
		*/
		int nr_local_cuts(Pointer<MinlpNode> node) const;
		/** The number of global cuts.
		    Does not count local cuts.
		*/
		int nr_global_cuts() const;
		/** The complete number of cuts, all local cuts plus the global ones.
		*/
		int nr_all_cuts() const { return cuts_size; }
};

class LinearizedConCutGenerator {
	private:
		Pointer<MinlpProblem> prob;
		
		Pointer<MINLPData> minlpdata;
		
		Pointer<Reformulation> reform;
		
		static const double tol;

	public:
		/** The maximum violation of a quadratic underestimator in a reference point.
		*/
		double max_violation;	
	
		LinearizedConCutGenerator(Pointer<MinlpProblem> prob_, Pointer<MINLPData> minlpdata_=NULL, Pointer<Reformulation> reform_=NULL);

		/** Returns a linearization of a constraint, which is defined over one block only.
				@param x The point, in which we linearize the constraint. block[block_nr].size()
				@param c The constraint number.
				@param block_nr The number of the block, for which the constraint is defined.
				@param val The value of the constraint in x, if available. INFINITY else.
		*/
		LinearizationCut get_cut(const dvector& x, int c, int block_nr, double val=INFINITY);

		int get_cuts(list<pair<LinearizationCut, pair<int, bool> > >& cuts, const MINLPData::ObjCon& objcon, int objcon_nr, const dvector& x, const dvector& lower, const dvector& upper, bool violated_polyest_only=false);
		int get_cuts(list<pair<LinearizationCut, pair<int, bool> > >& cuts, const MINLPData::Constraint& objcon, int objcon_nr, const dvector& x, const dvector& lower, const dvector& upper, bool violated_polyest_only=false);
		/** Computes linearization cuts by linearization of the convex relaxation in a given point.
		 * @param cuts A list of type { cut , { blocknr, is the function where the cut was derived from convex? } } to store the cuts.
		 * @param x The point where to compute the cut.
		 * @param lower The current lower bound on the variables.
		 * @param upper The current upper bound on the variables.
		 * @param violated_polyest_only If true, cuts are only generated for functions where x violates the polynomial under/over estimator. Useful if cut is ment to cuf off a solution of the relaxation.
		 */
		int get_cuts(list<pair<LinearizationCut, pair<int, bool> > >& cuts, const dvector& x, const dvector& lower, const dvector& upper, bool violated_polyest_only=false);
		
		/** Updates a cut, using a current box.
		 * @param lower The lower bounds on the variables, here only the lower bounds for the considered block.
		 * @param upper The upper boudns on the variables, here only the upper bounds for the considered block.
		 */
		Pointer<LinearizationCut> update_cut(LinearizationCut& cut, unsigned int block_nr, const dvector& lower, const dvector& upper);
		
		void set_problem(Pointer<MinlpProblem> prob_) { prob=prob_; }
		void set_reform(Pointer<Reformulation> reform_);
		void set_MINLPData(Pointer<MINLPData> minlpdata_) { minlpdata=minlpdata_; }
		
};

class LinearRelax;

#endif
