// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef DECOMP_H
#define DECOMP_H

#include "standard.h"
#include "usermatrix.h"
#include "problem.h"

class DecompGraph {
	/** Prints the graph.
	    @param out The ostream to print to.
	    @param g The DecompGraph to print.
	    @return The ostream out.
	*/
  friend ostream& operator<<(ostream& out, const DecompGraph& g);

  public:
		class Edge;
		class Node {
			friend class DecompGraph;
			private:
				int set_component(int comp);
			public:
				map<int, set<Edge>::iterator> adj;
				int weight;
				int component;

				Node(int weight_=0)
				: weight(weight_), component(-1)
				{ }
		};

		class Edge {
			public:
				map<int, Node>::iterator node1, node2;
				mutable int weight;

				Edge(int weight_=0)
				: weight(weight_)
				{ }

				Edge(map<int, Node>::iterator node1_, map<int, Node>::iterator node2_, int weight_=0)
				: node1(node1_), node2(node2_), weight(weight_)
				{ }

				/** Gives neighbour.
				    @param index Needs to be the index one of the nodes, which are adjacent to this edge.
				    @return The neighbour of node with number index.
			  */
//				map<int, Node>::iterator neighbour(int index) const {
//					return (index==node1->first) ? node2 : node1;
//				}

				Node& neighbour(const Node& node) const {
					return (&node==&node1->second) ? node2->second : node1->second;
				}

				bool operator<(const Edge& e) const {
					return pair<int,int>(node1->first, node2->first)<pair<int,int>(e.node1->first, e.node2->first);
				}
		};

		map<int, Node> nodes;
		set<Edge> edges;

		int nrcomp;

    /** Standard-Constructor.
    	  Constructs an empty graph.
    */
    DecompGraph()
		: nrcomp(0)
		{ }

		DecompGraph(const SparsityInfo& si);

		virtual ~DecompGraph() { }

		map<int, Node>::iterator add_node(int index, int weight=0);

		set<Edge>::iterator add_edge(map<int, Node>::iterator node1, map<int, Node>::iterator node2, int weight=0);
		set<Edge>::iterator add_edge(int node1, int node2, int weight=0);

  	/** The number of nodes.
  	    @return Then number of nodes: nodes.size().
  	*/
  	int n() const { return nodes.size(); }

  	/** The number of edges.
  	    @return The number of edges: edges.size().
  	*/
  	int m() const { return edges.size(); }

		void compute_connected_components();

		void get_component_members(vector<list<int> >& members);

		void compute_partition(int nparts);
};

/** Represents a function, which relates only to a subset of the variables of another function.
*/
class SplitFunc: public Func {
	private:
		/** Pointer to original function.
		*/
		const Pointer<Func> orig;

		/** The indices of the variables of original function, this function relates on.
		*/
		ivector indices;

    /** A point, we can use to evaluate the original function.
    */
    mutable dvector point;

		Func::CurvatureType curv_type;

	public:
    /** Variables, we relate to, but we ignore at evaluation.
		*/
    vector<bool> ignore;

		/** Constructor.
		    @param orig_ The original function.
		    @param indices_ The set of indices, for which this function represents the original one.
		    @param point_ The point, where we take the other values from.
		    @param ignore_ Indices of variables (of this function, not the original one) which should be ignored.
		*/
		SplitFunc(const Pointer<Func>& orig_, const ivector& indices_, const dvector& point_, const vector<int>& ignore_=vector<int>(0))
		: Func(indices_.dim()), orig(orig_), indices(indices_), point(point_), curv_type(orig_->get_curvature()), ignore(indices_.dim(), false)
		{ assert(orig!=NULL);
			for (vector<int>::const_iterator it(ignore_.begin()); it!=ignore_.end(); it++) ignore[*it]=true;
		}

		SplitFunc(const Pointer<Func>& orig_, const ivector& indices_, const dvector& point_, const set<int>& ignore_, Pointer<SparsityInfo> si=NULL)
		: Func(indices_.dim()), orig(orig_), indices(indices_), point(point_), curv_type(orig_->get_curvature()), ignore(indices_.dim(), false)
		{ assert(orig!=NULL);
			for (set<int>::const_iterator it(ignore_.begin()); it!=ignore_.end(); ++it) ignore[*it]=true;
			sparsity=si;
		}

		/** Sets point to value of x for variables, we represent and don't ignore.
		    @param x The UserVector to take to values from.
		*/
		void set_point(const UserVector<double>& x) const {
			for (int i=0; i<indices.size(); i++)
				if (!ignore[i]) point[indices[i]]=x(i);
		}

		double eval(const UserVector<double>& x) const {
			set_point(x);
			return orig->eval(point);
		}

	
		void grad(UserVector<double>& g, const UserVector<double>& x) const {
			set_point(x);
			dvector biggrad(orig->dim());
			orig->grad(biggrad, point);
			g=biggrad(indices);
			for (int i=0; i<indices.size(); i++)
				if (ignore[i]) g[i]=0.;
		}

		using Func::grad;

		void HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
			set_point(z);
			
			dvector big_x(point);
			for (int i=0; i<indices.size(); i++)
				if (!ignore[i]) big_x.SetElement(indices[i], x(i));
			
			dvector bighess(orig->dim());
			orig->HessMult(bighess, big_x, point);
			y=bighess(indices);
		}

		using Func::HessMult;

#ifdef FILIB_AVAILABLE
			bool is_interval_compliant() const { return orig->is_interval_compliant(); }

			void set_point(IntervalVector& bigx, const IntervalVector& x) const {
				for (int i=0; i<indices.size(); i++)
					if (!ignore[i]) bigx[indices[i]]=x(i);
			}

			interval<double> eval(const IntervalVector& x) const {
				IntervalVector bigx(point);
				set_point(bigx, x);
				return orig->eval(bigx);
			}

			void grad(IntervalVector& g, const IntervalVector& x) const {
				IntervalVector bigx(point);
				set_point(bigx, x);
				IntervalVector biggrad(bigx.dim());
				orig->grad(biggrad, bigx);
				for (int i=0; i<indices.size(); i++)
					g[i] = ignore[i] ? interval<double>(0.) : biggrad(indices[i]);
			}

			int valgrad(interval<double>& val, IntervalVector& y, const IntervalVector& x) const {
				IntervalVector bigx(point);
				set_point(bigx, x);
				IntervalVector biggrad(bigx.dim());
				int ret=orig->valgrad(val, biggrad, bigx);
				for (int i=0; i<indices.size(); i++)
					y[i] = ignore[i] ? interval<double>(0.) : biggrad(indices[i]);
				return ret;
			}

		using Func::valgrad;
#endif

		virtual void set_curvature(CurvatureType ct) { curv_type=ct; };
		virtual CurvatureType get_curvature() const  { return curv_type; };

		void print(ostream& out) const;
};


class SplittingScheme2 {
public:
	/** Gives position (or positions, if variable was copied) of an old variable in the new structure.
	    new_pos[i] gives a list of (blocknr, indexnr)-pairs, which are the same variables as variable i was in the old structure.
	    So it is the reverse of old_pos.
	*/
	vector<list<pair<int, int> > > new_pos;

	/** The new block structure.
	*/
	vector<ivector> newblock;

	SplittingScheme2(const DecompGraph& graph);
};

/** Decompose a function or problem.
    @class Decomposition
    Constructs also sparsity pattern and computes nonlinear and nonquadratic variable indices.
    @param Decomposition replace quadratic
    %options 0 or 1
    %default 1
    If 1, (s-)functions with constant hessians are replaced by their quadratic approximation.
    @param decompose method
		%level 1
    %options mincut or smooth
    %default mincut
    Determines, which method should be used to decompose the sparsity graph.
    @param decompose amount of copyvars
    %options double $\geq 0$
    %default 0
    Determines the amount of copy variables, which are allowed, when we split a problem, using the mincut-algorithm.
		@param Decomposition average blocksize
		% options integer $\geq 1$
		% default 5
		% level 2
		The average blocksize, we should try to partition the sparsity graph into.
*/
class Decomposition {
	private:
	  /** Indicates, whether a (s-)function should be replaced, if it seems to be quadratic (constant hessian).
	  */
		bool replace_if_quadratic;

		/** The average block size, we should try to achieve.
		*/
		int avg_blocksize;

		/** Parameters.
		*/
		Pointer<Param> param;

	public:
		Pointer<set<int> > I_var;

		/** (Default-)Constructor.
		    @param replace_if_quadratic_ Should we replace a s-function, if it's Hessian is constant over some sample point? Default is false.
		*/
		Decomposition(Pointer<Param> param_=NULL)
		: param(param_), replace_if_quadratic(param_ ? param_->get_i("Decomposition replace quadratic", 1) : 1),
		  avg_blocksize(param_ ? param_->get_i("Decomposition average blocksize", 5) : 5)
		{ }

	
	void set_sparsity_pattern(MinlpProblem& prob, const vector<vector<dvector> >& sample_set);

	
	/** Decompose a problem.
	    Computes sparsity patterns and indices of quadratic and nonquadratic variables.
	    Decompose sparsity-graphs.
	    Constructs SplittingScheme.
	    Decompose problem according to this information.
	    @param ss A reference to a Pointer, where we can store the SplittingScheme in.
	    @param prob The problem to decompose.
	    @param sample_set Points to compute the hessian for to get sparsity pattern, needed for functions, where s[k]!=NULL.
	    @return The decomposed problem.
	*/
	Pointer<MinlpProblem> decompose(MinlpProblem& prob, vector<vector<dvector> >& sample_set);

	Pointer<MinlpProblem> decompose(MinlpProblem& prob, vector<Pointer<SplittingScheme2> >& ss);
	void decompose(SepQcFunc& f, int block_offset, const SepQcFunc& old_f, int k, const SplittingScheme2& ss, Pointer<dvector> primal);

	Pointer<SplittingScheme2> get_splittingscheme(MinlpProblem& prob, int blocknr);

	Pointer<SepQcFunc> decompose(Pointer<Func> old_f, const dvector& primal);
};

#endif
