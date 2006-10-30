// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef GRAPH_H
#define GRAPH_H
#include "standard.h"

template <typename NodeDataType=int, typename EdgeDataType=int, bool multi_edges=false> class Edge;

// -------------------------------- Node ------------------------------------

template <typename NodeDataType=int, typename EdgeDataType=int, bool multi_edges=false> struct NodeAdj;
template <typename NodeDataType, typename EdgeDataType> struct NodeAdj<NodeDataType,EdgeDataType,false> {
	typedef Edge<NodeDataType, EdgeDataType, false> EdgeType;
	typedef map<int, typename set<EdgeType>::iterator> adj_type;
};
template <typename NodeDataType, typename EdgeDataType> struct NodeAdj<NodeDataType,EdgeDataType,true> {
	typedef Edge<NodeDataType, EdgeDataType, true> EdgeType;
	typedef multimap<int, typename set<EdgeType>::iterator> adj_type;
};

template <typename NodeDataType=int, typename EdgeDataType=int, bool multi_edges=false> class Node {
		friend ostream& operator<<(ostream& out, const Node<NodeDataType,EdgeDataType,multi_edges>& node) {
			out << node.index << '(' << node.data << ')' << ':';
			for (typename Node<NodeDataType,EdgeDataType,multi_edges>::adj_type::const_iterator it(node.adj.begin()); it!=node.adj.end(); ++it)
				out << ' ' << it->first << '(' << it->second->data << ')';
			return out;
		}

	public:
		typedef typename NodeAdj<NodeDataType,EdgeDataType,multi_edges>::adj_type adj_type;

	private:
		unsigned int index;

		mutable adj_type adj;

		typedef Node<NodeDataType, EdgeDataType, multi_edges> NodeType;
		typedef Edge<NodeDataType, EdgeDataType, multi_edges> EdgeType;

	public:
		mutable NodeDataType data;

		Node(unsigned int index_, const NodeDataType& data_=NodeDataType())
		: index(index_), data(data_)
		{ }

		unsigned int idx() const { return index; }

		bool operator<(const NodeType& node) const {
			return index<node.index;
		}

		void add_data(const NodeDataType& data2) const { data+=data2; }
		void add_neighbour(int neighbour_index, typename set<EdgeType>::iterator& edge) const {
			adj.insert(typename adj_type::value_type(neighbour_index, edge));
		}

		const adj_type& get_adj() const { return adj; }
};

// -------------------------------- Edge ------------------------------------

template <typename NodeDataType, typename EdgeDataType> class Edge<NodeDataType, EdgeDataType, false> {
	private:
		typedef Node<NodeDataType, EdgeDataType, false> NodeType;
		typedef Edge<NodeDataType, EdgeDataType, false> EdgeType;

		typename set<NodeType>::iterator node1, node2;

	public:
		mutable EdgeDataType data;

		Edge(const typename set<NodeType>::iterator& node1_, const typename set<NodeType>::iterator& node2_, const EdgeDataType& data_=EdgeDataType())
		: node1(node1_), node2(node2_), data(data_)
		{ }

		void add_data(const EdgeDataType& data2) const { data+=data2; }

		bool operator<(const EdgeType& e) const {
			return pair<unsigned int,unsigned int>(node1->idx(), node2->idx())<pair<unsigned int, unsigned int>(e.node1->idx(), e.node2->idx());
		}

		const NodeType& get_node1() const { return *node1; }
		const NodeType& get_node2() const { return *node2; }

/*
				Node& neighbour(const Node& node) const {
					return (&node==&node1->second) ? node2->second : node1->second;
				}
*/

};

// -------------------------------- MultiEdge ------------------------------------

template <typename NodeDataType, typename EdgeDataType> class Edge<NodeDataType, EdgeDataType, true> {
	private:
		typedef Node<NodeDataType, EdgeDataType, true> NodeType;
		typedef Edge<NodeDataType, EdgeDataType, true> EdgeType;

		typename set<NodeType>::iterator node1, node2;

	public:
		EdgeDataType data;

		Edge(const typename set<NodeType>::iterator& node1_, const typename set<NodeType>::iterator& node2_, const EdgeDataType& data_=EdgeDataType())
		: node1(node1_), node2(node2_), data(data_)
		{ }

		void add_data(const EdgeDataType& data2) const { data+=data2; }

		bool operator<(const EdgeType& e) const {
			if (node1->idx()==e.node1->idx() && node2->idx()==e.node2->idx()) return data<e.data;
			return pair<unsigned int,unsigned int>(node1->idx(), node2->idx())<pair<unsigned int, unsigned int>(e.node1->idx(), e.node2->idx());
		}

		const NodeType& get_node1() const { return *node1; }
		const NodeType& get_node2() const { return *node2; }
};

// -------------------------------- Graph ---------------------------------------

template <typename NodeDataType=int, typename EdgeDataType=int, bool directed=false, bool multi_edges=false> class Graph {
		friend ostream& operator<<(ostream& out, const Graph<NodeDataType,EdgeDataType,directed,multi_edges>& g) {
			out << g.n() << " nodes, " << g.m() << " edges:" << endl;
			for (typename set<Node<NodeDataType, EdgeDataType, multi_edges> >::const_iterator it(g.nodes.begin()); it!=g.nodes.end(); ++it)
				out << *it << endl;
			return out;
		}

	public:
		typedef Node<NodeDataType, EdgeDataType, multi_edges> NodeType;
		typedef Edge<NodeDataType, EdgeDataType, multi_edges> EdgeType;

		set<NodeType> nodes;
		set<EdgeType> edges;

		typename set<NodeType>::iterator add_node(unsigned int index, const NodeDataType& data=NodeDataType());
		const NodeType& get_node(unsigned int index) const;

		typename set<EdgeType>::iterator add_edge(typename set<NodeType>::iterator node1, typename set<NodeType>::iterator node2, const EdgeDataType& data=EdgeDataType());
		typename set<EdgeType>::iterator add_edge(unsigned int index1, unsigned int index2, const EdgeDataType& data=EdgeDataType());

  	/** The number of nodes.
  	    @return Then number of nodes: nodes.size().
  	*/
  	int n() const { return nodes.size(); }

  	/** The number of edges.
  	    @return The number of edges: edges.size().
  	*/
  	int m() const { return edges.size(); }

		void clear() { nodes.clear(); edges.clear(); }

};

template<typename NodeDataType, typename EdgeDataType, bool directed, bool multi_edges>
typename set<Node<NodeDataType, EdgeDataType, multi_edges> >::iterator Graph<NodeDataType,EdgeDataType,directed,multi_edges>::add_node(unsigned int index, const NodeDataType& data) {
	pair<typename set<NodeType>::iterator, bool> ret(nodes.insert(NodeType(index, data)));
	if (!ret.second) ret.first->add_data(data);
	return ret.first;
}

template<typename NodeDataType, typename EdgeDataType, bool directed, bool multi_edges>
const Node<NodeDataType,EdgeDataType,multi_edges>& Graph<NodeDataType,EdgeDataType,directed,multi_edges>::get_node(unsigned int index) const {
	typename set<NodeType>::const_iterator it(nodes.find(NodeType(index)));
	assert(it!=nodes.end());
	assert(index==it->idx());
	return *it;
}

template<typename NodeDataType, typename EdgeDataType, bool directed, bool multi_edges>
typename set<Edge<NodeDataType, EdgeDataType, multi_edges> >::iterator Graph<NodeDataType,EdgeDataType,directed,multi_edges>::add_edge(typename set<NodeType>::iterator node1, typename set<NodeType>::iterator node2, const EdgeDataType& data) {
	assert(node1->idx()!=node2->idx()); // do not allow loops
	if ((!directed) && node1->idx()>node2->idx()) return add_edge(node2, node1, data);

	pair<typename set<EdgeType>::iterator, bool> ret(edges.insert(EdgeType(node1, node2, data)));
	if (!ret.second) { // edge already there
		ret.first->add_data(data);
		return ret.first;
	}

	node1->add_neighbour(node2->idx(), ret.first);
	if (!directed) node2->add_neighbour(node1->idx(), ret.first);
	return ret.first;
}

template<typename NodeDataType, typename EdgeDataType, bool directed, bool multi_edges>
typename set<Edge<NodeDataType, EdgeDataType, multi_edges> >::iterator Graph<NodeDataType,EdgeDataType,directed,multi_edges>::add_edge(unsigned int index1, unsigned int index2, const EdgeDataType& data) {
	typename set<NodeType>::iterator node1(add_node(index1));
	typename set<NodeType>::iterator node2(add_node(index2));
	return add_edge(node1, node2, data);
}

#endif // GRAPH_H
