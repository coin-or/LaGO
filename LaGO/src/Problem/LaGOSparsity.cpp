// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSparsity.hpp"

namespace LaGO {

SparsityGraph::SparsityGraph(SparsityGraph& graph, vector<int>& indices_map)
: Cgc::DynNet<SparsityGraphNode, SparsityGraphEdge>(graph.size(), graph.arc_size())
{
	add(graph, indices_map);
}

SmartPtr<SparsityGraph> SparsityGraph::getComponent(int comp) const {
	SparsityGraph* graph=new SparsityGraph(size(), 2*size());

	map<Cgc::NodeId, const_iterator> map_nodes;
	for (const_iterator it_node(begin()); it_node!=end(); ++it_node)
		if ((**it_node).component==comp)
			map_nodes[getNodeId(it_node)]=graph->insert(**it_node);
	
	for (const_arc_iterator it_arc(arc_begin()); it_arc!=arc_end(); ++it_arc) {
		if ((**(*it_arc).tail()).component==comp)
			graph->arc_insert(
				map_nodes[getNodeId((*it_arc).tail())],
				**it_arc,
				map_nodes[getNodeId((*it_arc).head())]
			);		
	}

	return graph;
}

void SparsityGraph::add(SparsityGraph& graph, vector<int>& indices_map) {
//	clog << "graph before: " << *this;
//	clog << "graph adding: " << graph;
	bool allnew=(size()==0); // do not need to check for existing nodes or arcs if graph is empty
	map<Cgc::NodeId, const_iterator> map_nodes;
	for (const_iterator it_node(graph.begin()); it_node!=graph.end(); ++it_node) {
		SparsityGraphNode newnode(**it_node);
		newnode.varindex=indices_map[newnode.varindex];
		if (allnew) {
			map_nodes[getNodeId(it_node)]=insert(newnode);
		}	else {
			iterator it=find(newnode);
			map_nodes[getNodeId(it_node)]=it==end() ? insert(newnode) : it;
		}
	}
	
	for (const_arc_iterator it_arc(graph.arc_begin()); it_arc!=graph.arc_end(); ++it_arc) {
		const_iterator tail=map_nodes[getNodeId((*it_arc).tail())];
		const_iterator head=map_nodes[getNodeId((*it_arc).head())];
//		clog << "wanna add arc " << (**(*it_arc).tail()).varindex << ',' << (**(*it_arc).head()).varindex << endl;  
		if (allnew || (arc_find(tail, head)==arc_end())) {
			arc_insert(tail, **it_arc, head);
//			clog << "added arc " << (**tail).varindex << ',' << (**head).varindex << endl; 
		} else {
//			clog << "skiped arc " << (**tail).varindex << ',' << (**head).varindex << endl;
		}
	}
//	clog << "graph after: " << *this;
}

} // namespace LaGO
