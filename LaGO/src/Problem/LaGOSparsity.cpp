// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSparsity.hpp"

namespace LaGO {

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

} // namespace LaGO
