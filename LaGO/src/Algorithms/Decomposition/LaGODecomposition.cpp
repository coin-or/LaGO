// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGODecomposition.hpp"
#include "LaGOSampling.hpp"

// connected component finder from Cgc
#include "ConnComp.h"

namespace LaGO {

void Decomposition::decompose() {
	decompose(data.obj);
	for (int c=0; c<data.numConstraints(); ++c)
		decompose(data.con[c]); 
}

void Decomposition::decompose(MINLPData::ObjCon& objcon) {
	objcon.decompfuncConstant=objcon.origfuncConstant;
	if (IsNull(objcon.origfuncNL)) {
		objcon.decompfuncLin=objcon.origfuncLin;
		return;
	}
	if (IsValid(objcon.origfuncLin))
		objcon.decompfuncLin=new SparseVector(*objcon.origfuncLin);
	else
		objcon.decompfuncLin=new SparseVector();
	
	
	vector<int> nonzeros_dummy; // if the function does not know its sparsity pattern, we assume a dense function
	if (!objcon.origfuncNL->haveSparsity()) {
		nonzeros_dummy.reserve(data.numVariables());
		for (int i=0; i<data.numVariables(); ++i) nonzeros_dummy.push_back(i);
	}
	const vector<int>& nonzeros(objcon.origfuncNL->haveSparsity() ? objcon.origfuncNL->getSparsity() : nonzeros_dummy);
	 
	DenseVector lower(nonzeros.size());
	DenseVector upper(nonzeros.size());
	data.getBox(lower, upper, nonzeros);
	
	// generate sample points (w.r.t. origfuncNL->sparsity)
	list<DenseVector> samplepoints; 
	Sampling sampling;
	sampling.addVector(samplepoints, data.start_points);
	assert(!data.start_points.empty());
	sampling.monteCarlo(samplepoints, data.start_points.front(), nonzeros, lower, upper, 20); 
	
	// compute sparsity graph
	if (IsNull(objcon.sparsitygraph))
		computeSparsityGraph(objcon, samplepoints, nonzeros);

	// find connected components
	
	// distinguish between nonquadratic and quadratic components
	
	// decompose function

	
}


void Decomposition::computeSparsityGraph(MINLPData::ObjCon& objcon, list<DenseVector>& samplepoints, const vector<int>& nonzeros) {
	SmartPtr<SparsityGraph> graph=new SparsityGraph(nonzeros.size(),nonzeros.size());
	vector<SparsityGraph::iterator> nodeits(nonzeros.size(), graph->end());

	DenseVector e(data.numVariables()); // storage for hessian multiplier
	DenseVector hm(data.numVariables()); // storage for hessian-vector product
	for (unsigned int i=0; i<nonzeros.size(); ++i) {
		e[nonzeros[i]]=1.;
		list<DenseVector>::iterator it_sp(samplepoints.begin());
		while (it_sp!=samplepoints.end()) {
			try {
				objcon.origfuncNL->hessianVectorProduct(hm, *it_sp, e);
			} catch (FunctionEvaluationError error) {
				clog << "computeSparsityPattern for " << objcon.name << ": skip sample point due to " << error << endl;
				it_sp=samplepoints.erase(it_sp);
				continue;
			}
			
			for (unsigned int j=0; j<nonzeros.size(); ++j) {
				if (nonzeros[j]>nonzeros[i]) continue;
				if (hm[nonzeros[j]]==0) continue; // zero entry in hessian
				
				if (nodeits[i]==graph->end()) { // first time we find variable nonzeros[i] in the hessian
					nodeits[i]=graph->insert(SparsityGraphNode(nonzeros[i]));
				}
				if (nodeits[j]==graph->end()) { // first time we find variable nonzeros[j] in the hessian
					nodeits[j]=graph->insert(SparsityGraphNode(nonzeros[j]));
				}
//				clog << '(' << nonzeros[i] << ',' << nonzeros[j] << ") " << hm[nonzeros[j]] << endl;
				SparsityGraph::arc_iterator arc=graph->arc_find(nodeits[i], nodeits[j]);
				if (arc==graph->arc_end()) { // new arc
					graph->arc_insert(nodeits[i], SparsityGraphEdge(hm[nonzeros[j]]), nodeits[j]);
					if (i!=j) graph->arc_insert(nodeits[j], SparsityGraphEdge(0.), nodeits[i]);
				} else {
					if ((**arc).hessian_entry!=hm[nonzeros[j]]) {
						(const_cast<SparsityGraphNode*>(&**nodeits[i]))->isquad=false;
						(const_cast<SparsityGraphNode*>(&**nodeits[j]))->isquad=false;
					}					
				}				
			}
			++it_sp;
		} 
		e[nonzeros[i]]=0.;
	}
	
	clog << "Sparsity graph for constraint " << objcon.name << ": " << endl << *graph;
	
	// search for connected components
	findConnectedComponents(*graph);
	for (unsigned int i=0; i<nodeits.size(); ++i)
		if (nodeits[i]!=graph->end()) clog << **nodeits[i] << endl;
			
#if 0
	if (!quadratic->empty()) { // linear coefficients for quadratic variables
		dvector try_point(sample_set.front());
		for (map<int, QuadraticVariable>::iterator it(quadratic->begin()); it!=quadratic->end(); ++it)
			try_point[it->first]=0.;
		dvector grad(try_point.dim());
		f.grad(grad, try_point);
		for (map<int, QuadraticVariable>::iterator it(quadratic->begin()); it!=quadratic->end(); ++it)
			it->second.coeff_lin=grad(it->first);
	}
#endif	
	
	objcon.sparsitygraph=graph;
}


void Decomposition::findConnectedComponents(SparsityGraph& graph) {
	int nrcomp=0;
	SparsityGraph::iterator it_n(graph.begin());
	while (it_n!=graph.end()) {
		if((**it_n).component>=0) { ++it_n; continue; }
		bool isnonquad=setComponent(*it_n, nrcomp);
		
		if (isnonquad) {	// mark all variables in this component as nonquadratic; TODO: find a nicer way to do it 
			for(SparsityGraph::iterator it(it_n); it!=graph.end(); ++it)
				if ((**it).component==nrcomp && (**it).isquad)
					const_cast<SparsityGraphNode&>(**it).isquad=false;
		}
		
		++it_n;
		++nrcomp;
	}	
}

bool Decomposition::setComponent(const SparsityGraph::Node& node, int comp) {
	if ((*node).component==comp) return false; // already set
	const_cast<SparsityGraphNode&>(*node).component=comp;

	bool isnonquad=!(*node).isquad;
	for (SparsityGraph::Node::iterator it(node.begin()); it!=node.end(); ++it)		
		if (setComponent(*(*it).head(),comp)) isnonquad=true;

	return isnonquad;
}
	
} // namespace LaGO
