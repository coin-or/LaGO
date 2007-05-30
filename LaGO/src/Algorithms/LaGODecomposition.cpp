// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGODecomposition.hpp"
#include "LaGOSampling.hpp"
#include "LaGORestrictedFunction.hpp"
#include "LaGOSymSparseMatrix.hpp"

namespace LaGO {

void Decomposition::decompose() {
	decompose(data.obj);
	for (int c=0; c<data.numConstraints(); ++c)
		decompose(data.con[c]); 
}

void Decomposition::decompose(MINLPData::ObjCon& objcon) {
	objcon.decompfuncConstant=objcon.origfuncConstant;
	if (objcon.isLinear()) { // easy case
		objcon.decompfuncLin=objcon.origfuncLin;
		if (IsNull(objcon.sparsitygraph))
			objcon.sparsitygraph=new SparsityGraph(0,0);
		return;
	}
	
	// if the function does not know its sparsity pattern, we assume a dense function
	vector<int> nonzeros_dummy;
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
	
	// compute sparsity graph
	list<int> lin_nonzeros; // list of variables in nonzeros which turn out to be linear
	if (IsNull(objcon.sparsitygraph)) {
		sampling.monteCarlo(samplepoints, data.start_points.front(), nonzeros, lower, upper, 20); 
		computeSparsityGraph(objcon, lin_nonzeros, samplepoints, nonzeros);
	} else {
		cerr << "Not implemented yet: Filling of lin_nonzeros = nonzeros - nodes in graph" << endl;
		exit(EXIT_FAILURE);
	}
	SparsityGraph& graph(*objcon.sparsitygraph);
//	clog << "Sparsity graph for constraint " << objcon.name << ": " << endl << *objcon.sparsitygraph;

	// find connected components
	vector<bool> component_isnonquad;
	findConnectedComponents(component_isnonquad, graph);
	// mark all variables in nonquad. components as nonquadratic
	bool have_quadratic_component=false;
	for (SparsityGraph::iterator it_node(graph.begin()); it_node!=graph.end(); ++it_node) {
		const SparsityGraphNode& node(**it_node);
		if (node.isquad)
			if (component_isnonquad[node.component])
				const_cast<SparsityGraphNode&>(node).isquad=false;
			else
				have_quadratic_component=true;
//		clog << node << endl;
	}

	// decompose function
	createDecomposedFunctions(objcon, samplepoints.front(), nonzeros, lin_nonzeros, component_isnonquad, have_quadratic_component); 
}


void Decomposition::computeSparsityGraph(MINLPData::ObjCon& objcon, list<int>& lin_nonzeros, list<DenseVector>& samplepoints, const vector<int>& nonzeros) {
	SmartPtr<SparsityGraph> graph=new SparsityGraph(nonzeros.size(), 2*nonzeros.size());
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
	
	if (samplepoints.empty() || ++samplepoints.begin()==samplepoints.end()) {
		cerr << "Not enough sample points to generate a reliable sparsity graph." << endl;
		exit(EXIT_FAILURE);   		
	}

	// generate list of variables that are nonzero	
	if (graph->size()<nonzeros.size())
		for (unsigned int i=0; i<nodeits.size(); ++i)
			if (nodeits[i]==graph->end()) lin_nonzeros.push_back(i);
	
	objcon.sparsitygraph=graph;
}

void Decomposition::findConnectedComponents(vector<bool>& component_isnonquad, SparsityGraph& graph) {
	int nrcomp=0;
	SparsityGraph::iterator it_n(graph.begin());
	while (it_n!=graph.end()) {
		if((**it_n).component>=0) { ++it_n; continue; }
		component_isnonquad.push_back(setComponent(*it_n, nrcomp));
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

void Decomposition::createDecomposedFunctions(MINLPData::ObjCon& objcon, DenseVector& refpoint, const vector<int>& nonzeros, const list<int>& lin_nonzeros, const vector<bool>& component_isnonquad, bool have_quadratic_component) {
	assert(IsValid(objcon.sparsitygraph));
	SparsityGraph& graph(*objcon.sparsitygraph);
	int nr_components=component_isnonquad.size();

	objcon.decompfuncConstant=objcon.origfuncConstant;

	SparseVectorCreator decompfuncLin;
	if (IsValid(objcon.origfuncLin)) decompfuncLin.add(*objcon.origfuncLin);

	// setup reference point: 0 for linear and quadratic variables 
	for (list<int>::const_iterator it(lin_nonzeros.begin()); it!=lin_nonzeros.end(); ++it)
		refpoint[nonzeros[*it]]=0.;
	if (have_quadratic_component)
		for (SparsityGraph::iterator it_node(graph.begin()); it_node!=graph.end(); ++it_node) {
			const SparsityGraphNode& node(**it_node);
			if (node.isquad) refpoint[node.varindex]=0.;
		}

	// coefficients of linear part for linear and quadratic variables
	DenseVector grad;
	if (have_quadratic_component || !lin_nonzeros.empty()) {
		grad.resize(refpoint.size());
		objcon.origfuncNL->gradient(grad, refpoint);
		
		for (list<int>::const_iterator it(lin_nonzeros.begin()); it!=lin_nonzeros.end(); ++it)
			decompfuncLin[nonzeros[*it]]+=grad[nonzeros[*it]];
	}
	
	objcon.decompfuncNL.reserve(nr_components);
	objcon.decompmapping.resize(data.numVariables());
	
	int nr_nonquad_components=0;
	for (int comp=0; comp<nr_components; ++comp) {
		SmartPtr<BlockFunction> blockfunc=new BlockFunction();
		objcon.decompfuncNL.push_back(blockfunc);		

		// sparsity graph for block
		blockfunc->sparsitygraph=graph.getComponent(comp);
		SparsityGraph& comp_graph(*blockfunc->sparsitygraph);

		// indices for block and linear coefficients for quad. variables
		blockfunc->indices.reserve(comp_graph.size());
		for (SparsityGraph::iterator it(comp_graph.begin()); it!=comp_graph.end(); ++it) {
			const SparsityGraphNode& node(**it);
			blockfunc->indices.push_back(node.varindex);
			if (!component_isnonquad[comp]) decompfuncLin[node.varindex]+=grad[node.varindex];						

			objcon.decompmapping[node.varindex].push_back(pair<int,int>(comp, blockfunc->indices.size()-1));
			const_cast<SparsityGraphNode&>(node).varindex=blockfunc->indices.size()-1;
		}
		
		// function for block
		if (component_isnonquad[comp]) { // nonquadratic
			blockfunc->nonquad=new RestrictedFunction(GetRawPtr(objcon.origfuncNL), blockfunc->indices, refpoint);
			++nr_nonquad_components;
						
		}	else { // quadratic
			SymSparseMatrixCreator creator(blockfunc->indices.size());
			for (SparsityGraph::arc_iterator it_arc(comp_graph.arc_begin()); it_arc!=comp_graph.arc_end(); ++it_arc) {
				const SparsityGraphNode& node1(**(*it_arc).head());
				const SparsityGraphNode& node2(**(*it_arc).tail());
				double coeff=(**it_arc).hessian_entry/2.;
				if (coeff) 
					creator.insert(node1.varindex, node2.varindex, coeff);
			}
			blockfunc->quad=new SymSparseMatrix(creator);
		}
	}
	objcon.decompfuncLin=new SparseVector(decompfuncLin);
	
	// correct constant term
	if (nr_nonquad_components!=1) {
		double val=objcon.origfuncNL->eval(refpoint);
		if (nr_nonquad_components==0) { // function is quadratic
			objcon.decompfuncConstant+=val;
		} else { // function has more then one nonquadratic block
			objcon.decompfuncConstant+=(nr_nonquad_components-1)*val;
		}
	}
	
}
	
} // namespace LaGO
