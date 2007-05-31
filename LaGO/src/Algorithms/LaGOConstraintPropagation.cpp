// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOConstraintPropagation.hpp"

namespace LaGO {
	
ostream& operator<<(ostream& out, const ConstraintPropagation::EdgeInfo& edge) {
	for (list<CoinTriple<int, ConstraintPropagation::BoundType, double> >::const_iterator it(edge.coninfo.begin()); it!=edge.coninfo.end(); ++it) {
		out << '(' << it->first << ',';
		switch (it->second) {
			case ConstraintPropagation::LOWER: out << "low"; break;
			case ConstraintPropagation::UPPER: out << "up"; break;
			case ConstraintPropagation::ANY: out << "both"; break;						
		}
		out << ',' << it->third << ')';
	}
	return out; 
}

int ConstraintPropagation::funceval_limit=10000;
double ConstraintPropagation::min_impr=0.01;

void ConstraintPropagation::initDependencyGraph() {
//	reduction_by_block.resize(prob->block.size());

	assert(depgraph.size()==0);

	vector<DependencyGraph::iterator> nodes(data.numVariables());
	for (int i=0; i<data.numVariables(); ++i)
		nodes[i]=depgraph.insert(NodeInfo(i));
	DependencyGraph::iterator dummy_node(depgraph.insert(NodeInfo(data.numVariables())));

	for (int c=0; c<data.numConstraints(); ++c) {
		const MINLPData::Constraint& con(data.getConstraint(c));

		const int* linvar=con.decompfuncLin->getIndices();
		const double* linvar_el=con.decompfuncLin->getElements();
		
		if (con.decompfuncLin->getNumElements()==1 && con.decompfuncNL.empty()) { // func. is of form a*x+const
			// add dummy edge from dummy node to the single variable
			DependencyGraph::arc_iterator arc=depgraph.arc_find(dummy_node, nodes[*linvar]);
			if (arc==depgraph.arc_end()) // new arc
				depgraph.arc_insert(dummy_node, EdgeInfo(c, ANY, *linvar_el), nodes[*linvar]);
			else
				(**arc).add(c, ANY, *linvar_el);
			continue;
		} 
		assert(con.decompfuncNL.empty() || (int)con.decompmapping.size()==data.numVariables());
		
		int i_linvar=0; // runs over linear variables in con.decompfuncLin
		for (int i=0; i<data.numVariables(); ++i) { 
			if (linvar[i_linvar]<i && i_linvar<con.decompfuncLin->getNumElements()-1) {
				++i_linvar;
				assert(linvar[i_linvar-1]<linvar[i_linvar]);
			}
			if (linvar[i_linvar]!=i && (con.decompfuncNL.empty() || con.decompmapping[i].empty())) continue; // variable not in function
			
			for (int i_linvar2=0; i_linvar2<con.decompfuncLin->getNumElements(); ++i_linvar2) {
				if (linvar[i_linvar2]==i) continue;
				
				BoundType wb;
				// if var i is nonlinear, I cannot say what effect a change on one of its bounds has on var i_linvar2 can have
				if ((!con.decompfuncNL.empty()) && !con.decompmapping[i].empty()) wb=ANY;
				else if (con.lower>-getInfinity() && con.upper<getInfinity()) wb=ANY;
				else if (con.lower>-getInfinity()) wb=(linvar_el[i_linvar]<0 ? LOWER : UPPER);
				else if (con.upper< getInfinity()) wb=(linvar_el[i_linvar]<0 ? UPPER : LOWER);

				DependencyGraph::arc_iterator arc=depgraph.arc_find(nodes[i], nodes[linvar[i_linvar2]]);
				if (arc==depgraph.arc_end()) // new arc
					depgraph.arc_insert(nodes[i], EdgeInfo(c, wb, linvar_el[i_linvar2]), nodes[linvar[i_linvar2]]);
				else
					(**arc).add(c, wb, linvar_el[i_linvar2]);
			}			
		}
	}
//	clog << depgraph;

	clog << "ConstraintPropagation dependency graph nr. edges: " << depgraph.arc_size() << " (" << dummy_node->size() << ')' << endl;
}

bool ConstraintPropagation::evalConstraint(interval<double>& val, const MINLPData::Constraint& con, IntervalVector& box, int index) {
	val=con.decompfuncConstant;
	try {
		for (unsigned int k=0; k<con.decompfuncNL.size(); ++k) {
			if (!con.decompfuncNL[k]->canIntervalEvaluation()) return false;
			IntervalVector block(box, con.decompfuncNL[k]->indices);
			val+=con.decompfuncNL[k]->eval(block);					
		}
	} catch (FunctionEvaluationError error) {
		cerr << "Error evaluating block of constraint " << con.name << " over an interval: " << error;
		return false;					
	} 
	box[index]=interval<double>::ZERO();
	val+=box* *con.decompfuncLin;
	return true;
}

double tiny_tol=1E-10;

BoxReductionStatistics ConstraintPropagation::run(DenseVector& newlow, DenseVector& newup, const DenseVector& oldlow, const DenseVector& oldup, set<pair<const DependencyGraph::Node*, BoundType> >& nodeset) {
	BoxReductionStatistics statistics;	
	
	int funceval=0;
//	reduced_integer.clear();
	bool box_changed=false;
#ifdef COIN_HAS_FILIB
	IntervalVector box(oldlow, oldup);

	interval<double> oldbounds, newbounds, val;
	while (!nodeset.empty() && funceval<funceval_limit) {
		const DependencyGraph::Node& node(*nodeset.begin()->first);
		BoundType wb=nodeset.begin()->second;
		nodeset.erase(nodeset.begin());
		if (!nodeset.empty() && (**nodeset.begin()->first).varindex==(*node).varindex) { // same node two time in nodeset, so both bounds have changed
			nodeset.erase(nodeset.begin());
			wb=ANY;
		}

		for (DependencyGraph::Node::const_iterator adj(node.begin()); adj!=node.end() && funceval<funceval_limit; ++adj) { // check each adjacent edge
			int index=(**(*adj).head()).varindex;
			for (list<CoinTriple<int, BoundType, double> >::const_iterator it((**adj).coninfo.begin()); it!=(**adj).coninfo.end() && funceval<funceval_limit; ++it) {
				if ((it->second & wb)==0) continue; // required and provided bound change do not match 

				interval<double> oldbounds=box(index);
				if (oldbounds.diam()<tiny_tol) break;
		
				const MINLPData::Constraint& con(data.getConstraint(it->first));
				double coeff=it->third;

				++funceval;
				if (!evalConstraint(val, con, box, index)) continue; // skip this constraint 	
//				val/=-coeff;
				
				// bounds given by interval evaluation
				interval<double> eval_bounds(translateNInftyCoin2Filib(con.lower), translatePInftyCoin2Filib(con.upper));
				eval_bounds-=val;
				eval_bounds/=-coeff;
				
				newbounds=oldbounds.intersect(eval_bounds);
				if (newbounds.isEmpty()) { // maybe rounding error, adding some tolerance
					eval_bounds.blow(tiny_tol);
					newbounds=oldbounds.intersect(eval_bounds);
					if (!newbounds.isEmpty()) newbounds=newbounds.mid(); // make a point interval out of it (since it was empty before) 
				}
				
				if (data.getVariable(index).isDiscrete()) {
					if (newbounds.inf()>oldbounds.inf()+tiny_tol) {
						newbounds=interval<double>(roundUp(newbounds.inf()), newbounds.sup());
						++statistics.shrinked_integer_var;
					} else if (newbounds.sup()<oldbounds.sup()-tiny_tol) {
						newbounds=interval<double>(newbounds.inf(), roundDown(newbounds.sup()));
						++statistics.shrinked_integer_var;
					} else newbounds=oldbounds; // no tiny_tol's in bounds of discrete variables please
				} 

				box[index]=newbounds;
	
				if (oldbounds==newbounds) continue; // no change of box
				box_changed=true;
	
				clog << con.name << ": Reduce box of " << data.getVariable(index).getName() << " from " << oldbounds << " to " << newbounds << "\t (" << val << ')';
				if (newbounds.isEmpty()) {
					statistics.empty_box=true;
					clog << endl;
					break;
				}
	
				bool low_changed, up_changed;
				if (filib::fp_traits<double>::IsInf(oldbounds.inf())) low_changed=!filib::fp_traits<double>::IsInf(newbounds.inf());
				else if (filib::fp_traits<double>::IsInf(oldbounds.sup())) low_changed=newbounds.inf()-oldbounds.inf()>tiny_tol;
				else low_changed=newbounds.inf()-oldbounds.inf()>min_impr*(oldbounds.diam()+1);
				if (filib::fp_traits<double>::IsInf(oldbounds.sup())) up_changed=!filib::fp_traits<double>::IsInf(newbounds.sup());
				else if (filib::fp_traits<double>::IsInf(oldbounds.inf())) up_changed=oldbounds.sup()-newbounds.sup()>tiny_tol;
				else up_changed=oldbounds.sup()-newbounds.sup()>min_impr*(oldbounds.diam()+1);
				clog << "\t low changed: " << low_changed << "\t up changed: " << up_changed << endl;
	
				if (low_changed || up_changed)
					nodeset.insert(
					pair<const DependencyGraph::Node*, BoundType>(&*(*adj).head(),
						low_changed && up_changed ? ANY : (low_changed ? LOWER : UPPER) ));
			}
		}
	}

	clog << "IntervalReduction: " << funceval << " function evaluations";
	if (!box_changed) {
		clog << "\t no reduction" << endl;
		newlow=oldlow;
		newup=oldup;
		return statistics;
	}

	if (statistics.shrinked_integer_var) clog << "\t reduced integers: " << statistics.shrinked_integer_var;
	if (statistics.empty_box) {
		clog << "\t empty box!" << endl;
		statistics.avg_reduction=getInfinity();
		statistics.max_reduction=getInfinity();
		return statistics;
	}

	for (int i=0; i<data.numVariables(); ++i) {
		newlow[i]=translateNInftyFilib2Coin(box(i).inf());
		newup[i]=translatePInftyFilib2Coin(box(i).sup());
		double olddiam=oldup(i)-oldlow(i);
		double red=olddiam ? box(i).diam()/olddiam : 1.;
		if (red>statistics.max_reduction) statistics.max_reduction=red;
		statistics.avg_reduction+=red;
	}
	statistics.avg_reduction/=data.numVariables();

	clog << "\t avg_reduction: " << statistics.avg_reduction;
	clog << endl;
#else
	newlow=oldlow;
	newup=oldup;
#endif // COIN_HAS_FILIB
	return statistics;
}

BoxReductionStatistics ConstraintPropagation::reduceBox() {
	set<pair<const DependencyGraph::Node*, BoundType> > nodeset;
	for (DependencyGraph::iterator it(depgraph.begin()); it!=depgraph.end(); ++it)
		nodeset.insert(pair<const DependencyGraph::Node*, BoundType>(&*it, ANY));
		
	DenseVector oldlow, oldup;
	data.getBox(oldlow, oldup);
	DenseVector newlow(oldlow.getNumElements());
	DenseVector newup(oldlow.getNumElements());

	return run(newlow, newup, oldlow, oldup, nodeset);
}


	
} // namespace LaGO
