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

void ConstraintPropagation::initDependencyGraph() {
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
			if (i_linvar<con.decompfuncLin->getNumElements()-1 && linvar[i_linvar]<i) {
				++i_linvar;
				assert(linvar[i_linvar-1]<linvar[i_linvar]);
			}
			if ((i_linvar>=con.decompfuncLin->getNumElements() || linvar[i_linvar]!=i) && (con.decompfuncNL.empty() || con.decompmapping[i].empty())) continue; // variable not in function
			
			for (int i_linvar2=0; i_linvar2<con.decompfuncLin->getNumElements(); ++i_linvar2) {
				if (linvar[i_linvar2]==i) continue;
				
				BoundType wb;
				// if var i is nonlinear, I cannot say what effect a change on one of its bounds has on var i_linvar2 can have
				if ((!con.decompfuncNL.empty()) && !con.decompmapping[i].empty()) wb=ANY;
				else if (con.lower>-getInfinity() && con.upper<getInfinity()) wb=ANY;
				else if (con.lower>-getInfinity()) wb=(linvar_el[i_linvar]<0 ? LOWER : UPPER);
				else if (con.upper< getInfinity()) wb=(linvar_el[i_linvar]<0 ? UPPER : LOWER);
				else continue; // free constraint

				DependencyGraph::arc_iterator arc=depgraph.arc_find(nodes[i], nodes[linvar[i_linvar2]]);
				if (arc==depgraph.arc_end()) // new arc
					depgraph.arc_insert(nodes[i], EdgeInfo(c, wb, linvar_el[i_linvar2]), nodes[linvar[i_linvar2]]);
				else
					(**arc).add(c, wb, linvar_el[i_linvar2]);
			}			
		}
	}
	if (print_level>4)
		clog << depgraph;

	if (print_level)
		clog << "ConstraintPropagation dependency graph nr. edges: " << depgraph.arc_size() << " (" << dummy_node->size() << ')' << endl;
}

//bool ConstraintPropagation::evalConstraint(interval<double>& val, const MINLPData::Constraint& con, IntervalVector& box, int index) {
//	val=con.decompfuncConstant;
//	try {
//		for (unsigned int k=0; k<con.decompfuncNL.size(); ++k) {
//			if (!con.decompfuncNL[k]->canIntervalEvaluation()) return false;
//			IntervalVector block(box, con.decompfuncNL[k]->indices);
//			val+=con.decompfuncNL[k]->eval(block);					
//		}
//	} catch (FunctionEvaluationError error) {
//		cerr << "Error evaluating block of constraint " << con.name << " over an interval: " << error;
//		return false;					
//	} 
//	box[index]=interval<double>::ZERO();
//	val+=box* *con.decompfuncLin;
//	return true;
//}

bool ConstraintPropagation::evalConstraint(interval<double>& val, const MINLPData::Constraint& con, IntervalVector& box, int index) {
	val=con.decompfuncConstant;
	if (IsNull(con.decompfuncLin) || con.decompfuncLin->getNumElements()==0) { // easy way
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
//		box[index]=interval<double>::ZERO();
		return true;
	}

	DenseVector lin(data.numVariables(), *con.decompfuncLin);
	lin[index]=0.;	
	try {
		for (unsigned int k=0; k<con.decompfuncNL.size(); ++k) {
			const BlockFunction& block_func(*con.decompfuncNL[k]);
			if (!block_func.canIntervalEvaluation()) return false;
			IntervalVector block(box, block_func.indices);
			if (IsValid(block_func.quad)) { // we let evaluate xAx+bx here to obtain tighter interval
				DenseVector lin_block(lin, block_func.indices);
				val+=block_func.quad->xAx_bx(block, lin_block);
//					
//				clog << "quad. over block: val= " << val << endl;
//				clog << "x: " << block << endl;
//				clog << *block_func.quad << endl;
//				clog << "lin: " << lin_block << endl;

				for (int i=0; i<lin_block.getNumElements(); ++i)
					lin[block_func.indices[i]]=0.;
			}					
			if (IsValid(block_func.nonquad))
				val+=block_func.nonquad->eval(block);
		}
	} catch (FunctionEvaluationError error) {
		cerr << "Error evaluating block of constraint " << con.name << " over an interval: " << error;
		return false;					
	}
	val+=box*lin;
//
//	interval<double> val_save=val;	
//	val=con.decompfuncConstant;
//	try {
//		for (unsigned int k=0; k<con.decompfuncNL.size(); ++k) {
//			if (!con.decompfuncNL[k]->canIntervalEvaluation()) return false;
//			IntervalVector block(box, con.decompfuncNL[k]->indices);
//			val+=con.decompfuncNL[k]->eval(block);					
//		}
//	} catch (FunctionEvaluationError error) {
//		cerr << "Error evaluating block of constraint " << con.name << " over an interval: " << error;
//		return false;					
//	} 
//	box[index]=interval<double>::ZERO();
//	val+=box* *con.decompfuncLin;
//	
//	clog << con.name << ' ' << index << ": " << val_save << '\t' << val << endl; 
////	clog << box << endl;
////	clog << lin;
//	assert(val==val_save);
	
	return true;
}

BoxReductionStatistics ConstraintPropagation::run(DenseVector& newlow, DenseVector& newup, const DenseVector& oldlow, const DenseVector& oldup, set<pair<const DependencyGraph::Node*, BoundType> >& nodeset) {
	BoxReductionStatistics statistics;	
	
	int funceval=0;
	bool box_changed=false;
#ifdef COIN_HAS_FILIB
	IntervalVector box(oldlow, oldup);

	interval<double> oldbounds, newbounds, val;
	while (!nodeset.empty() && funceval<funceval_limit && !statistics.empty_box) {
		const DependencyGraph::Node& node(*nodeset.begin()->first);
		BoundType wb=nodeset.begin()->second;
		nodeset.erase(nodeset.begin());
		if (!nodeset.empty() && (**nodeset.begin()->first).varindex==(*node).varindex) { // same node two time in nodeset, so both bounds have changed
			nodeset.erase(nodeset.begin());
			wb=ANY;
		}

		for (DependencyGraph::Node::const_iterator adj(node.begin()); adj!=node.end() && funceval<funceval_limit; ++adj) { // check each adjacent edge
			int index=(**(*adj).head()).varindex;
			bool low_changed=false, up_changed=false;
			for (list<CoinTriple<int, BoundType, double> >::const_iterator it((**adj).coninfo.begin()); it!=(**adj).coninfo.end() && funceval<funceval_limit; ++it) {
				if ((it->second & wb)==0)
					continue; // required and provided bound change do not match

				oldbounds=box(index);
				if (oldbounds.diam()<getTinyTol()) break;

				const MINLPData::Constraint& con(data.getConstraint(it->first));
				double coeff=it->third;


				// assume the constraint has the form a*x + h(x,y) \in [l,u], where x is given by index and a by coeff
				// a box for h(x,y) is evaluated by evalConstraint(val, con, box, index)
				// then a new box for x is given by eval_bounds := ([l, u] - h(x,y))/a.

				++funceval;
				if (!evalConstraint(val, con, box, index)) continue; // skip this constraint 	
				
				// bounds given by interval evaluation
				interval<double> eval_bounds(translateNInftyCoin2Filib(con.lower), translatePInftyCoin2Filib(con.upper));
				eval_bounds-=val;
				eval_bounds/=coeff;

				if (print_level>2)
					clog << "Con. " << con.name << " Var. " << data.getVariable(index).getName() << ": value: " << val << " oldbox: " << oldbounds << " prop. box: " << eval_bounds << endl;
				
				newbounds=oldbounds.intersect(eval_bounds);
				if (newbounds.isEmpty()) { // maybe rounding error, adding some tolerance
					eval_bounds.blow(getTinyTol());
					newbounds=oldbounds.intersect(eval_bounds);
					if (!newbounds.isEmpty()) newbounds=newbounds.mid(); // make a point interval out of it (since it was empty before) 
				}
				
				if (data.getVariable(index).isDiscrete()) {
					if (newbounds.inf()>oldbounds.inf()+getTinyTol()) {
						newbounds=interval<double>(roundUp(newbounds.inf()), newbounds.sup());
						++statistics.shrinked_integer_var;
					} else if (newbounds.sup()<oldbounds.sup()-getTinyTol()) {
						newbounds=interval<double>(newbounds.inf(), roundDown(newbounds.sup()));
						++statistics.shrinked_integer_var;
					} else newbounds=oldbounds; // no tiny tol's in bounds of discrete variables please
				} 

				box[index]=newbounds;
	
				if (oldbounds==newbounds) continue; // no change of box
				box_changed=true;
	
				if (print_level>1)
					clog << con.name << ": Reduce box of " << data.getVariable(index).getName() << " from " << oldbounds << " to " << newbounds << "\t (val=" << val << " eval_bounds=" << eval_bounds << ')' << endl;
				if (newbounds.isEmpty()) {
					statistics.empty_box=true;
					low_changed=up_changed=true;
					break;
				}
	
				if (!low_changed)
					if (filib::fp_traits<double>::IsInf(oldbounds.inf())) low_changed=!filib::fp_traits<double>::IsInf(newbounds.inf());
					else if (filib::fp_traits<double>::IsInf(oldbounds.sup())) low_changed=newbounds.inf()-oldbounds.inf()>getTinyTol();
					else low_changed=newbounds.inf()-oldbounds.inf()>min_impr*(oldbounds.diam()+1);
				if (!up_changed)
					if (filib::fp_traits<double>::IsInf(oldbounds.sup())) up_changed=!filib::fp_traits<double>::IsInf(newbounds.sup());
					else if (filib::fp_traits<double>::IsInf(oldbounds.inf())) up_changed=oldbounds.sup()-newbounds.sup()>getTinyTol();
					else up_changed=oldbounds.sup()-newbounds.sup()>min_impr*(oldbounds.diam()+1);
			}

			if (statistics.empty_box) break;
			if (low_changed || up_changed) {
//				clog << data.getVariable(index).getName() << " low changed: " << low_changed << "\t up changed: " << up_changed << endl;
				nodeset.insert(
				pair<const DependencyGraph::Node*, BoundType>(&*(*adj).head(),
					low_changed && up_changed ? ANY : (low_changed ? LOWER : UPPER) ));
			}
		}
	}

	if (print_level)
		clog << "IntervalReduction: " << funceval << " function evaluations";
	if (!box_changed) {
		if (print_level)
			clog << "\t no reduction" << endl;
		newlow=oldlow;
		newup=oldup;
		return statistics;
	}

	if (print_level && statistics.shrinked_integer_var) clog << "\t reduced integers: " << statistics.shrinked_integer_var;
	if (statistics.empty_box) {
		if (print_level)
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

	if (print_level)
		clog << "\t avg_reduction: " << statistics.avg_reduction << "\t max_reduction: " << statistics.max_reduction << endl;
#else
	newlow=oldlow;
	newup=oldup;
#endif // COIN_HAS_FILIB
	return statistics;
}

BoxReductionStatistics ConstraintPropagation::run(DenseVector& newlow, DenseVector& newup, const DenseVector& oldlow, const DenseVector& oldup) {
	set<pair<const DependencyGraph::Node*, BoundType> > nodeset;
	for (DependencyGraph::iterator it(depgraph.begin()); it!=depgraph.end(); ++it)
		nodeset.insert(pair<const DependencyGraph::Node*, BoundType>(&*it, ANY));

	return run(newlow, newup, oldlow, oldup, nodeset);
}


BoxReductionStatistics ConstraintPropagation::reduceBox() {
	set<pair<const DependencyGraph::Node*, BoundType> > nodeset;
	for (DependencyGraph::iterator it(depgraph.begin()); it!=depgraph.end(); ++it)
		nodeset.insert(pair<const DependencyGraph::Node*, BoundType>(&*it, ANY));
		
	DenseVector oldlow, oldup;
	data.getBox(oldlow, oldup);
	DenseVector newlow(oldlow.getNumElements());
	DenseVector newup(oldlow.getNumElements());

	BoxReductionStatistics statistics(run(newlow, newup, oldlow, oldup, nodeset));
	
	if (statistics.avg_reduction>0) { // some reduction achieved
		for (int i=0; i<data.numVariables(); ++i) {
			data.var[i].lower=newlow[i];
			data.var[i].upper=newup[i];
		}		
	}
	
	return statistics;
}


	
} // namespace LaGO
