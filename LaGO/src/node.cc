// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "node.h"

MinlpNode::MinlpNode(const dvector& lower_, const dvector& upper_/*, const map<int, set<int> >& already_fixed*/)
: low_bound(-INFINITY), ref_point(lower_.dim()),
  lower(lower_), upper(upper_), update_subdiv_bound_called(false)
{ }

MinlpNode::MinlpNode(MinlpNode& node)
: low_bound(node.low_bound), ref_point(node.ref_point), dual_point(node.dual_point),
	/*bcp_fixed_var(node.bcp_fixed_var), part_con(node.part_con), */lower(node.lower), upper(node.upper),
	fix_branch_var(node.fix_branch_var), update_subdiv_bound_called(false), lagprob_solutions(node.lagprob_solutions)
{ }

pair<double, pair<int,int> > MinlpNode::bcp_rho(const dvector& x, const vector<ivector>& block, const vector<bool>& discr) {
	pair<double, pair<int,int> > ret(INFINITY, pair<int,int>(-1,-1));
	for (int k=0; k<block.size(); k++) {
		for (int i=0; i<block[k].size(); i++) {
			int i0=block[k][i];
			if (i0>=discr.size() || (!discr[i0])) continue; // not discrete variable
			if (upper(i0)-lower(i0)<rtol) continue; // already fixed
			double rho=fabs(x(i0)-.5*(upper(i0)+lower(i0))) / (upper(i0)-lower(i0));
			if (rho<ret.first) {
				ret.first=rho;
				ret.second=pair<int,int>(k,i);
			}
		}
	}

	return ret;
}

set<SolCandidate>::const_iterator MinlpNode::outside_part_set(const set<SolCandidate>& points) {
	if (points.empty()) return points.end();
	set<SolCandidate>::const_iterator it(points.begin());
	bool feas;
	do {
		feas=true;
		for (int i=0; feas && i<it->second.dim(); ++i)
			feas=(lower(i)<=it->second(i)+rtol) && (upper(i)>=it->second(i)-rtol);
//		for (int k=0; feas && k<part_con.size(); k++)
//			for (list<Pointer<SepQcFunc> >::iterator it_con(part_con[k].begin()); feas && it_con!=part_con[k].end(); it_con++)
//				feas=(*it_con)->eval(it->second)<1E-4;
	} while ((!feas) && (++it)!=points.end());
	return it;
}

bool MinlpNode::inside_part_set(dvector &point, int k, const vector<ivector>& block) {
	bool feas=true;

	for (int i=0; feas && i<point.size(); ++i)
		feas=(lower(block[k][i])<=point(i)+rtol) && (upper(block[k][i])>=point(i)-rtol);
//	if (k<part_con.size())
//		for (list<Pointer<SepQcFunc> >::iterator it(part_con[k].begin()); feas && it!=part_con[k].end(); it++)
//			feas=(*it)->eval(point, k)+(*it)->c<1E-4;

	return feas;
}

double MinlpNode::key(const vector<int>& i_discr) {
	bool allfixed=true;
	for (int i=0; i<i_discr.size(); ++i)
		if (lower(i_discr[i])!=upper(i_discr[i])) {
			allfixed=false;
			break;
		} 
	return key() + (allfixed ? rtol : 0.);
}
