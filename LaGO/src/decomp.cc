// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "decomp.h"

#define METIS_AVAILABLE
#ifdef METIS_AVAILABLE
namespace metis {
extern "C" {
#include "metis.h"
}
}
#endif



// ----------------------------------- DecompGraph ---------------------------------------

ostream& operator<<(ostream& out, const DecompGraph& g) {
	out << g.n() << " nodes, " << g.m() << " edges:" << endl;
	for (map<int, DecompGraph::Node>::const_iterator it(g.nodes.begin()); it!=g.nodes.end(); ++it) {
		out << it->first;
		if (it->second.weight) out << " (" << it->second.weight << ')';
		out << ':';
		for (map<int, set<DecompGraph::Edge>::iterator>::const_iterator it2(it->second.adj.begin()); it2!=it->second.adj.end(); ++it2) {
			out << ' ' << it2->first;
			if (it2->second->weight) out << '(' << it2->second->weight << ')';
		}
		out << endl;
	}
	return out;
}

DecompGraph::DecompGraph(const SparsityInfo& si)
: nrcomp(0)
{	for (VariableIterator it(si); it; ++it) add_node(it());
	for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it)
		add_edge(it->first.first, it->first.second);
}

map<int, DecompGraph::Node>::iterator DecompGraph::add_node(int index, int weight) {
	pair<map<int, Node>::iterator, bool> ret(nodes.insert(pair<int, Node>(index, Node(weight))));
	if (!ret.second) ret.first->second.weight+=weight;
	return ret.first;
}

set<DecompGraph::Edge>::iterator DecompGraph::add_edge(map<int, Node>::iterator node1, map<int, Node>::iterator node2, int weight) {
	assert(node1->first!=node2->first);
	if (node1->first>node2->first) return add_edge(node2, node1, weight);

	pair<set<Edge>::iterator, bool> ret(edges.insert(Edge(node1, node2, weight)));
	if (!ret.second) { // edge already there
		ret.first->weight+=weight;
		return ret.first;
	}

	node1->second.adj.insert(pair<int, set<Edge>::iterator>(node2->first, ret.first));
	node2->second.adj.insert(pair<int, set<Edge>::iterator>(node1->first, ret.first));
	return ret.first;
}

set<DecompGraph::Edge>::iterator DecompGraph::add_edge(int node1, int node2, int weight) {
	map<int, Node>::iterator n1(add_node(node1));
	map<int, Node>::iterator n2(add_node(node2));
	return add_edge(n1, n2, weight);
}

int DecompGraph::Node::set_component(int comp) {
	if (component==comp) return 0; // already set
	component=comp;
	int nr=1;
	for (map<int, set<Edge>::iterator>::iterator it(adj.begin()); it!=adj.end(); ++it)
		nr+=it->second->neighbour(*this).set_component(comp);

	return nr;
}

void DecompGraph::compute_connected_components() {
	nrcomp=0;
	map<int, Node>::iterator it_n(nodes.begin());
	while (it_n!=nodes.end()) {
		if (it_n->second.component>=0) { ++it_n; continue; }
		int size=it_n->second.set_component(nrcomp);
//		out_log << nrcomp << ": " << size << endl;

		++it_n;
		++nrcomp;
	}
}

void DecompGraph::get_component_members(vector<list<int> >& members) {
	members.clear();
	members.resize(nrcomp);
	for (map<int, Node>::iterator it(nodes.begin()); it!=nodes.end(); ++it)
		members[it->second.component].push_back(it->first);
}

void DecompGraph::compute_partition(int nparts) {
	nrcomp=0;

	if (nparts==1) {
		nrcomp=1;
		for (map<int, Node>::iterator it(nodes.begin()); it!=nodes.end(); ++it)
			it->second.component=0;
		return;
	}

#ifdef METIS_AVAILABLE
	// compute METIS representation
	int n_=n();
	Pointer<metis::idxtype> xadj(new metis::idxtype[n()+1]);
	Pointer<metis::idxtype> adjncy(new metis::idxtype[2*m()]);

	Pointer<metis::idxtype> vwgt(new metis::idxtype[n()]);
	Pointer<metis::idxtype> adjwgt(new metis::idxtype[2*m()]);

	metis::idxtype complete_edgewgt=0;

	metis::idxtype i=0;
	for (map<int, Node>::iterator it_n(nodes.begin()); it_n!=nodes.end(); ++it_n) {
		vwgt[it_n->first]=it_n->second.weight;
		xadj[it_n->first]=i;
		for (map<int, set<Edge>::iterator>::iterator it_e(it_n->second.adj.begin()); it_e!=it_n->second.adj.end(); ++it_e, ++i) {
			adjncy[i]=it_e->first;
			adjwgt[i]=it_e->second->weight;
			complete_edgewgt+=adjwgt[i];
		}
	}
	xadj[n()]=i;

	int wgtflag=3; // weights on edges and vertices
	int numflag=0; // numbering-style = "C"

	Pointer<int> options(new int[5]); options[0]=0; // default options
	int edgecut;

	Pointer<metis::idxtype> part(new metis::idxtype[n()]); // partition-vector

	metis::METIS_PartGraphKway(&n_, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, part);
//	metis::METIS_PartGraphVKway(&n_, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, part);

	out_log << "Edgecut: " << edgecut << "\t Complete weight of edges: " << complete_edgewgt << endl;

	for (map<int, Node>::iterator it(nodes.begin()); it!=nodes.end(); ++it)
		it->second.component=part[it->first];
	nrcomp=nparts;

#else
	out_err << "Sorry, METIS needed. Aborting." << endl;
	exit(-1);
#endif
}

// ---------------------------------------- SplitFunc --------------------------

void SplitFunc::print(ostream& out) const {
	out << "SplitFunc for indices:";
	for (int i=0; i<indices.size(); i++)
		if (ignore[i]) out << " (" << indices[i] << ")";
		else out << " " << indices[i];
	out << endl << "Original function: " << *orig;
}



// --------------------------------------- SplittingScheme2 --------------------------------------

SplittingScheme2::SplittingScheme2(const DecompGraph& g)
: new_pos(g.n()), newblock(g.nrcomp)
{ for (map<int, DecompGraph::Node>::const_iterator it(g.nodes.begin()); it!=g.nodes.end(); ++it) {
		int comp=it->second.component;
		int new_index=newblock[comp].size();
		newblock[comp].resize(new_index+1);
		newblock[comp][new_index]=it->first;
		new_pos[it->first].push_back(pair<int, int>(comp, new_index));
	}
}

// ------------------------------------- Decomposition ------------------------------------------


void Decomposition::set_sparsity_pattern(MinlpProblem& prob, const vector<vector<dvector> >& sample_set) {
	out_out << "Computing sparsity pattern ";

	double sign_changed=prob.con.size()+1; // percentage number of functions, which are similar to quadratic forms

//	out_log << "Add sparsity pattern for objective:";
	out_log << ".";

	bool sign_changed_=false;
	for (int k=0; k<prob.block.size(); k++)
		sign_changed_=sign_changed_ | prob.obj->compute_sparsity(k, sample_set[k], replace_if_quadratic);
	if (sign_changed_) --sign_changed;

	for (int c=0; c<prob.con.size(); c++) {
		out_log << ".";
		bool sign_changed_=false;
		for (int k=0; k<prob.block.size(); k++) {
			sign_changed_=sign_changed_ | prob.con[c]->compute_sparsity(k, sample_set[k], replace_if_quadratic);
//			out_log << "con. " << c << " block " << k << '\t' << prob.con[c]->get_sparsity(k);
		}
		if (sign_changed_) --sign_changed;
	}

	out_out << endl << "Percentage number of quadratic forms: " << 100*sign_changed/(double)(prob.con.size()+1) << endl;
}


void Decomposition::decompose(SepQcFunc& f, int block_offset, const SepQcFunc& old_f, int k, const SplittingScheme2& ss, Pointer<dvector> primal) {
	const SparsityInfo& si(old_f.get_sparsity(k));

	vector<set<int> > linquad(ss.newblock.size());  // linear and quadratic variables, which should be ignored by SplitFunc, because they were added above
	vector<Pointer<SparsityInfo> > si_new(ss.newblock.size()); // sparsity pattern of new s-functions

	Pointer<dvector> refpoint;
	if (si.nonlinear->size()>si.quadratic->size()) {
		assert(primal);
		refpoint=new dvector(*primal);
		for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
			si_new[new_k]=new SparsityInfo(2);
	}
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
		f.sparsity_block[block_offset+new_k]=new SparsityInfo(2);

	// linear parts
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
		f.b[block_offset+new_k]=new SparseVector<double>(ss.newblock[new_k].size());
	for (map<int, SparsityInfo::LinearVariable>::iterator it(si.linear->begin()); it!=si.linear->end(); ++it) {
		int new_k=ss.new_pos[it->first].front().first;
		int i0=ss.new_pos[it->first].front().second;
		f.b[block_offset+new_k]->SetElement(i0, .5*it->second.coeff);
		linquad[new_k].insert(i0);
		if (refpoint) (*refpoint)[it->first]=0.;
		f.sparsity_block[block_offset+new_k]->add_linear(i0, it->second.coeff);
	}

	// quadratic parts
	if (!si.quadratic->empty()) {
		vector<Pointer<SparseMatrix2> > A_new(ss.newblock.size());
		for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
			A_new[new_k]=new SparseMatrix2(ss.newblock[new_k].size());

		for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it) {
			int i=it->first.first;
			int j=it->first.second;
			int new_k=ss.new_pos[i].front().first;
			assert(new_k==ss.new_pos[j].front().first); // both in same block
			int i0=ss.new_pos[i].front().second;
			int j0=ss.new_pos[j].front().second;
			f.sparsity_block[block_offset+new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
			if ((!si.quadratic->count(i)) || (!si.quadratic->count(j))) {
//				if (si.quadratic->count(i) || si.quadratic->count(j)) out_log << si;
				assert((!si.quadratic->count(i)) && (!si.quadratic->count(j)));
				si_new[new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
				continue;
			}
			A_new[new_k]->AddElement(i0, j0, .5*it->second.coeff);
			A_new[new_k]->AddElement(j0, i0, .5*it->second.coeff);
		}
		for (map<int, SparsityInfo::QuadraticVariable>::iterator it(si.quadratic->begin()); it!=si.quadratic->end(); ++it) {
			int new_k=ss.new_pos[it->first].front().first;
			int i0=ss.new_pos[it->first].front().second;
			A_new[new_k]->AddElement(i0, i0, it->second.coeff_quad);
			if (it->second.coeff_lin) (*f.b[block_offset+new_k])[i0]+=.5*it->second.coeff_lin;
			linquad[new_k].insert(i0);
			if (refpoint) (*refpoint)[it->first]=0.;
			f.sparsity_block[block_offset+new_k]->add_quadratic(i0, it->second.coeff_lin, it->second.coeff_quad);
		}
		for (int new_k=0; new_k<ss.newblock.size(); ++new_k) {
			if (!A_new[new_k]->nonzeros()) continue;
			A_new[new_k]->finish();
			f.A[block_offset+new_k]=A_new[new_k];
		}
	} else { // compute sparsity pattern
		for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it) {
			int i=it->first.first;
			int j=it->first.second;
			int new_k=ss.new_pos[i].front().first;
			assert(new_k==ss.new_pos[j].front().first); // both in same block
			int i0=ss.new_pos[i].front().second;
			int j0=ss.new_pos[j].front().second;
			si_new[new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
			f.sparsity_block[block_offset+new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
		}
	}

	// clear empty linear parts
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
		if (*f.b[block_offset+new_k]==0.) f.b[block_offset+new_k]=NULL;

	// nonlinear, nonquadratic parts
	if (si.nonlinear->size()>si.quadratic->size()) {
		// check, for which blocks exist nonlinear/nonquadratic variables
		vector<bool> nonquadlin(ss.newblock.size(), false);
		for (map<int, SparsityInfo::NonlinearVariable>::iterator it(si.nonlinear->begin()); it!=si.nonlinear->end(); ++it) {
			int new_k=ss.new_pos[it->first].front().first;
			int i0=ss.new_pos[it->first].front().second;
			f.sparsity_block[block_offset+new_k]->add_nonlinear(i0);
			if (si.quadratic->count(it->first)) continue;
			nonquadlin[new_k]=true;
//			out_log << "nonlinquad: " << k << "\t " << new_k << endl;
			si_new[new_k]->add_nonlinear(i0);
		}

		int nr=0;
		double val;
		for (int new_k=0; new_k<ss.newblock.size(); ++new_k) {
			if (nonquadlin[new_k]) {
				f.s[block_offset+new_k]=new SplitFunc(old_f.s[k], ss.newblock[new_k], *refpoint, linquad[new_k], si_new[new_k]);
				++nr;
			}
		}
		if (nr>1) f.c-=(nr-1)*old_f.s[k]->eval(*refpoint);
	} else {
		if (old_f.s[k]) f.c+=old_f.s[k]->eval(SparseVector<double>(old_f.block[k].dim())); // value in 0
		for (map<int, SparsityInfo::NonlinearVariable>::iterator it(si.nonlinear->begin()); it!=si.nonlinear->end(); ++it)
			f.sparsity_block[block_offset+ss.new_pos[it->first].front().first]->add_nonlinear(ss.new_pos[it->first].front().second);
	}

	// curvature
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
		if ((!f.A[block_offset+new_k]) && (!f.s[block_offset+new_k]))
			f.set_curvature(block_offset+new_k, Func::LINEAR);
		else {
			Func::CurvatureType oldcurv(old_f.get_curvature(k));
			f.set_curvature(block_offset+new_k, (oldcurv&Func::LINEAR) ? oldcurv : Func::UNKNOWN); // convex or concave or both
		}
		
	// ensure, that variables in linear blocks are known to be linear
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k) {
		if (f.b[block_offset+new_k] && (!f.A[block_offset+new_k]) && (!f.s[block_offset+new_k])) {
			assert(f.sparsity_block[block_offset+new_k]->sparsity_pattern->empty());
			for (VariableIterator it(*f.sparsity_block[block_offset+new_k], false, false, true); it; ++it) {
//				out_log << "changed type of variable " << it() << " in new block " << block_offset+new_k << endl;
				f.sparsity_block[block_offset+new_k]->add_linear(it(), it.coeff_lin());
			}
			f.sparsity_block[block_offset+new_k]->quadratic->clear();
			f.sparsity_block[block_offset+new_k]->nonlinear->clear();
		}
	}
		
}



Pointer<MinlpProblem> Decomposition::decompose(MinlpProblem& prob, vector<Pointer<SplittingScheme2> >& ss) {
	Pointer<MinlpProblem> dprob(new MinlpProblem(prob.problem_type, prob.prob_name, prob.out_problem_p, prob.out_problem_log_p));
	bool nonquadlin=false; // decomposed problem is not quadratic/linear ?

	vector<Pointer<dvector> > primal(prob.block.size());
	if (prob.problem_type==MinlpProblem::MINLP)
		for (int k=0; k<prob.block.size(); ++k) primal[k]=new dvector(prob.primal_point, prob.block[k]);

	int new_blocknr=0;
	for (int k=0; k<ss.size(); ++k) new_blocknr+=ss[k]->newblock.size();

	int acc_blocknr=0; // accumalated (new) block nr
	// create new variable structure and set discrete variables, lower and upper bounds
	for (int k=0; k<ss.size(); ++k) // old blocks
		for (int new_k=0; new_k<ss[k]->newblock.size(); ++new_k, ++acc_blocknr) // new blocks inside old blocks
			for (int i=0; i<ss[k]->newblock[new_k].size(); ++i) { // variables inside new blocks
				int oldindex=prob.block[k][ss[k]->newblock[new_k][i]];
				dprob->add_var(ss[k]->newblock[new_k][i], acc_blocknr, prob.discr[oldindex], prob.lower[oldindex], prob.upper[oldindex], prob.var_names[oldindex]);
			}

	out_out << "Decompose objective and constraints "; out_log << ".";

	// decompose objective
	Pointer<SepQcFunc> obj(new SepQcFunc(dprob->block));
	obj->c=prob.obj->c;
	acc_blocknr=0;
	for (int k=0; k<ss.size(); ++k) {
		decompose(*obj, acc_blocknr, *prob.obj, k, *ss[k], primal[k]);
		acc_blocknr+=ss[k]->newblock.size();
	}
	dprob->add_obj(obj);

	for (int k=0; (!nonquadlin) && k<obj->block.size(); k++) nonquadlin|=(bool)dprob->obj->s[k];

	// add decomposed constraints
	for (int c=0; c<prob.con.size(); c++) {
		out_log << ".";
//		out_log << "Decompose constraint " << prob.con_names[c] << endl;
		Pointer<SepQcFunc> con(new SepQcFunc(dprob->block));
		con->c=prob.con[c]->c;
		acc_blocknr=0;
		for (int k=0; k<ss.size(); ++k) {
			decompose(*con, acc_blocknr, *prob.con[c], k, *ss[k], primal[k]);
			acc_blocknr+=ss[k]->newblock.size();
		}
		dprob->add_con(con, prob.con_eq[c], prob.con_names[c]);
		for (int k=0; (!nonquadlin) && k<con->block.size(); k++) nonquadlin|=(bool)con->s[k];
	}
	out_out << "done" << endl;

	// set problem type, if its now quadratic
	if (!nonquadlin) dprob->problem_type=MinlpProblem::QQP;

	dprob->primal_point=prob.primal_point;
//out_log << *dprob;
	return dprob;
}



Pointer<MinlpProblem> Decomposition::decompose(MinlpProblem& prob, vector<vector<dvector> >& sample_set) {
	// compute sparsity-pattern
	set_sparsity_pattern(prob, sample_set);

	// decompose sparsity-graphs and compute splitting schemes
	vector<Pointer<SplittingScheme2> > ss2(prob.block.size());
	for (int k=0; k<prob.block.size(); k++) ss2[k]=get_splittingscheme(prob, k);

	return decompose(prob, ss2);
}

// generate "local" splitting scheme
Pointer<SplittingScheme2> Decomposition::get_splittingscheme(MinlpProblem& prob, int blocknr) {
	Pointer<SparsityInfo> si(prob.get_sparsity(blocknr));

	DecompGraph g(*si);
	for (int i=0; i<prob.block[blocknr].size(); ++i) g.add_node(i);
	g.compute_connected_components();
	out_out << "Number of small blocks: " << g.nrcomp << endl;

	// merge small blocks again, which are close together
	vector<list<int> > components;
	g.get_component_members(components);

//	out_log << "Before: " << endl;
	DecompGraph h; int i=0;
	double avgsize=0;
	for (vector<list<int> >::iterator it(components.begin()); it!=components.end(); ++it, ++i) {
		int size=it->size();
		avgsize+=size;
		h.add_node(i, size);
//		out_log << i << ':'; for (list<int>::iterator it2(it->begin()); it2!=it->end(); ++it2) out_log << ' ' << *it2; out_log << endl;
	}
	avgsize/=h.n();
	out_log << "Average size: " << avgsize << endl;
	for (int c=0; c<=prob.con.size(); ++c) {
		SparsityInfo& si((c ? prob.con[c-1] : prob.obj)->get_sparsity(blocknr));
		set<pair<int,int> > connects;
		for (VariableIterator it1(si); it1; ++it1)
			for (VariableIterator it2(si); it2; ++it2) {
				if (it1()>=it2()) continue;
				int comp1=g.nodes[it1()].component;
				int comp2=g.nodes[it2()].component;
				if (comp1!=comp2) connects.insert(pair<int,int>(comp1, comp2));
			}
		for (set<pair<int,int> >::iterator it(connects.begin()); it!=connects.end(); ++it)
			h.add_edge(it->first, it->second, (c && prob.con_eq[c-1]) ? 2 : 1);
	}
//	out_log << h;

	int nparts=(int)MAX(1, prob.block[blocknr].size()/(avg_blocksize*avgsize));
	h.compute_partition(nparts);

	vector<list<int> > comp2;
	h.get_component_members(comp2);

	// clear empty components
	while (comp2.back().empty()) {
		comp2.pop_back();
		--h.nrcomp;
	}
	i=0;
	for (vector<list<int> >::iterator it(comp2.begin()); it!=comp2.end(); ++it, ++i)
		if (it->empty()) { // empty component
			*it=comp2.back();
			for (list<int>::iterator it2(it->begin()); it2!=it->end(); ++it2)
				h.nodes[*it2].component=i;
			comp2.pop_back();
			while (comp2.back().empty()) {
				comp2.pop_back();
				--h.nrcomp;
			}
			--h.nrcomp;
		}

	// set new components
	for (map<int, DecompGraph::Node>::iterator it(g.nodes.begin()); it!=g.nodes.end(); ++it)
		it->second.component=h.nodes[it->second.component].component;
	g.nrcomp=h.nrcomp;

	out_out << "Number of bigger blocks: " << g.nrcomp << endl;
//	out_log << "Now: " << endl;
	g.get_component_members(components);
	i=0;
	avgsize=0;
	int maxsize=0, minsize=INF;
	for (vector<list<int> >::iterator it(components.begin()); it!=components.end(); ++it, ++i) {
		int size=it->size();
		avgsize+=size;
		if (size>maxsize) maxsize=size;
		if (size<minsize) minsize=size;
//		out_log << i << ':'; for (list<int>::iterator it2(it->begin()); it2!=it->end(); ++it2) out_log << '\t' << *it2; out_log << endl;
	}
	avgsize/=g.nrcomp;
	out_log << "Avg. block size: " << avgsize << endl;
	out_log << "Max. block size: " << maxsize << endl;
	out_log << "Min. block size: " << minsize << endl;

	return new SplittingScheme2(g);
}

Pointer<SepQcFunc> Decomposition::decompose(Pointer<Func> old_f, const dvector& primal) {
	const SparsityInfo& si(((const Func*)(Func*)old_f)->get_sparsity());
	DecompGraph g(si);
	g.compute_connected_components();
	bool addedone=false;
	for (int i=0; i<old_f->dim(); ++i) { // add nonappearing variables in extra component
		if (si.linear->count(i) || si.nonlinear->count(i)) continue;
		g.add_node(i)->second.component=g.nrcomp;
		addedone=true;
	}
	if (addedone) ++g.nrcomp;
	SplittingScheme2 ss(g);

	Pointer<SepQcFunc> ret(new SepQcFunc(ss.newblock));
	SepQcFunc& f(*ret);

	vector<set<int> > linquad(ss.newblock.size());  // linear and quadratic variables, which should be ignored by SplitFunc, because they were added above
	vector<Pointer<SparsityInfo> > si_new(ss.newblock.size()); // sparsity pattern of new s-functions

	dvector refpoint(primal);
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k) {
		si_new[new_k]=new SparsityInfo(2);
		f.sparsity_block[new_k]=new SparsityInfo(2);
	}

	// linear parts
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
		f.b[new_k]=new SparseVector<double>(ss.newblock[new_k].size());
	for (map<int, SparsityInfo::LinearVariable>::iterator it(si.linear->begin()); it!=si.linear->end(); ++it) {
		int new_k=ss.new_pos[it->first].front().first;
		int i0=ss.new_pos[it->first].front().second;
		f.b[new_k]->SetElement(i0, .5*it->second.coeff);
		linquad[new_k].insert(i0);
		refpoint[it->first]=0.;
		f.sparsity_block[new_k]->add_linear(i0, it->second.coeff);
	}

	// quadratic parts
	if (!si.quadratic->empty()) {
		vector<Pointer<SparseMatrix2> > A_new(ss.newblock.size());
		for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
			A_new[new_k]=new SparseMatrix2(ss.newblock[new_k].size());

		for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it) {
			int i=it->first.first;
			int j=it->first.second;
			int new_k=ss.new_pos[i].front().first;
			assert(new_k==ss.new_pos[j].front().first); // both in same block
			int i0=ss.new_pos[i].front().second;
			int j0=ss.new_pos[j].front().second;
			f.sparsity_block[new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
			if ((!si.quadratic->count(i)) || (!si.quadratic->count(j))) {
				si_new[new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
				continue;
			}
			A_new[new_k]->AddElement(i0, j0, .5*it->second.coeff);
			A_new[new_k]->AddElement(j0, i0, .5*it->second.coeff);
		}
		for (map<int, SparsityInfo::QuadraticVariable>::iterator it(si.quadratic->begin()); it!=si.quadratic->end(); ++it) {
			int new_k=ss.new_pos[it->first].front().first;
			int i0=ss.new_pos[it->first].front().second;
			A_new[new_k]->AddElement(i0, i0, it->second.coeff_quad);
			if (it->second.coeff_lin) (*f.b[new_k])[i0]+=.5*it->second.coeff_lin;
			linquad[new_k].insert(i0);
			refpoint[it->first]=0.;
			f.sparsity_block[new_k]->add_quadratic(i0, it->second.coeff_lin, it->second.coeff_quad);
		}
		for (int new_k=0; new_k<ss.newblock.size(); ++new_k) {
			if (!A_new[new_k]->nonzeros()) continue;
			A_new[new_k]->finish();
			f.A[new_k]=A_new[new_k];
		}
	} else { // compute sparsity pattern
		for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it) {
			int i=it->first.first;
			int j=it->first.second;
			int new_k=ss.new_pos[i].front().first;
			assert(new_k==ss.new_pos[j].front().first); // both in same block
			int i0=ss.new_pos[i].front().second;
			int j0=ss.new_pos[j].front().second;
			si_new[new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
			f.sparsity_block[new_k]->add_sparsity_pattern(i0, j0, it->second.coeff);
		}
	}

	// clear empty linear parts
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
		if (*f.b[new_k]==0.) f.b[new_k]=NULL;

	// nonlinear, nonquadratic parts
	if (si.nonlinear->size()>si.quadratic->size()) {
		// check, for which blocks exist nonlinear/nonquadratic variables
		vector<bool> nonquadlin(ss.newblock.size(), false);
		for (map<int, SparsityInfo::NonlinearVariable>::iterator it(si.nonlinear->begin()); it!=si.nonlinear->end(); ++it) {
			int new_k=ss.new_pos[it->first].front().first;
			int i0=ss.new_pos[it->first].front().second;
			f.sparsity_block[new_k]->add_nonlinear(i0);
			if (si.quadratic->count(it->first)) continue;
			nonquadlin[new_k]=true;
			si_new[new_k]->add_nonlinear(i0);
		}

		int nr=0;
		double val;
		for (int new_k=0; new_k<ss.newblock.size(); ++new_k) {
			if (nonquadlin[new_k]) {
				f.s[new_k]=new SplitFunc(old_f, ss.newblock[new_k], refpoint, linquad[new_k], si_new[new_k]);
				++nr;
			}
		}
		if (nr>1) f.c-=(nr-1)*old_f->eval(refpoint);
	} else {
		f.c+=old_f->eval(SparseVector<double>(old_f->dim())); // value in 0
		for (map<int, SparsityInfo::NonlinearVariable>::iterator it(si.nonlinear->begin()); it!=si.nonlinear->end(); ++it)
			f.sparsity_block[ss.new_pos[it->first].front().first]->add_nonlinear(ss.new_pos[it->first].front().second);
	}

	// curvature
	for (int new_k=0; new_k<ss.newblock.size(); ++new_k)
		if ((!f.A[new_k]) && (!f.s[new_k]))
			f.set_curvature(new_k, Func::LINEAR);
		else {
			Func::CurvatureType oldcurv(old_f->get_curvature());
			f.set_curvature(new_k, (oldcurv&Func::LINEAR) ? oldcurv : Func::INDEFINITE); // convex or concave or both
		}
		
	return ret;
}
