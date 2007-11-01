// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "cuts.h"
#include "node.h"

#include <map>

// ---------------------- IntervalGradientCut ------------------------------------

Pointer<IntervalGradientCut> IntervalGradientCutGenerator::get_cuts(const dvector& x, int k, const dvector& low, const dvector& up) {
	Pointer<IntervalGradientCut> cut;
#ifdef FILIB_AVAILABLE
	IntervalVector X(low, up); // X = [low x, up x]

	dvector x_ref(x); // reference point; might be moved a bit, if close to some bounds
	ivector omit(x.dim(), 0); // the indices to omit, if >0, and the reason, why to do so: 1: on lower bound, 2: on upper bound, 3: linear variable
	int size=0;
	for (int i=0; i<x.dim(); ++i) {
		if (low(i)<=-INFINITY) X[i]=interval<double>(filib::fp_traits<double>::ninfinity(), up(i));
		if (up(i)>=INFINITY) X[i]=interval<double>(X(i).inf(), filib::fp_traits<double>::infinity());

		if (sparsity[k] && !sparsity[k]->nonlinear->count(i)) omit[i]=3;
		else if (x(i)-low(i)<=1E-4*(up(i)-low(i))) { omit[i]=1; x_ref[i]=low(i); }
		else if (up(i)-x(i)<=1E-4*(up(i)-low(i))) { omit[i]=2; x_ref[i]=up(i); }
		else ++size; // reference point not on bounds
	}

	if (!size) return NULL;

	dvector w_z_up(2*size);
	ivector indices(size);
	for (int i=0, j=0; i<x.dim(); ++i)
		if (!omit[i]) {
			w_z_up[j]=x(i)-low(i); // w(i)
			w_z_up[size+j]=up(i)-x(i); // z(i)
			indices[j]=i;
			++j;
		}
//	out_log << "x: " << x << "X: " << X;
//	out_log << "indices: " << indices << "omit: " << omit;

	cut=new IntervalGradientCut(x_ref, w_z_up, indices);

	IntervalVector intgrad(X.dim()); // \nabla f(X)
	for (int c=0; c<prob->con.size(); c++) {
		if ((!prob->con[c]->A[k]) && (!prob->con[c]->s[k])) continue;
		if (!prob->con[c]->is_interval_compliant()) {
			out_log << prob->con_names[c] << " is not interval compliant. Skipping IntervalGradientCut." << endl;
			continue;
		}
		Func::CurvatureType ct=prob->con[c]->get_curvature(k);
		if (ct&Func::CONVEX && ((!prob->con_eq[c]) || ct&Func::CONCAVE)) continue;

//		out_log << "Computing IntervalGradientCut for con " << prob->con_names[c] << endl;
//		prob->con[c]->print(*out_log_p, prob->var_names);

		double val=prob->con[c]->eval(x_ref, k)+prob->con[c]->c;
		if (!finite(val)) continue;

		prob->con[c]->grad(intgrad, X, k);
		bool finite=true;
		for (int i=0; i<X.dim(); ++i)
			if (intgrad[i].isInfinite()) { finite=false; break; }
		if (!finite) continue;

		dvector grad(x_ref.dim());
		prob->con[c]->grad(grad, x_ref, k);

//		for (int i=0; i<X.dim(); ++i)
//			if (omit[i]==3 && !intgrad[i].isPoint()) { out_err << intgrad[i] << '\t' << grad(i) << endl; }

//		out_log << "grad: " << grad << "intgrad: " << intgrad << "value: " << val << endl;

		Pointer<UserVector<double> > b_w(new dvector(size));
		Pointer<UserVector<double> > b_z(new dvector(size));

		int i0;
		for (int i=0; i<indices.dim(); ++i) {
			i0=indices(i);
			b_w->SetElement(i, -(intgrad(i0).sup()-grad(i0)));
			b_z->SetElement(i, (intgrad(i0).inf()-grad(i0)));
		}

		if (!(ct&Func::CONVEX)) { // nonconvex function
			Pointer<UserVector<double> > b_x(new dvector(grad)); // \nabla f(\hat x)
			double constant=val;
			for (int i=0, j=0; i<x.dim(); ++i) {
				switch(omit[i]) {
					case 0:
						constant-=grad(i)*x_ref(i);
						++j;
						break;
					case 1: // w_i fixed to 0
						b_x->SetElement(i, intgrad(i).inf());
						constant-=intgrad(i).inf()*x_ref(i);
						break;
					case 2: // z_i fixed to 0
					case 3: // linear variable
						b_x->SetElement(i, intgrad(i).sup());
						constant-=intgrad(i).sup()*x_ref(i);
						break;
				}
			}

			cut->coninfos.push_back(IntervalGradientCut::ConstraintInfo(b_x, b_w, b_z, constant));
			cut->coninfos_size++;
//		out_log << prob->con_names[c] << ":\t b_x: " << *b_x << "b_w: " << *b_w << "b_z: " << *b_z << "c: " << constant << endl;
		}

		if (prob->con_eq[c] && !(ct&Func::CONCAVE)) { // nonconcave function
			Pointer<UserVector<double> > b_x=new dvector(grad); *b_x*=-1;
			double constant=-val;
			for (int i=0; i<x.dim(); ++i) {
				switch(omit[i]) {
					case 0:
						constant+=grad(i)*x_ref(i);
						break;
					case 1: // w_i fixed to 0
						b_x->SetElement(i, -intgrad(i).sup());
						constant+=intgrad(i).sup()*x_ref(i);
						break;
					case 2: // z_i fixed to 0
					case 3: // linear variable
						b_x->SetElement(i, -intgrad(i).inf());
						constant+=intgrad(i).inf()*x_ref(i);
						break;
				}
			}
			cut->coninfos.push_back(IntervalGradientCut::ConstraintInfo(b_x, b_z, b_w, constant));
			cut->coninfos_size++;
//		out_log << '-' << prob->con_names[c] << ":\t b_x: " << *b_x << "b_w: " << *b_w << "b_z: " << *b_z << "c: " << constant << endl;
		}
	}
	if (cut->coninfos.empty()) cut=NULL;
#endif
	return cut;
}

// ---------------------------------------- SimpleCut ----------------------------------------

void SimpleCut::scale(double max_coeff) {
	double infnorm=max_coeff; // well, this is cheating, but we are interested in the inf norm only if its larger than max_coeff 
	for (int i=0; i<coeff->dim(); ++i) {
		double coeff_i=(*coeff)(i);
		if (fabs(coeff_i)>infnorm) infnorm=fabs(coeff_i);
	}
	if (infnorm>max_coeff) {
		*coeff/=infnorm;
		constant/=infnorm;
	}
}

ostream& operator<<(ostream& out, const SimpleCut& cut) {
	if (cut.coeff)
		for (int i=0; i<cut.coeff->dim(); ++i)
			if ((*cut.coeff)(i)) out << " +" << 2*(*cut.coeff)(i) << "*x" << i;
	out << " + " << cut.constant << " <=0" << endl;
	return out;	
}

// ---------------------------------------- Cut ----------------------------------------------

template <class CutType>
Cut<CutType>::~Cut() { }

template <class CutType>
void Cut<CutType>::duplicate_nodeinfo(Pointer<MinlpNode> oldnode, Pointer<MinlpNode> newnode) {
	if (global) return;
	typename map<Pointer<MinlpNode>, NodeInfo>::iterator it(nodes.find(oldnode));
	if (it!=nodes.end()) nodes.insert(pair<Pointer<MinlpNode>, NodeInfo>(newnode, NodeInfo(it->second)));
}

template <class CutType>
bool Cut<CutType>::remove_node(Pointer<MinlpNode> node) {
	if (global) return false;
	typename map<Pointer<MinlpNode>, NodeInfo>::iterator it(nodes.find(node));
	if (it!=nodes.end()) nodes.erase(it);
	return nodes.empty();
}

template <class CutType>
pair<bool, bool> Cut<CutType>::set_inactivetime(Pointer<MinlpNode> node, bool increase, int limit) {
	pair<bool, bool> ret(false, false);
	if (global) {
		if (increase) {
			if (++inactive_time_global>limit) ret.first=ret.second=true;
		} else
			inactive_time_global=0;
		return ret;
	}
	// rest is for local cuts
	
	typename map<Pointer<MinlpNode>, NodeInfo>::iterator it(nodes.find(node));
	assert(it!=nodes.end());
//	if (it==nodes.end()) return ret; // cut doesn't belong to this node.

	if (increase) { // remove local cuts immediately
		if (++it->second.inactive_time>3) {
			ret.second=true;
			nodes.erase(it);
			ret.first=nodes.empty();
			return ret;
		}
	} else
		it->second.inactive_time=0;
	return ret;
}

// -------------------------------------- CutPool ------------------------------------------

pair<bool, bool> CutPool::CutInfo::set_inactivetime(Pointer<MinlpNode> node, int inactivetime_limit) {
	switch (type) {
		case SIMPLE: return it_simplecuts->set_inactivetime(node, inactive, inactivetime_limit);
		case LINEARIZATION: return it_linearizationcuts->set_inactivetime(node, inactive, inactivetime_limit);
		case INTERVALGRADIENT: return it_intgradcuts->set_inactivetime(node, inactive, inactivetime_limit);
	}
	out_err << "CutType unknown. Aborting." << endl;
	exit(-1);
	return pair<bool,bool>(false, false);
}

CutPool::CutPool(const CutPool& cutpool, Pointer<MinlpNode> node)
: inactivetime_limit_global(cutpool.inactivetime_limit_global), inactivetime_limit_local(cutpool.inactivetime_limit_local), cuts_size(0)
{	for (list<Cut<SimpleCut> >::const_iterator it(cutpool.simplecuts.begin()); it!=cutpool.simplecuts.end(); ++it)
		if (it->valid(node)) {
			simplecuts.push_back(Cut<SimpleCut>(it->cut));
			simplecuts.back().tagged=true;
			++cuts_size;
		}
	for (list<Cut<LinearizationCut> >::const_iterator it(cutpool.linearizationcuts.begin()); it!=cutpool.linearizationcuts.end(); ++it)
		if (it->valid(node)) {
			linearizationcuts.push_back(Cut<LinearizationCut>(it->cut));
			linearizationcuts.back().tagged=true;
			++cuts_size;
		}
	for (list<Cut<IntervalGradientCut> >::const_iterator it(cutpool.intgradcuts.begin()); it!=cutpool.intgradcuts.end(); ++it)
		if (it->valid(node)) {
			intgradcuts.push_back(Cut<IntervalGradientCut>(it->cut));
			intgradcuts.back().tagged=true;
			++cuts_size;
		}
}

CutPool::CutInfo CutPool::add_cut(Pointer<SimpleCut> simplecut, Pointer<MinlpNode> node, int blocknr) {
	simplecuts.push_back(Cut<SimpleCut>(simplecut, node));
	++cuts_size;
	return CutInfo(--simplecuts.end(), blocknr);
}

CutPool::CutInfo CutPool::add_cut(Pointer<LinearizationCut> linearizationcut, Pointer<MinlpNode> node, int blocknr) {
	linearizationcuts.push_back(Cut<LinearizationCut>(linearizationcut, node));
	++cuts_size;
	return CutInfo(--linearizationcuts.end(), blocknr);
}

CutPool::CutInfo CutPool::add_cut(Pointer<IntervalGradientCut> intgradcut, Pointer<MinlpNode> node, int blocknr) {
	intgradcuts.push_back(Cut<IntervalGradientCut>(intgradcut, node));
	++cuts_size;
	return CutInfo(--intgradcuts.end(), blocknr);
}

void CutPool::integrate(const CutPool& cutpool, Pointer<MinlpNode> node) {
	for (list<Cut<SimpleCut> >::const_iterator it(cutpool.simplecuts.begin()); it!=cutpool.simplecuts.end(); ++it)
		if (it->global && (!it->tagged)) {
			simplecuts.push_back(Cut<SimpleCut>(it->cut, node));
			++cuts_size;
		}
	for (list<Cut<LinearizationCut> >::const_iterator it(cutpool.linearizationcuts.begin()); it!=cutpool.linearizationcuts.end(); ++it)
		if (it->global && (!it->tagged)) {
			linearizationcuts.push_back(Cut<LinearizationCut>(it->cut, node));
			++cuts_size;
		}
	for (list<Cut<IntervalGradientCut> >::const_iterator it(cutpool.intgradcuts.begin()); it!=cutpool.intgradcuts.end(); ++it)
		if (it->global && (!it->tagged)) {
			intgradcuts.push_back(Cut<IntervalGradientCut>(it->cut, node));
			++cuts_size;
		}
}

void CutPool::duplicate_nodeinfo(Pointer<MinlpNode> oldnode, Pointer<MinlpNode> newnode) {
	for (list<Cut<SimpleCut> >::iterator it(simplecuts.begin()); it!=simplecuts.end(); ++it)
		it->duplicate_nodeinfo(oldnode, newnode);
	for (list<Cut<LinearizationCut> >::iterator it(linearizationcuts.begin()); it!=linearizationcuts.end(); ++it)
		it->duplicate_nodeinfo(oldnode, newnode);
	for (list<Cut<IntervalGradientCut> >::iterator it(intgradcuts.begin()); it!=intgradcuts.end(); ++it)
		it->duplicate_nodeinfo(oldnode, newnode);
}

void CutPool::remove_node(Pointer<MinlpNode> node) {
	list<Cut<SimpleCut> >::iterator it(simplecuts.begin());
	while (it!=simplecuts.end())
		if (it->remove_node(node)) {
			it=simplecuts.erase(it);
			--cuts_size;
		}
		else ++it;

	list<Cut<LinearizationCut> >::iterator it3(linearizationcuts.begin());
	while (it3!=linearizationcuts.end())
		if (it3->remove_node(node)) {
			it3=linearizationcuts.erase(it3);
			--cuts_size;
		}
		else ++it3;


	list<Cut<IntervalGradientCut> >::iterator it2(intgradcuts.begin());
	while (it2!=intgradcuts.end())
		if (it2->remove_node(node)) {
			it2=intgradcuts.erase(it2);
			--cuts_size;
		}
		else ++it2;
}

void CutPool::update_nodeinfo(list<CutInfo>& cutinfos, int block_nr, Pointer<MinlpNode> node) {
	for (list<CutInfo>::iterator it(cutinfos.begin()); it!=cutinfos.end(); it++) {
		if (it->block_nr!=block_nr) continue;

		int inactivetime_limit;
		switch(it->type) {
			case SIMPLE: inactivetime_limit=it->it_simplecuts->global ? inactivetime_limit_global : inactivetime_limit_local; break;
			case LINEARIZATION: inactivetime_limit=it->it_linearizationcuts->global ? inactivetime_limit_global : inactivetime_limit_local; break;
			case INTERVALGRADIENT: inactivetime_limit=it->it_intgradcuts->global ? inactivetime_limit_global : inactivetime_limit_local; break;
			default: out_err << "CutType unknown. Aborting." << endl; exit(-1);
		}
		pair<bool, bool> ret(it->set_inactivetime(node, inactivetime_limit));

		it->removed=ret.second; // node invalid for the cut now?
		if (ret.first) {  // cut can be removed
			--cuts_size;
			switch (it->type) {
				case SIMPLE: simplecuts.erase(it->it_simplecuts); break;
				case LINEARIZATION: linearizationcuts.erase(it->it_linearizationcuts); break;
				case INTERVALGRADIENT: intgradcuts.erase(it->it_intgradcuts); break;
			}
		}
	}
}

void CutPool::get_cuts(list<CutInfo>& cutinfos, int blocknr, Pointer<MinlpNode> node) {
	for (list<Cut<SimpleCut> >::iterator it(simplecuts.begin()); it!=simplecuts.end(); ++it)
		if (it->valid(node))
			cutinfos.push_back(CutInfo(it, blocknr));
	for (list<Cut<LinearizationCut> >::iterator it(linearizationcuts.begin()); it!=linearizationcuts.end(); ++it)
		if (it->valid(node))
			cutinfos.push_back(CutInfo(it, blocknr));
	for (list<Cut<IntervalGradientCut> >::iterator it(intgradcuts.begin()); it!=intgradcuts.end(); ++it)
		if (it->valid(node))
			cutinfos.push_back(CutInfo(it, blocknr));
}

void CutPool::update_cuts(Pointer<MinlpNode> node, int blocknr, const dvector& low, const dvector& up, IntervalGradientCutGenerator& generator, LinearizedConCutGenerator& linconcutgen) {
	list<Pointer<IntervalGradientCut> > newcuts;
	int oldcuts_nr=0;
	for (list<Cut<IntervalGradientCut> >::iterator it(intgradcuts.begin()); it!=intgradcuts.end();)
		if (it->valid(node)) {
			oldcuts_nr+=it->cut->coninfos_size;
			Pointer<IntervalGradientCut> newcut=generator.update_cuts(it->cut, blocknr, low, up);
			if (it->remove_node(node)) {
				it=intgradcuts.erase(it);
				--cuts_size;
			} else ++it;
			if (newcut) newcuts.push_back(newcut);
		} else ++it;

	if (oldcuts_nr) {
		int cuts_nr=0;
		for (list<Pointer<IntervalGradientCut> >::iterator it(newcuts.begin()); it!=newcuts.end(); ++it) {
			cuts_nr+=(*it)->coninfos_size;
			add_cut(*it, node, blocknr);
		}
//		out_log << "Replaced " << oldcuts_nr << " IntervalGradientCuts by " << cuts_nr << " updated ones." << endl;
	}
	newcuts.clear();
	
	list<Pointer<LinearizationCut> > newcuts2;
	oldcuts_nr=0;
	for (list<Cut<LinearizationCut> >::iterator it(linearizationcuts.begin()); it!=linearizationcuts.end();)
		if (it->valid(node) && it->cut->objcon_nr>-2) {
			++oldcuts_nr;
			Pointer<LinearizationCut> newcut=linconcutgen.update_cut(*it->cut, blocknr, low, up);
			if (it->remove_node(node)) {
				it=linearizationcuts.erase(it);
				--cuts_size;
			} else ++it;
			if (newcut) newcuts2.push_back(newcut);			
		} else ++it;
		
	if (oldcuts_nr) {
		for (list<Pointer<LinearizationCut> >::iterator it(newcuts2.begin()); it!=newcuts2.end(); ++it)
			add_cut(*it, node, blocknr);
//		out_log << "Replaced " << oldcuts_nr << " LinearizationCuts by " << newcuts2.size() << " updated ones." << endl;
	}
}

bool CutPool::feasible(Pointer<MinlpNode> node, const dvector& x, double tol) const {
	for (list<Cut<SimpleCut> >::const_iterator it(simplecuts.begin()); it!=simplecuts.end(); ++it) {
		if (!it->valid(node)) continue; // cut does not belong to this node
		if (!it->cut->feasible(x, tol)) return false;
	}
	for (list<Cut<LinearizationCut> >::const_iterator it(linearizationcuts.begin()); it!=linearizationcuts.end(); ++it) {
		if (!it->valid(node)) continue; // cut does not belong to this node
		if (!it->cut->feasible(x, tol)) return false;
	}
	if (!intgradcuts.empty() && node) out_log << "CutPool::feasible: Warning, check of satisfaction of intervalgradient cuts not implemented currently." << endl;
//	assert(intgradcuts.empty()); // TODO
	return true;
}

int CutPool::nr_local_cuts(Pointer<MinlpNode> node) const {
	int count=0;
	for (list<Cut<SimpleCut> >::const_iterator it(simplecuts.begin()); it!=simplecuts.end(); ++it)
		if ((!it->global) && it->nodes.count(node)) count++;
	for (list<Cut<LinearizationCut> >::const_iterator it(linearizationcuts.begin()); it!=linearizationcuts.end(); ++it)
		if ((!it->global) && it->nodes.count(node)) count++;
	for (list<Cut<IntervalGradientCut> >::const_iterator it(intgradcuts.begin()); it!=intgradcuts.end(); ++it)
		if ((!it->global) && it->nodes.count(node)) count++;
	return count;
}

int CutPool::nr_global_cuts() const {
	int count=0;
	for (list<Cut<SimpleCut> >::const_iterator it(simplecuts.begin()); it!=simplecuts.end(); it++)
		if (it->global) count++;
	for (list<Cut<LinearizationCut> >::const_iterator it(linearizationcuts.begin()); it!=linearizationcuts.end(); ++it)
		if (it->global) count++;
	for (list<Cut<IntervalGradientCut> >::const_iterator it(intgradcuts.begin()); it!=intgradcuts.end(); ++it)
		if (it->global) count++;
	return count;
}

// ------------------------------------ LinearizedConCut -------------------------------------------------

const double LinearizedConCutGenerator::tol=1E-4;

LinearizedConCutGenerator::LinearizedConCutGenerator(Pointer<MinlpProblem> prob_, Pointer<MINLPData> minlpdata_, Pointer<Reformulation> reform_)
: prob(prob_), minlpdata(minlpdata_), reform(reform_)
{ }


LinearizationCut LinearizedConCutGenerator::get_cut(const dvector& x, int c, int block_nr, double val) {
	if (val==INFINITY) val=prob->con[c]->eval(x);
//	out_log << "value of quadr. approx.: " << val << endl;

	LinearizationCut cut;
	if (!finite(val)) return cut;
	
	cut.constant=val;
	if (prob->con[c]->A[block_nr] || prob->con[c]->s[block_nr]) {
		cut.coeff=new SparseVector<double>(prob->block[block_nr].size());
		*cut.coeff=prob->con[c]->grad(x, block_nr);
		*cut.coeff*=.5; // :-(
	} else {
		assert(prob->con[c]->b[block_nr]);
		cut.coeff=prob->con[c]->b[block_nr]->getcopy();
	}
	cut.constant-=2*(*cut.coeff*x);

	cut.scale();
	return cut;
}

int LinearizedConCutGenerator::get_cuts(list<pair<LinearizationCut, pair<int, bool> > >& cuts, const MINLPData::ObjCon& objcon, int objcon_nr, const dvector& x, const dvector& lower, const dvector& upper, bool violated_polyest_only) {
	int nr=0;
	for (int k=0; k<objcon.func->block.size(); ++k) {
		if (objcon.func->type(k)==SepQcFunc::LINEAR || objcon.func->type(k)==SepQcFunc::CONSTANT) continue;
		
		dvector xk(x, objcon.func->block[k]);
		double val=objcon.func->eval(xk, k);
		if (!finite(val)) {
			out_err << "Generation of LinearizedConCut for constraint " << objcon.name << " failed due to function evaluation error." << endl;
			continue;
		}

		map<int,int>::const_iterator reformit(objcon.reformulation_constraints_lower.find(k));
		int tindex=(reformit==objcon.reformulation_constraints_lower.end() ? -1 : reform->related_t[reformit->second]);
		if (tindex<0) val+=objcon.func->c; // hasn't been reformulated

		// not active
		if (tindex<0 && val<-tol) continue;
		if (tindex>=0 && val-x(prob->block[k][tindex])<-tol) continue;

		const SepQcFunc& quadfunc(objcon.polynomial_underestimator ? *objcon.polynomial_underestimator : *objcon.func);
		if (objcon.polynomial_underestimator) {
			val=quadfunc.eval(xk, k);
			if (tindex<0) val+=quadfunc.c;
			else {
				map<int,double>::const_iterator approxit(objcon.polynomial_approx_constants_lower.find(k));
				if (approxit!=objcon.polynomial_approx_constants_lower.end()) val+=approxit->second; //add constant part introduced by computing quadratic underestimator
			}
		}
		if (violated_polyest_only) { // skip constraint where poly.uest. not violated
			if (tindex<0 && val<-tol) continue;
			if (tindex>=0 && val-x(prob->block[k][tindex])<-tol) continue;
		}
//		if (tindex>=0) out_log << "value of quadr. approx.: " << val-x(prob->block[k][tindex]) << endl;

		LinearizationCut cut(new SparseVector<double>(quadfunc.grad(xk, k)), val, objcon_nr, true, xk, tindex>=0 ? x(prob->block[k][tindex]) : 0);
		*cut.coeff*=.5;

		// and now the convexification....
		for (int i=0; i<cut.coeff->size(); ++i) {
			int i0=objcon.func->block[k][i];
			if (lower(i0)==upper(i0)) continue;
			map<int,double>::const_iterator shiftit(objcon.convexification_characteristica_lower.find(i0));
			if (shiftit==objcon.convexification_characteristica_lower.end()) continue;
			double shift=shiftit->second;
			if (!shift) continue;
			(*cut.coeff)[i]+=shift*(xk(i)-.5*lower(i0)-.5*upper(i0));
			cut.constant+=(xk(i)-lower(i0))*shift*(xk(i)-upper(i0));
		}
		
		double violation=cut.constant;
		if (tindex>=0) violation-=x(prob->block[k][tindex]);
		if (violation<tol && violated_polyest_only) continue; // skip constraints where the convex relax. is not violated

		double maxcoeff=1;
		for (int i=0; i<cut.coeff->size(); ++i) if (2*fabs((*cut.coeff)(i))>maxcoeff) maxcoeff=2*fabs((*cut.coeff)(i));
//		out_log << "Constraint " << objcon.name << " violated by " << violation << '\t' << maxcoeff << '\t' << violation/maxcoeff;
		violation/=maxcoeff;
		if (violation>max_violation) max_violation=violation;
		
		cut.constant-=2*(*cut.coeff*xk);
		cut.coeff->resize(prob->block[k].size()); // resize to extended problem size

		// if had been reformulated, add -t_k variable to gradient
		if (tindex>=0) cut.coeff->SetElement(tindex, -.5);

//		out_log << objcon.name << " block " << k << " const.= " << cut.constant << "\t coeff=" << *cut.coeff;
		cut.scale();
		cuts.push_back(pair<LinearizationCut, pair<int, bool> >(cut, pair<int, bool>(k, objcon.func->get_curvature(k)==Func::CONVEX)));
		++nr;

//		double valcheck=2* (*cut.coeff * x(prob->block[k]))+cut.constant;
//		out_log << '\t' << valcheck;
//		out_log << endl;
	}

	return nr;	
}

int LinearizedConCutGenerator::get_cuts(list<pair<LinearizationCut, pair<int, bool> > >& cuts, const MINLPData::Constraint& con, int objcon_nr, const dvector& x, const dvector& lower, const dvector& upper, bool violated_polyest_only) {
	int nr=get_cuts(cuts, (const MINLPData::ObjCon&)con, objcon_nr, x, lower, upper, violated_polyest_only);
	if (!con.equality) return nr;
	
	// cut for upper part of equality constraint
	for (int k=0; k<con.func->block.size(); ++k) {
		if (con.func->type(k)==SepQcFunc::LINEAR || con.func->type(k)==SepQcFunc::CONSTANT) continue;
		dvector xk(x, con.func->block[k]); // this skips the reformulation variables
		double val=-con.func->eval(xk, k);
		if (!finite(val)) {
			out_err << "Generation of LinearizedConCut for constraint -" << con.name << " failed due to function evaluation error." << endl;
			continue;
		}

		map<int,int>::const_iterator reformit(con.reformulation_constraints_upper.find(k));
		int tindex=(reformit==con.reformulation_constraints_upper.end() ? -1 : reform->related_t[reformit->second]);
		if (tindex<0) val-=con.func->c; // hasn't been reformulated

		// not active
		if (tindex<0 && val<-tol) continue;
		if (tindex>=0 && val+x(prob->block[k][tindex])<-tol) continue;
		
		Pointer<UserVector<double> > coeff;
//		const SepQcFunc& quadfunc(con.polynomial_underestimator ? *con.polynomial_overestimator : *con.func);
		if (con.polynomial_overestimator) {
			val=con.polynomial_overestimator->eval(xk, k);
			if (tindex<0) val+=con.polynomial_overestimator->c;
			else {
				map<int,double>::const_iterator approxit(con.polynomial_approx_constants_upper.find(k));
				if (approxit!=con.polynomial_approx_constants_upper.end()) val+=approxit->second; //add constant part introduced by computing quadratic underestimator
			}

			if (violated_polyest_only) { // skip constraint where poly.oest. not violated
				if (tindex<0 && val<-tol) continue;
				if (tindex>=0 && val+x(prob->block[k][tindex])<-tol) continue;
			}
			
			coeff=new SparseVector<double>(con.polynomial_overestimator->grad(xk, k));
			*coeff*=.5;
		} else {
			coeff=new SparseVector<double>(con.func->grad(xk, k));
			*coeff*=-.5;
		}

		LinearizationCut cut(coeff, val, objcon_nr, false, xk, tindex>=0 ? x(prob->block[k][tindex]) : 0);

		// and now the convexification....
		for (int i=0; i<cut.coeff->size(); ++i) {
			int i0=con.func->block[k][i];
			if (lower(i0)==upper(i0)) continue;
			map<int,double>::const_iterator shiftit(con.convexification_characteristica_upper.find(i0));
			if (shiftit==con.convexification_characteristica_upper.end()) continue;
			double shift=shiftit->second;
			if (!shift) continue;
			(*cut.coeff)[i]+=shift*(xk(i)-.5*lower(i0)-.5*upper(i0));
			cut.constant+=(xk(i)-lower(i0))*shift*(xk(i)-upper(i0));
		}

		double violation=cut.constant;
		if (tindex>=0) violation+=x(prob->block[k][tindex]);

		if (violation<tol && violated_polyest_only) continue; // skip constraints where the convex relax. is not violated

		double maxcoeff=1;
		for (int i=0; i<cut.coeff->size(); ++i) if (2*fabs((*cut.coeff)(i))>maxcoeff) maxcoeff=2*fabs((*cut.coeff)(i));
//		out_log << "Constraint " << con.name << " violated by " << violation << '\t' << maxcoeff << '\t' << violation/maxcoeff;
		violation/=maxcoeff;
		if (violation>max_violation) max_violation=violation;

		cut.constant-=2*(*cut.coeff*xk);
		cut.coeff->resize(prob->block[k].size()); // resize to extended problem size

		// if had been reformulated, add t_k variable
		if (tindex>=0) cut.coeff->SetElement(tindex, .5);
		
//		out_log << '-' << con.name << " block " << k  << "\t const.= " << cut.constant << "\t coeff=" << *cut.coeff;
		cut.scale();
		cuts.push_back(pair<LinearizationCut, pair<int, bool> >(cut, pair<int, bool>(k, con.func->get_curvature(k)==Func::CONCAVE)));
		++nr;

//		double valcheck=2* (*cut.coeff * x(prob->block[k]))+cut.constant;
//		out_log << '\t' << valcheck;
//		out_log << endl;
	}
	
	return nr;
}


int LinearizedConCutGenerator::get_cuts(list<pair<LinearizationCut, pair<int, bool> > >& cuts, const dvector& x, const dvector& lower, const dvector& upper, bool violated_polyest_only) {
	int nr=0;
//	out_log << "Generating cuts for objective: " << endl;
	nr+=get_cuts(cuts, minlpdata->obj, -1, x, lower, upper, violated_polyest_only);
	
	for (int c=0; c<minlpdata->con.size(); ++c) {
		if (minlpdata->con[c].functype==SepQcFunc::QUADRATIC || minlpdata->con[c].functype==SepQcFunc::NONQUAD) {
//			out_log << "Generating cuts for constraint " << minlpdata->con[c].name << ':' << endl;
			nr+=get_cuts(cuts, minlpdata->con[c], c, x, lower, upper, violated_polyest_only);
		}
	}
	return nr;
}


Pointer<LinearizationCut> LinearizedConCutGenerator::update_cut(LinearizationCut& cut, unsigned int block_nr, const dvector& lower, const dvector& upper) {
//	return new LinearizationCut(cut);
	
	Pointer<LinearizationCut> newcut=new LinearizationCut(NULL, 0., cut.objcon_nr, cut.derived_from_lower_part, cut.x, cut.t_value);
	for (int i=0; i<newcut->x.dim(); ++i) {
		if (newcut->x(i)<lower(i)) newcut->x[i]=lower(i);
		else if (newcut->x(i)>upper(i)) newcut->x[i]=upper(i);
	}

	if (cut.derived_from_lower_part) {
		const MINLPData::ObjCon& objcon(cut.objcon_nr==-1 ? (MINLPData::ObjCon&)minlpdata->obj : (MINLPData::ObjCon&)minlpdata->con[cut.objcon_nr]);
	
		map<int,int>::const_iterator reformit(objcon.reformulation_constraints_lower.find(block_nr));
		int tindex=(reformit==objcon.reformulation_constraints_lower.end() ? -1 : reform->related_t[reformit->second]);
		if (tindex>0) {
			if (newcut->t_value<lower(tindex)) newcut->t_value=lower(tindex);
			else if (newcut->t_value>upper(tindex)) newcut->t_value=upper(tindex);
		}

		double val=objcon.func->eval(newcut->x, block_nr);
		if (!finite(val)) {
			out_err << "Update of LinearizedConCut failed." << endl;
			return NULL;
		}
		if (tindex<0) val+=objcon.func->c; // hasn't been reformulated
	
		// not active
		if (tindex<0 && val<-tol) return NULL;
		if (tindex>=0 && val-newcut->t_value<-tol) return NULL;
	
		const SepQcFunc& quadfunc(objcon.polynomial_underestimator ? *objcon.polynomial_underestimator : *objcon.func);
		if (objcon.polynomial_underestimator) {
			val=quadfunc.eval(newcut->x, block_nr);
			if (tindex<0) val+=quadfunc.c;
			else {
				map<int,double>::const_iterator approxit(objcon.polynomial_approx_constants_lower.find(block_nr));
				if (approxit!=objcon.polynomial_approx_constants_lower.end()) val+=approxit->second; //add constant part introduced by computing quadratic underestimator
			}
		}
	
		newcut->coeff=new SparseVector<double>(quadfunc.grad(newcut->x, block_nr));
		*newcut->coeff*=.5;
		newcut->constant=val;
	
		// and now the convexification....
		for (int i=0; i<newcut->coeff->size(); ++i) {
			if (lower(i)==upper(i)) continue;
			map<int,double>::const_iterator shiftit(objcon.convexification_characteristica_lower.find(objcon.func->block[block_nr][i]));
			if (shiftit==objcon.convexification_characteristica_lower.end()) continue;
			double shift=shiftit->second;
			if (!shift) continue;
			(*newcut->coeff)[i]+=shift*(newcut->x(i)-.5*lower(i)-.5*upper(i));
			newcut->constant+=(newcut->x(i)-lower(i))*shift*(newcut->x(i)-upper(i));
		}
		
		double violation=newcut->constant;
		if (tindex>=0) violation-=newcut->t_value;
		if (violation<tol) return NULL;		
		
		newcut->constant-=2*(*newcut->coeff*newcut->x);
	
		newcut->coeff->resize(cut.coeff->size()); // resize to extended problem size
	
		// if had been reformulated, add -t_k variable to gradient
		if (tindex>=0) newcut->coeff->SetElement(tindex, -.5);
	} else {
		MINLPData::Constraint& con(minlpdata->con[cut.objcon_nr]);
	
		map<int,int>::const_iterator reformit(con.reformulation_constraints_upper.find(block_nr));
		int tindex=(reformit==con.reformulation_constraints_upper.end() ? -1 : reform->related_t[reformit->second]);
		if (tindex>0) {
			if (newcut->t_value<lower(tindex)) newcut->t_value=lower(tindex);
			else if (newcut->t_value>upper(tindex)) newcut->t_value=upper(tindex);
		}

		double val=-con.func->eval(newcut->x, block_nr);
		if (!finite(val)) {
			out_err << "Update of LinearizedConCut failed." << endl;
			return NULL;
		}
		if (tindex<0) val-=con.func->c; // hasn't been reformulated
	
		// not active
		if (tindex<0 && val<-tol) return NULL;
		if (tindex>=0 && val+newcut->t_value) return NULL;
	
		if (con.polynomial_overestimator) {
			val=con.polynomial_overestimator->eval(newcut->x, block_nr);
			if (tindex<0) val+=con.polynomial_overestimator->c;
			else {
				map<int,double>::const_iterator approxit(con.polynomial_approx_constants_upper.find(block_nr));
				if (approxit!=con.polynomial_approx_constants_upper.end()) val+=approxit->second; //add constant part introduced by computing quadratic underestimator
			}
			
			newcut->coeff=new SparseVector<double>(con.polynomial_overestimator->grad(newcut->x, block_nr));
			*newcut->coeff*=.5;
		} else {
			newcut->coeff=new SparseVector<double>(con.func->grad(newcut->x, block_nr));
			*newcut->coeff*=-.5;
		}
		newcut->constant=val;
	
		// and now the convexification....
		for (int i=0; i<newcut->coeff->size(); ++i) {
			if (lower(i)==upper(i)) continue;
			map<int,double>::const_iterator shiftit(con.convexification_characteristica_upper.find(con.func->block[block_nr][i]));
			if (shiftit==con.convexification_characteristica_upper.end()) continue;
			double shift=shiftit->second;
			if (!shift) continue;
			(*newcut->coeff)[i]+=shift*(newcut->x(i)-.5*lower(i)-.5*upper(i));
			newcut->constant+=(newcut->x(i)-lower(i))*shift*(newcut->x(i)-upper(i));
		}
		
		double violation=newcut->constant;
		if (tindex>=0) violation+=newcut->t_value;
		if (violation<tol) return NULL;
				
		newcut->constant-=2*(*newcut->coeff*newcut->x);
	
		newcut->coeff->resize(cut.coeff->size()); // resize to extended problem size
	
		// if had been reformulated, add t_k variable
		if (tindex>=0) newcut->coeff->SetElement(tindex, .5);
	}
	//		out_log << objcon.name << " block " << k << " const.= " << cut.constant << "\t coeff=" << *cut.coeff;
	if (newcut) newcut->scale();
	return newcut;
}

void LinearizedConCutGenerator::set_reform(Pointer<Reformulation> reform_) { reform=reform_; }
