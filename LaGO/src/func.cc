// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "func.h"

// -------------------------------- SparsityInfo -----------------------------------------

void SparsityInfo::add_quadratic(int index, double coeff_lin, double coeff_quad) {
	pair<map<int, QuadraticVariable>::iterator, bool> ret(quadratic->insert(pair<int, QuadraticVariable>(index, QuadraticVariable(coeff_lin, coeff_quad))));
	if (!ret.second) {
		ret.first->second.coeff_lin+=coeff_lin;
		ret.first->second.coeff_quad+=coeff_quad;
	}
}

void SparsityInfo::add_sparsity_pattern(int index1, int index2, double coeff) {
	if (index2<index1) add_sparsity_pattern(index2, index1, coeff);
	pair<map<pair<int, int>, NonlinearConnection>::iterator, bool> ret(sparsity_pattern->insert(pair<pair<int, int>, NonlinearConnection>(pair<int, int>(index1, index2), NonlinearConnection(coeff))));
	if (!ret.second) ret.first->second.coeff+=coeff;
}

SparsityInfo::SparsityInfo(int init_level) {
	if (init_level>=1) {
		linear=new map<int, LinearVariable>();
		nonlinear=new map<int, NonlinearVariable>();
		if (init_level>=2) {
			quadratic=new map<int, QuadraticVariable>();
			sparsity_pattern=new map<pair<int, int>, NonlinearConnection>();
		}
	}
}

SparsityInfo::SparsityInfo(const SparsityInfo& si) {
	if (si.linear) linear=new map<int, LinearVariable>(*si.linear);
	if (si.nonlinear) nonlinear=new map<int, NonlinearVariable>(*si.nonlinear);
	if (si.quadratic) quadratic=new map<int, QuadraticVariable>(*si.quadratic);
	if (si.sparsity_pattern) sparsity_pattern=new map<pair<int, int>, NonlinearConnection>(*si.sparsity_pattern);
}

SparsityInfo::SparsityInfo(const SparsityInfo& si1, const SparsityInfo& si2) {
	if (si1.linear) linear=new map<int, LinearVariable>(*si1.linear);
	if (si1.nonlinear) nonlinear=new map<int, NonlinearVariable>(*si1.nonlinear);
	if (si1.quadratic) quadratic=new map<int, QuadraticVariable>(*si1.quadratic);
	if (si1.sparsity_pattern) sparsity_pattern=new map<pair<int, int>, NonlinearConnection>(*si1.sparsity_pattern);
	add(si2);
}

ostream& operator<<(ostream& out, const SparsityInfo& si) {
	if (si.linear) {
		out << "Linear: ";
		for (map<int, SparsityInfo::LinearVariable>::iterator it(si.linear->begin()); it!=si.linear->end(); ++it)
			out << it->first << '(' << it->second << ')' << ' ';
		out << endl;
	}
	if (si.nonlinear) {
		out << "Nonlinear: ";
		for (map<int, SparsityInfo::NonlinearVariable>::iterator it(si.nonlinear->begin()); it!=si.nonlinear->end(); ++it)
			out << it->first << ' ';
		out << endl;
	}
	if (si.quadratic) {
		out << "Quadratic: ";
		for (map<int, SparsityInfo::QuadraticVariable>::iterator it(si.quadratic->begin()); it!=si.quadratic->end(); ++it)
			out << it->first << '(' << it->second << ')' << ' ';
		out << endl;
	}
	if (si.sparsity_pattern) {
		out << "Sparsity pattern: ";
		for (map<pair<int,int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it)
			out << '(' << it->first.first << ',' << it->first.second << ';' << it->second << ')' << ' ';
		out << endl;
	}

	return out;
}

bool SparsityInfo::compute_sparsity_pattern(const Func& f, const vector<dvector>& sample_set) {
	// we need more than one sample point, if s-function is set and we want to know, which variables are quadratic and nonlinear/nonquadratic.
	assert(sample_set.size()>1);

	quadratic=new map<int, SparsityInfo::QuadraticVariable>();
	sparsity_pattern=new map<pair<int, int>, SparsityInfo::NonlinearConnection>();
	bool sign_changed=false;

	Pointer<dvector> hm, hm_old;
	SparseVector<double> e(f.dim());
	for (map<int, SparsityInfo::NonlinearVariable>::iterator it_index(nonlinear->begin()); it_index!=nonlinear->end(); ++it_index) {
		e.SetElement(it_index->first, 1., false);
		bool quad=true; // indicates, wheter it_index->first is quadratic
		double coeff_quad=0.;
		hm=NULL;
		for (vector<dvector>::const_iterator it_sp(sample_set.begin()); it_sp!=sample_set.end(); it_sp++) {
			hm_old=hm;
//			hm=new SparseVector<double>(func.block[block_nr].size());
			hm=new dvector(f.dim());

			f.HessMult(*hm, *it_sp, e); // computes it_index-th column of hessian for sample point it_sp
			for (int i=0; i<it_index->first; i++) { // todo: use SparseVector-properties here
				assert(finite((*hm)(i)));
				if ((*hm)(i)==0.) continue;
				if ((!hm_old) || (*hm_old)(i)==0.)  // appears first
					add_sparsity_pattern(i, it_index->first, (*hm)(i));
				if (hm_old && (*hm_old)(i)!=(*hm)(i)) {  // value changed since last run -> change to nonquadratic variable
					quad=false;
					quadratic->erase(i);
					sign_changed|=(hm_old && (*hm)(i)*(*hm_old)(i)<0);
				}
			}
			coeff_quad=(*hm)(it_index->first);
			if (coeff_quad && hm_old && coeff_quad!=(*hm_old)(it_index->first)) { // value changed since last run -> change to nonquadratic variable
				quad=false;
				sign_changed|=(hm_old && (*hm)(it_index->first)*(*hm_old)(it_index->first)<-rtol);
			}
		}

		e.DelElement(it_index->first);

		if (quad) add_quadratic(it_index->first, 0., coeff_quad/2.);
	}

	// make all quadratic variables, which are (indirectly) connected to nonquadratic variables, nonquadratic
	// this might be removed later
	set<int> to_check;
	for (map<int, SparsityInfo::NonlinearVariable>::iterator it(nonlinear->begin()); it!=nonlinear->end(); ++it)
		if (!quadratic->count(it->first)) to_check.insert(it->first);
	while (to_check.size()) {
		int i=*to_check.begin();
		to_check.erase(to_check.begin());
		map<pair<int,int>, SparsityInfo::NonlinearConnection>::iterator it(sparsity_pattern->begin());
		while (it!=sparsity_pattern->end() && (it->first.first<i)) {
			if (it->first.second==i) {
				map<int, SparsityInfo::QuadraticVariable>::iterator it2(quadratic->find(it->first.first));
				if (it2!=quadratic->end()) {
					quadratic->erase(it2);
					to_check.insert(it->first.first);
				}
			}
			++it;
		}
		while (it!=sparsity_pattern->end() && (it->first.first==i)) {
			map<int, SparsityInfo::QuadraticVariable>::iterator it2(quadratic->find(it->first.second));
			if (it2!=quadratic->end()) {
				quadratic->erase(it2);
				to_check.insert(it->first.second);
			}
			++it;
		}
	}

	if (!quadratic->empty()) { // linear coefficients for quadratic variables
		dvector try_point(sample_set.front());
		for (map<int, QuadraticVariable>::iterator it(quadratic->begin()); it!=quadratic->end(); ++it)
			try_point[it->first]=0.;
		dvector grad(try_point.dim());
		f.grad(grad, try_point);
		for (map<int, QuadraticVariable>::iterator it(quadratic->begin()); it!=quadratic->end(); ++it)
			it->second.coeff_lin=grad(it->first);
	}

  return sign_changed;
}

void SparsityInfo::add(const UserVector<double>& b) {
	double val;
	map<int, QuadraticVariable>::iterator it_q;
	for (int i=0; i<b.dim(); ++i) { // todo: use probable sparsity of b
		if (!(val=2*b(i))) continue;
		it_q=quadratic->find(i);
		if (it_q!=quadratic->end())
			it_q->second.coeff_lin+=val;
		else if (!nonlinear->count(i)) {
			pair<map<int, LinearVariable>::iterator, bool> ret(linear->insert(pair<int, LinearVariable>(i, val)));
			if (!ret.second) ret.first->second.coeff+=val;
		}
	}
}

void SparsityInfo::add(const UserVector<double>& b, const ivector& block) {
	double val;
	map<int, QuadraticVariable>::iterator it_q;
	for (int i=0; i<b.dim(); ++i) { // todo: use probable sparsity of b
		if (!(val=2*b(i))) continue;
		it_q=quadratic->find(block(i));
		if (it_q!=quadratic->end())
			it_q->second.coeff_lin+=val;
		else if (!nonlinear->count(block(i))) {
			pair<map<int, LinearVariable>::iterator, bool> ret(linear->insert(pair<int, LinearVariable>(block(i), val)));
			if (!ret.second) ret.first->second.coeff+=val;
		}
	}
}

void SparsityInfo::add(const UserMatrix& A) {
	const ExtUserMatrix* EA=NULL;
	const SparseMatrix2* SA=NULL;
	bool delete_SA=false;

	SA=dynamic_cast<const SparseMatrix2*>(&A);
	if (!SA) {
		EA=dynamic_cast<const ExtUserMatrix*>(&A);
		if (!EA) {
			SA=new SparseMatrix2(A);
			delete_SA=true;
		}
	}

	if (SA) {
		const int* col_ptr=SA->GetColPtr();
		const int* row_ind=SA->GetRowInd();
		for (int col=0; col<SA->dim(); col++) {
			if (col_ptr[col]==col_ptr[col+1]) continue; // empty column
			const int* row=row_ind+col_ptr[col];
			const double* val=SA->GetVal()+col_ptr[col];
			while (row<row_ind+col_ptr[col+1] && *row<col) {
				if (sparsity_pattern) add_sparsity_pattern(*row, col, 2**val);
				++row; ++val;
			}
			double coeff=0.;
			if (row<row_ind+col_ptr[col+1] && *row==col) { coeff=*val; ++row; ++val; }
			add_quadratic(col, 0., coeff);
			add_nonlinear(col);
		}
	} else {
		double val;
		for (int row=0; row<A.dim(); ++row) {
			bool is_quad=false;
			if (sparsity_pattern) {
				for (int col=0; col<row; ++col) {
					val=(*EA)(col, row);
					if (val) {
						add_sparsity_pattern(row, col, 2*val);
						is_quad=true;
					}
				}
				val=(*EA)(row, row);
				if (val) is_quad=true;
				if (!is_quad) // check rest of column
					for (int col=row+1; col<A.dim(); ++col)
						if ((*EA)(col, row)) { is_quad=true; break; }
				if (is_quad) add_quadratic(row, 0., 2*val);
				add_nonlinear(row);
			} else {
				for (int col=0; col<row; ++col)
					if ((*EA)(col, row)) { is_quad=true; break; }
				val=(*EA)(row, row);
				if (val) is_quad=true;
				if (!is_quad)
					for (int col=row+1; col<A.dim(); ++col)
						if ((*EA)(col, row)) { is_quad=true; break; }
				if (is_quad) {
					add_quadratic(row, 0., 2*val);
					add_nonlinear(row);
				}
			}
		}
	}
	if (delete_SA) delete SA;
}

void SparsityInfo::add(const UserMatrix& A, const ivector& block) {
	SparsityInfo si(2);
	si.add(A);
	add(si, block);
}

void SparsityInfo::add(const SparsityInfo& si) {
	if ((!si.linear) || (!si.nonlinear)) return; // not enough information to generate something useful

	for (map<int, LinearVariable>::iterator it(si.linear->begin()); it!=si.linear->end(); ++it) {
		pair<map<int, LinearVariable>::iterator, bool> ret(linear->insert(*it));
		if (!ret.second) ret.first->second.coeff+=it->second.coeff;
	}

	nonlinear->insert(si.nonlinear->begin(), si.nonlinear->end());

	// construct sparsity pattern
	if (sparsity_pattern && si.sparsity_pattern) {
		for (map<pair<int, int>, NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it) {
			pair<map<pair<int, int>, NonlinearConnection>::iterator, bool> ret(sparsity_pattern->insert(*it));
			if (!ret.second) ret.first->second.coeff+=it->second.coeff;
		}
	}

	if (quadratic && si.quadratic) {
		// nonquadratic variable in si
		for (map<int, QuadraticVariable>::iterator it(quadratic->begin()); it!=quadratic->end();)
			if (si.nonlinear->count(it->first) && (!si.quadratic->count(it->first))) {
				map<int, QuadraticVariable>::iterator next(it); ++next;
				quadratic->erase(it);
				it=next;
			}	else ++it;
		for (map<int, QuadraticVariable>::iterator it(si.quadratic->begin()); it!=si.quadratic->end(); ++it) {
			pair<map<int, QuadraticVariable>::iterator, bool> ret(quadratic->insert(*it));
			if (!ret.second) {
				ret.first->second.coeff_lin+=it->second.coeff_lin;
				ret.first->second.coeff_quad+=it->second.coeff_quad;
			} else if (nonlinear->count(it->first)) quadratic->erase(ret.first); // nonquadratic variable
		}
		for (map<int, QuadraticVariable>::iterator it(quadratic->begin()); it!=quadratic->end(); ++it) {
			map<int, LinearVariable>::iterator it_lin(linear->find(it->first));
			if (it_lin!=linear->end()) {
				it->second.coeff_lin+=it_lin->second.coeff;
				linear->erase(it_lin);
			}
		}
	}

	for (map<int, NonlinearVariable>::iterator it(nonlinear->begin()); it!=nonlinear->end(); ++it)
		linear->erase(it->first);
}

void SparsityInfo::add(const SparsityInfo& si, const ivector& block) {
	SparsityInfo sitmp;
	if (si.linear) {
		sitmp.linear=new map<int, SparsityInfo::LinearVariable>;
		for (map<int, SparsityInfo::LinearVariable>::iterator it(si.linear->begin()); it!=si.linear->end(); ++it)
			sitmp.linear->insert(pair<int,SparsityInfo::LinearVariable>(block[it->first], it->second));
	}

	if (si.nonlinear) {
		sitmp.nonlinear=new map<int, SparsityInfo::NonlinearVariable>;
		for (map<int, SparsityInfo::NonlinearVariable>::iterator it(si.nonlinear->begin()); it!=si.nonlinear->end(); ++it)
			sitmp.nonlinear->insert(pair<int,SparsityInfo::NonlinearVariable>(block[it->first], it->second));
	}

	if (si.quadratic) {
		sitmp.quadratic=new map<int, SparsityInfo::QuadraticVariable>;
		for (map<int, SparsityInfo::QuadraticVariable>::iterator it(si.quadratic->begin()); it!=si.quadratic->end(); ++it)
			sitmp.quadratic->insert(pair<int,SparsityInfo::QuadraticVariable>(block[it->first], it->second));
	}

	if (si.sparsity_pattern) {
		sitmp.sparsity_pattern=new map<pair<int,int>, SparsityInfo::NonlinearConnection>;
		for (map<pair<int,int>, SparsityInfo::NonlinearConnection>::iterator it(si.sparsity_pattern->begin()); it!=si.sparsity_pattern->end(); ++it)
			sitmp.sparsity_pattern->insert(pair<pair<int,int>,SparsityInfo::NonlinearConnection>(pair<int,int>(block[it->first.first], block[it->first.second]), it->second));
	}

	add(sitmp);
}

// ---------------------------------- Func -----------------------------------------------

double Func::min_eig_hess(const UserVector<double>& x, Param* param) const {
  HessMatrix H(*this, x);
  dvector eigvec(dim()); double eigval;
  H.eig(eigvec, eigval, param);  // and what are we going to do with the return code?

  return eigval;
}

Func::CurvatureType Func::add_curvatures(double a1, CurvatureType ct1, double a2, CurvatureType ct2) {
	if (!a1) return mult_curvature(a2, ct2);
	if (!a2) return mult_curvature(a1, ct1);

	if (ct1==Func::UNKNOWN || ct2==Func::UNKNOWN) return Func::UNKNOWN; // one type unknown
	if (ct1==Func::INDEFINITE || ct2==Func::INDEFINITE) return Func::INDEFINITE;

	ct1=mult_curvature(a1, ct1);
	ct2=mult_curvature(a2, ct2);

	if ((ct1==Func::LINEAR) && (ct2==Func::LINEAR)) return Func::LINEAR;
	if ((ct1&Func::CONVEX) && (ct2&Func::CONVEX)) return Func::CONVEX;
	if ((ct1&Func::CONCAVE) && (ct2&Func::CONCAVE)) return Func::CONCAVE;

	return Func::INDEFINITE;
}

Func::CurvatureType Func::mult_curvature(double a, CurvatureType ct) {
	if (a>0) return ct;
	if ((!a) || ct==Func::LINEAR) return Func::LINEAR;
	if (ct&Func::UNKNOWN) return ct; // unknown or indefinite
	return ct==Func::CONVEX ? Func::CONCAVE : Func::CONVEX;
}

ostream& operator<<(ostream& out, const Func::CurvatureType& ct) {
	switch (ct) {
		case Func::CONVEX : out << "convex"; break;
		case Func::CONCAVE : out << "concave"; break;
		case Func::LINEAR : out << "linear"; break;
		case Func::UNKNOWN : out << "unknown"; break;
		case Func::INDEFINITE : out << "indefinite"; break;
	}
	return out;	
}

// -------------------------------------- VariableIterator ------------------------------------------

VariableIterator::VariableIterator(const Func& f_, bool linear_, bool nonlinear_, bool quadratic_)
: VariableIterator_Type(linear_, nonlinear_, quadratic_), sparsity(f_.get_sparsity()), whereami(LINEAR)
{ init();
}

VariableIterator::VariableIterator(const SparsityInfo& sparsity_, bool linear_, bool nonlinear_, bool quadratic_)
: VariableIterator_Type(linear_, nonlinear_, quadratic_), sparsity(sparsity_), whereami(LINEAR)
{ init();
}

void VariableIterator::init() {
	if (linear) {
		assert(sparsity.linear);
		it_linear=sparsity.linear->begin();
		if (it_linear==sparsity.linear->end()) whereami=nonlinear ? NONLINEAR : (quadratic ? QUADRATIC : END);
	} else {
		if (nonlinear) whereami=NONLINEAR;
		else if (quadratic) whereami=QUADRATIC;
		else whereami=END;
	}
	if (nonlinear) {
		assert(sparsity.nonlinear);
		it_nonlinear=sparsity.nonlinear->begin();
		if (whereami==NONLINEAR && it_nonlinear==sparsity.nonlinear->end()) whereami=END;
	} else if (quadratic) {
		assert(sparsity.quadratic);
		it_quadratic=sparsity.quadratic->begin();
		if (whereami==QUADRATIC && it_quadratic==sparsity.quadratic->end()) whereami=END;
	}
}

int VariableIterator::operator()() const {
	assert(whereami!=END);
	switch(whereami) {
		case LINEAR: return it_linear->first;
		case NONLINEAR: return it_nonlinear->first;
		case QUADRATIC: return it_quadratic->first;
	}
	return -1;
}

double VariableIterator::coeff_lin() const {
	assert(whereami==LINEAR || whereami==QUADRATIC);
	if (whereami==LINEAR) return it_linear->second.coeff;
	return it_quadratic->second.coeff_lin;
}

double VariableIterator::coeff_quad() const {
	assert(whereami==QUADRATIC);
	return it_quadratic->second.coeff_quad;
}

void VariableIterator::operator++() {
	switch (whereami) {
		case LINEAR:
			if (++it_linear==sparsity.linear->end())
				whereami = (nonlinear && it_nonlinear!=sparsity.nonlinear->end()) ? NONLINEAR :
					((quadratic && it_quadratic!=sparsity.quadratic->end()) ? QUADRATIC :
					 END);
			break;
		case NONLINEAR:
			if (++it_nonlinear==sparsity.nonlinear->end()) whereami=END;
			break;
		case QUADRATIC:
			if (++it_quadratic==sparsity.quadratic->end()) whereami=END;
			break;
	}
}

//----------------------------------------

bool MinusFunc::compute_sparsity_pattern(const vector<dvector>& sample_set) {
	if (sparsity) return sparsity->compute_sparsity_pattern(*this, sample_set);
	bool ret=f->compute_sparsity_pattern(sample_set);
	sparsity=new SparsityInfo(f->get_sparsity());
	for (map<int, SparsityInfo::LinearVariable>::iterator it(sparsity->linear->begin()); it!=sparsity->linear->end(); ++it)
		it->second.coeff*=-1.;
	for (map<int, SparsityInfo::QuadraticVariable>::iterator it(sparsity->quadratic->begin()); it!=sparsity->quadratic->end(); ++it) {
		it->second.coeff_lin*=-1.;
		it->second.coeff_quad*=-1.;
	}
	for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(sparsity->sparsity_pattern->begin()); it!=sparsity->sparsity_pattern->end(); ++it)
		it->second.coeff*=-1.;

	return ret;
}

//----------------------------------------

bool SumFunc::compute_sparsity_pattern(const vector<dvector>& sample_set) {
	if (sparsity) return sparsity->compute_sparsity_pattern(*this, sample_set);
	assert(!g); // f and g defined, but no sparsity set
	bool ret=f->compute_sparsity_pattern(sample_set);
	if (a==1.) return ret;
	sparsity=new SparsityInfo(f->get_sparsity());
	for (map<int, SparsityInfo::LinearVariable>::iterator it(sparsity->linear->begin()); it!=sparsity->linear->end(); ++it)
		it->second.coeff*=a;
	for (map<int, SparsityInfo::QuadraticVariable>::iterator it(sparsity->quadratic->begin()); it!=sparsity->quadratic->end(); ++it) {
		it->second.coeff_lin*=a;
		it->second.coeff_quad*=a;
	}
	for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(sparsity->sparsity_pattern->begin()); it!=sparsity->sparsity_pattern->end(); ++it)
		it->second.coeff*=a;
	return ret;
}

void SumFunc::grad(UserVector<double>& y, const UserVector<double>& x) const {
  if (f) {
  	f->grad(y,x);
  	y*=a;
  }	else {
		g->grad(y,x);
	 	y*=b;
  	return;
  }
	if (g) {
	 	Pointer<UserVector<double> > z(y.getemptycopy());
		g->grad(*z, x);
		y.AddMult(b,*z);
  }
}

int SumFunc::valgrad(double& val, UserVector<double>& y, const UserVector<double>& x) const {
  int ret=0;
	if (f) {
		ret=f->valgrad(val, y, x);
		val*=a;
		y*=a;
	} else {
		ret+=g->valgrad(val, y, x);
		val*=b;
		y*=b;
		return ret;
	}
	if (g) {
		Pointer<UserVector<double> > y0(y.getemptycopy());
		double val0;
		ret+=g->valgrad(val0, *y0, x);
		val+=b*val0;
		y.AddMult(b,*y0);
	}
	return ret;
}

void SumFunc::HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
	if (f) {
		f->HessMult(y,x,z);
		y*=a;
	} else {
		g->HessMult(y,x,z);
		y*=b;
		return;
	}
	if (g) {
		Pointer<UserVector<double> > y0(y.getemptycopy());
		g->HessMult(*y0, x, z);
		y.AddMult(b,*y0);
	}
}

#ifdef FILIB_AVAILABLE
void SumFunc::grad(IntervalVector& y, const IntervalVector& x) const {
  if (f) {
  	f->grad(y,x);
  	y*=a;
  }	else {
		g->grad(y,x);
	 	y*=b;
  	return;
  }
	if (g) {
	 	IntervalVector z(y.dim());
		g->grad(z, x);
		y.AddMult(b, z);
  }
}

int SumFunc::valgrad(interval<double>& val, IntervalVector& y, const IntervalVector& x) const {
  int ret=0;
	if (f) {
		ret=f->valgrad(val, y, x);
		val*=a;
		y*=a;
	} else {
		ret=g->valgrad(val, y, x);
		val*=b;
		y*=b;
		return ret;
	}
	if (g) {
		IntervalVector y0(y.dim());
		interval<double> val0;
		ret+=g->valgrad(val0, y0, x);
		val+=b*val0;
		y.AddMult(b,y0);
	}
	return ret;
}
#endif

// -------------------------------- SepQcFunc::VariableIterator ----------------------------------

SepQcFunc::VariableIterator::VariableIterator(const SepQcFunc& f_, bool linear_, bool nonlinear_, bool quadratic_)
: VariableIterator_Type(linear_, nonlinear_, quadratic_), f(f_), whereami(END), blocknr(0)
{ find_start(); }

void SepQcFunc::VariableIterator::find_start() {
	while(blocknr<f.block.size()) {
		const SparsityInfo& si(f.get_sparsity(blocknr));
		if (linear) {
			assert(si.linear);
			it_linear=si.linear->begin();
			if (it_linear==si.linear->end()) whereami=nonlinear ? NONLINEAR : (quadratic ? QUADRATIC : END);
			else whereami=LINEAR;
		} else {
			if (nonlinear) whereami=NONLINEAR;
			else if (quadratic) whereami=QUADRATIC;
			else whereami=END;
		}
		if (nonlinear) {
			assert(si.nonlinear);
			it_nonlinear=si.nonlinear->begin();
			if (whereami==NONLINEAR && it_nonlinear==si.nonlinear->end()) whereami=quadratic ? QUADRATIC : END;
		}
		if (quadratic && (!nonlinear)) {
			assert(si.quadratic);
			it_quadratic=si.quadratic->begin();
			if (whereami==QUADRATIC && it_quadratic==si.quadratic->end()) whereami=END;
		}
		if (whereami!=END) break;
		++blocknr;
	}
}

int SepQcFunc::VariableIterator::bindex() const {
	assert(whereami!=END);
	switch(whereami) {
		case LINEAR: return it_linear->first;
		case NONLINEAR: return it_nonlinear->first;
		case QUADRATIC: return it_quadratic->first;
	}
	return -1;
}

double SepQcFunc::VariableIterator::coeff_lin() const {
	assert(whereami==LINEAR || whereami==QUADRATIC);
	if (whereami==LINEAR) return it_linear->second.coeff;
	return it_quadratic->second.coeff_lin;
}

double SepQcFunc::VariableIterator::coeff_quad() const {
	assert(whereami==QUADRATIC);
	return it_quadratic->second.coeff_quad;
}

void SepQcFunc::VariableIterator::operator++() {
	const SparsityInfo& si(f.get_sparsity(blocknr));
	switch (whereami) {
		case LINEAR:
			if (++it_linear==si.linear->end())
				whereami = (nonlinear && it_nonlinear!=si.nonlinear->end()) ? NONLINEAR :
					((quadratic && it_quadratic!=si.quadratic->end()) ? QUADRATIC :
				 	END);
			break;
		case NONLINEAR:
			if (++it_nonlinear==si.nonlinear->end())
				whereami = (quadratic && it_quadratic!=si.quadratic->end()) ? QUADRATIC : END;
			break;
		case QUADRATIC:
			if (++it_quadratic==si.quadratic->end()) whereami=END;
			break;
	}
	if (whereami==END) {
		++blocknr;
		find_start();
	}
}

// -------------------------------- SepQcFunc ----------------------------------

ostream& operator<<(ostream& out, const SepQcFunc::ftype& ft) {
	switch (ft) {
		case SepQcFunc::CONSTANT : out << "constant"; break;
		case SepQcFunc::LINEAR : out << "linear"; break;
		case SepQcFunc::QUADRATIC : out << "quadratic"; break;
		case SepQcFunc::NONQUAD : out << "nonquadratic";
	}
	return out;
}

bool SepQcFunc::compute_sparsity(int block_nr, const vector<dvector>& sample_set, bool replace_if_quadratic) {
	bool sign_changed=false;
	Pointer<SparsityInfo>& si(sparsity_block[block_nr]);

	if (s[block_nr]) {
		sign_changed=s[block_nr]->compute_sparsity_pattern(sample_set);
		si=new SparsityInfo(((const Func*)s[block_nr])->get_sparsity());
	} else {
		si=new SparsityInfo();
		si->linear=new map<int, SparsityInfo::LinearVariable>;
		si->nonlinear=new map<int, SparsityInfo::NonlinearVariable>;
		si->quadratic=new map<int, SparsityInfo::QuadraticVariable>;
		si->sparsity_pattern=new map<pair<int, int>, SparsityInfo::NonlinearConnection>;
	}

	if (A[block_nr]) si->add(*A[block_nr]);

	if (b[block_nr]) si->add(*b[block_nr]);

	if (!s[block_nr]) return sign_changed;

	const SparsityInfo& ssi(((const Func*)s[block_nr])->get_sparsity());
	if (replace_if_quadratic && ssi.quadratic->size()==ssi.nonlinear->size()) {
		Pointer<UserVector<double> > sb(new SparseVector<double>(block[block_nr].size()));
		Pointer<SparseMatrix2> sA(new SparseMatrix2(block[block_nr].size()));
		if (!ssi.quadratic->empty()) {
			for (map<pair<int,int>, SparsityInfo::NonlinearConnection>::iterator it(ssi.sparsity_pattern->begin()); it!=ssi.sparsity_pattern->end(); ++it) {
				sA->AddElement(it->first.first, it->first.second, .5*it->second.coeff);
				sA->AddElement(it->first.second, it->first.first, .5*it->second.coeff);
			}
			for (map<int, SparsityInfo::QuadraticVariable>::iterator it(ssi.quadratic->begin()); it!=ssi.quadratic->end(); ++it) {
				sA->AddElement(it->first, it->first, it->second.coeff_quad);
				sb->SetElement(it->first, .5*it->second.coeff_lin);
			}
			sA->finish();
			if (A[block_nr]) A[block_nr]=new SumMatrix((Pointer<const UserMatrix>)A[block_nr], (Pointer<const UserMatrix>)sA);
			else A[block_nr]=sA;
		} else sA=NULL;

		for (map<int, SparsityInfo::LinearVariable>::iterator it(ssi.linear->begin()); it!=ssi.linear->end(); ++it)
			sb->SetElement(it->first, .5*it->second.coeff);
		if (*sb==0.) sb=NULL;
		else if (b[block_nr]) {
			b[block_nr]=b[block_nr]->getcopy();
			*b[block_nr]+=*sb;
		} else b[block_nr]=sb;

		double constant=s[block_nr]->eval(sample_set[0]);
		if (sA) constant-=sA->xAx(sample_set[0]);
		if (sb) constant-=2*(*sb*sample_set[0]);			
		assert(finite(constant));
		c+=constant;

		s[block_nr]=NULL;
	}

  return sign_changed;
}

void SepQcFunc::set_sparsity() {
	sparsity=new SparsityInfo(2);

	for (int k=0; k<block.size(); ++k) {
		if (sparsity_block[k]) {
			sparsity->add(*sparsity_block[k], block[k]);
			continue;
		}
		if (s[k]) {
			assert(s[k]->sparsity_available());
			sparsity->add(((const Func*)(Func*)s[k])->get_sparsity(), block[k]);
		}
		if (A[k]) sparsity->add(*A[k], block[k]);
		if (b[k]) sparsity->add(*b[k], block[k]);
	}
}

SepQcFunc::SepQcFunc(const SepQcFunc &f, bool minus_func)
: SepFunc(f), A(f.block.size()), b(f.block.size()), s(f.block.size()), c(minus_func ? -f.c : f.c),
  curv_type(f.curv_type), sparsity_block(f.sparsity_block)
{ for (int i=block.size()-1; i>=0; --i) {
		if (f.A[i]) A[i]=(minus_func ? Pointer<UserMatrix>(new MinusMatrix(f.A[i])) : f.A[i]);
		if (f.b[i])
			if (minus_func) {
				b[i]=f.b[i]->getcopy();
				*b[i]*=-1;
			} else
				b[i]=f.b[i];
		if (f.s[i])
			if (minus_func) s[i]=new MinusFunc(f.s[i], f.s[i]->out_func_p, f.s[i]->out_func_log_p);
			else s[i]=f.s[i];
		if (minus_func) {
			if (curv_type[i]==Func::CONVEX) curv_type[i]=Func::CONCAVE;
			else if (curv_type[i]==Func::CONCAVE) curv_type[i]=Func::CONVEX;
			if (sparsity_block[i]) {
				sparsity_block[i]=new SparsityInfo(*sparsity_block[i]);
				if (sparsity_block[i]->linear)
					for (map<int, SparsityInfo::LinearVariable>::iterator it(sparsity_block[i]->linear->begin()); it!=sparsity_block[i]->linear->end(); ++it)
						it->second.coeff*=-1.;
				if (sparsity_block[i]->quadratic)
					for (map<int, SparsityInfo::QuadraticVariable>::iterator it(sparsity_block[i]->quadratic->begin()); it!=sparsity_block[i]->quadratic->end(); ++it) {
						it->second.coeff_lin*=-1.;
						it->second.coeff_quad*=-1.;
					}
				if (sparsity_block[i]->sparsity_pattern)
					for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(sparsity_block[i]->sparsity_pattern->begin()); it!=sparsity_block[i]->sparsity_pattern->end(); ++it)
						it->second.coeff*=-1.;
			}
		}
	}
	sparsity=f.sparsity;
	if (minus_func && sparsity) {
		sparsity=new SparsityInfo(*sparsity);
		if (sparsity->linear)
			for (map<int, SparsityInfo::LinearVariable>::iterator it(sparsity->linear->begin()); it!=sparsity->linear->end(); ++it)
				it->second.coeff*=-1.;
		if (sparsity->quadratic)
			for (map<int, SparsityInfo::QuadraticVariable>::iterator it(sparsity->quadratic->begin()); it!=sparsity->quadratic->end(); ++it) {
				it->second.coeff_lin*=-1.;
				it->second.coeff_quad*=-1.;
			}
		if (sparsity->sparsity_pattern)
			for (map<pair<int, int>, SparsityInfo::NonlinearConnection>::iterator it(sparsity->sparsity_pattern->begin()); it!=sparsity->sparsity_pattern->end(); ++it)
				it->second.coeff*=-1.;
	}
}

SepQcFunc::SepQcFunc(const SepQcFunc& f, const UserVector<double>& point, const int degree)
: SepFunc(f), c(f.c), curv_type(f.curv_type), sparsity_block(f.block.size())
{	switch (degree) {
		case 2: {
			s.resize(block.size());
			A=f.A;
			b=f.b;
			for (int k=0; k<block.size(); k++) {
				if (!f.s[k]) continue;
		    Pointer<const UserVector<double> > x0(point.getcopy(block[k]));
		    Pointer<UserVector<double> > grad(point.getemptycopy(x0->size()));
			  double val;
	    	f.s[k]->valgrad(val, *grad, *x0);

				c+=val - *grad**x0;
				if (b[k]) { b[k]=b[k]->getcopy(); b[k]->AddMult(2., *grad); }
				else { *grad*=.5; b[k]=grad; }

	    	Pointer<const UserMatrix> H(new HessMatrix(Pointer<const Func>(f.s[k]), x0));
				c+=.5*H->xAx(*x0);
				b[k]->AddMult(.5,*H * *x0);

				A[k]=new SumMatrix(H, A[k] ? Pointer<const UserMatrix>(A[k]) : Pointer<const UserMatrix>(NULL), .5, 1.);
			}
		} break;

		case 1: {
			s.resize(block.size());
			A.resize(block.size());
			b=f.b;
	  	for (int k=0; k<block.size(); k++) {
			  curv_type[k]=Func::LINEAR;
				if ((!f.s[k]) && (!f.A[k])) continue;
				Pointer<UserVector<double> > x0(point.getcopy(block[k]));
				Pointer<UserVector<double> > grad(point.getemptycopy(x0->size()));
				double val=0;
				if (f.s[k]) f.s[k]->valgrad(val, *grad, *x0);
				if (f.A[k]) {
					val+=f.A[k]->xAx(*x0);
					grad->AddMult(2.,(*f.A[k])* *x0);
				}
	 			c+=val - *grad * *x0;

				if (b[k]) { b[k]=b[k]->getcopy(); b[k]->AddMult(.5, *grad); }
				else { *grad*=.5;	b[k]=grad; }
			}
		} break;

		case 0: {
			s.resize(block.size());
			A.resize(block.size());
			b.resize(block.size());
			for (int k=0; k<curv_type.size(); ++k) curv_type[k]=Func::LINEAR;
			c=f.eval(point);
		} break;

		default:
			out_err << "SepQcFunc: Taylor approximation of degree " << degree << " not supported. Aborting." << endl;
			exit(-1);
	}
}

void SepQcFunc::add_var(int index, int bnum) {
	assert(bnum>=block.size() || (!A[bnum] && !s[bnum]));
	
	if (bnum>=block.size()) { // new block
		block.resize(bnum+1);
		A.resize(bnum+1);
		b.resize(bnum+1);
		s.resize(bnum+1);
		block[bnum].resize(1);
		block[bnum][0]=index;
		sparsity_block.resize(bnum+1);
		curv_type.resize(bnum+1);
		curv_type[bnum]=Func::UNKNOWN;
	} else { // old block
	  block[bnum].resize(block[bnum].size()+1);
  	block[bnum][block[bnum].size()-1]=index;
		if (b[bnum]) {
			b[bnum]=b[bnum]->getcopy();
			b[bnum]->resize(b[bnum]->size()+1);
		}
	}
	
	dim_++;
}

void SepQcFunc::addmult(const SepQcFunc& f, double a_, double b_) {
	assert(block.size()==f.block.size());
	for (int k=0; k<block.size(); k++) {
		assert(block[k]==f.block[k]);
		if (A[k] || f.A[k]) A[k]=new SumMatrix((UserMatrix*)A[k], (UserMatrix*)f.A[k], a_, b_);
		if (b[k]) {
			if (f.b[k] || a_!=1.) {
				b[k]=b[k]->getcopy();
				if (a_!=1.) *b[k]*=a_;
				if (f.b[k]) {
					if (b_!=1.) b[k]->AddMult(b_, *f.b[k]);
					else *b[k]+=*f.b[k];
				}
			}
		} else if (f.b[k]) { // no b[k]
		 	b[k]=f.b[k]->getcopy();
		 	if (b_!=1.) *b[k]*=b_;
		}
		if (s[k] || f.s[k]) s[k]=new SumFunc((Func*)s[k], (Func*)f.s[k], a_, b_);
	}

}

double SepQcFunc::eval(const UserVector<double>& x) const {
	double val(c);
	for (int k=block.size()-1; k>=0; --k) {
		if ((!A[k]) && (!b[k]) && (!s[k])) continue;
		Pointer<UserVector<double> > xb(x.getcopy(block[k]));
		if (A[k]) val+=A[k]->xAx(*xb);
		if (b[k]) val+=2*(*b[k]**xb);
		if (s[k]) val+=s[k]->eval(*xb);
	}
	return val;
}

void SepQcFunc::grad(UserVector<double>& g, const UserVector<double>& x) const {
	g=0;
  for (int k=block.size()-1; k>=0; --k) {
		if ((!A[k]) && (!b[k]) && (!s[k])) continue;
    Pointer<UserVector<double> > xb(x.getcopy(block[k]));
    Pointer<UserVector<double> > gb(g.getemptycopy(block[k].size()));
    if (s[k]) {
			s[k]->grad(*gb, *xb);
			if (A[k]) A[k]->AddMult(*gb, *xb, 2.);
	    if (b[k]) gb->AddMult(2., *b[k]);
		} else {
			if (A[k]) {
				A[k]->MultV(*gb, *xb);
				if (b[k]) *gb+=*b[k];
				*gb*=2.;
			} else
				if (b[k]) {
					*gb=*b[k];
					*gb*=2.;
				}
		}
    g.set_block(*gb, block[k]);
  }
}

void SepQcFunc::grad(UserVector<double>& g, const UserVector<double>& x, int bnum) const {
	if (s[bnum]) {
		s[bnum]->grad(g, x);
		if (A[bnum]) A[bnum]->AddMult(g, x, 2.);
		if (b[bnum]) g.AddMult(2., *b[bnum]);
	} else {
	  if (A[bnum]) A[bnum]->MultV(g,x);
		else g=0;
  	if (b[bnum]) g+=*b[bnum];
  	g*=2;
	}
}

int SepQcFunc::valgrad(double& val, UserVector<double>& g, const UserVector<double>& x) const {
	g=0;
	val=c;
	int ret=0;
	for (int k=block.size()-1; k>=0; --k) {
		if ((!A[k]) && (!b[k]) && (!s[k])) continue;
		Pointer<UserVector<double> > xb(x.getcopy(block[k]));
		Pointer<UserVector<double> > gb(g.getemptycopy(block[k].size()));
		if (s[k]) {
			double valb;
			ret+=s[k]->valgrad(valb, *gb, *xb);
			val+=valb;
			if (A[k]) {
				A[k]->AddMult(*gb, *xb, 2.);
				val+=A[k]->xAx(*xb);
			}
			if (b[k]) {
				gb->AddMult(2., *b[k]);
				val+=2*(*b[k]**xb);
			}
		} else {
			if (A[k]) {
				A[k]->MultV(*gb, *xb);
				val+=*xb**gb;
				if (b[k]) {
					*gb+=*b[k];
					val+=(*b[k])**xb;
				}
				*gb*=2.;
			} else
				if (b[k]) {
					*gb=*b[k];
					*gb*=2.;
					val+=*gb**xb;
				}
		}
	}
	return ret;
}

void SepQcFunc::HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z) const {
	y=0;
  for (int k=block.size()-1; k>=0; --k) {
		if ((!A[k]) && (!s[k])) continue;
    Pointer<UserVector<double> > yb(y.getemptycopy(block[k].size()));
    Pointer<UserVector<double> > xb(x.getcopy(block[k]));
    Pointer<UserVector<double> > zb(z.getcopy(block[k]));
    if (s[k]) s[k]->HessMult(*yb, *xb, *zb);
    if (A[k]) A[k]->AddMult(*yb, *zb, 2.); // yb += 2 A[i] * zb;
    y.set_block(*yb, block[k]);
  }
}

void SepQcFunc::HessMult(UserVector<double>& y, const UserVector<double>& x, const UserVector<double>& z, int k) const {
	y=0;
  if (s[k]) s[k]->HessMult(y, x, z);
  if (A[k]) A[k]->AddMult(y, z, 2.);
}

#ifdef FILIB_AVAILABLE
bool SepQcFunc::is_interval_compliant() const {
	bool itis=true;
	for (int k=0; k<block.size() && itis; k++)
		itis&=(!s[k]) || s[k]->is_interval_compliant();
	return itis;
}

interval<double> SepQcFunc::eval(const IntervalVector& x) const {
	interval<double> val(c);
	for (int k=0; k<block.size(); k++)
		if (b[k] || A[k] || s[k]) val+=eval(IntervalVector(x, block[k]), k);
	return val;
}

interval<double> SepQcFunc::eval(const IntervalVector& x, int k) const {
	interval<double> val;
	if (b[k]) val=2.*(x*(*b[k]));
	if (A[k]) {
		assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A[k]));
		val+=((const IntervalCompliantMatrix*)(const UserMatrix*)A[k])->xAx(x);
	}
	if (s[k]) val+=s[k]->eval(x);
	return val;
}

void SepQcFunc::grad(IntervalVector& g, const IntervalVector& x) const {
	g=interval<double>(0.);
	for (int k=0; k<block.size(); k++)
		if (b[k] || A[k] || s[k]) {
			IntervalVector gk(block[k].size());
			grad(gk, IntervalVector(x, block[k]), k);
			g.set_block(gk, block[k]);
		}
}

void SepQcFunc::grad(IntervalVector& g, const IntervalVector& x, int k) const {
	if (A[k]) {
		assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A[k]));
		((const IntervalCompliantMatrix*)(const UserMatrix*)A[k])->MultV(g, x);
		g*=2.;
		if (s[k]) {
			IntervalVector g0(g.dim());
			s[k]->grad(g0, x);
			g+=g0;
		}
	} else if (s[k]) s[k]->grad(g, x);
	else g=interval<double>(0.);
	if (b[k]) g.AddMult(2., *b[k]);
}

int SepQcFunc::valgrad(interval<double>& val, IntervalVector& g, const IntervalVector& x) const {
	val=c;
	g=interval<double>(0.);
	int ret=0;
	for (int k=0; k<block.size(); k++) {
		if ((!A[k]) && (!b[k]) && (!s[k])) continue;
		IntervalVector xk(x, block[k]);
		IntervalVector gk(xk.dim());
		if (A[k]) {
			assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A[k]));
			((const IntervalCompliantMatrix*)(const UserMatrix*)A[k])->MultV(gk, x);
			val+=xk*gk;
			gk*=2.;
			if (s[k]) {
				interval<double> valk;
				IntervalVector gk2(xk.dim());
				ret+=s[k]->valgrad(valk, gk2, xk);
				gk+=gk2;
				val+=valk;
			}
		} else if (s[k]) {
			interval<double> valk;
			ret+=s[k]->valgrad(valk, gk, xk);
			val+=valk;
		}
		if (b[k]) {
			val+=2.*(xk*(*b[k]));
			gk.AddMult(2, *b[k]);
		}
		g.set_block(gk, block[k]);
	}
	return ret;
}

#endif

Func::CurvatureType SepQcFunc::get_curvature() const {
	Func::CurvatureType ct(Func::LINEAR);
	for (int k=0; k<block.size(); ++k) ct=Func::add_curvatures(1., ct, 1., curv_type[k]);
  return ct;
}

void SepQcFunc::set_curvature(Func::CurvatureType ct) {
	if (block.size()==1) curv_type[0]=ct;
	else {
		out_err << "SepQcFunc::set_curvature: Setting one curvature for more than one block not supported. You need to specify a block number! Aborting." << endl;
		exit(-1);
	}
}

void SepQcFunc::print(ostream &out) const {
  out << "SepQcfunc: dim=" << dim() << " c= " << c << endl;
  for (int i=0; i<block.size(); i++) {
		if ((!A[i]) && (!b[i]) && (!s[i])) continue;
  	out << "block " << i << ": " << block[i];
		if (A[i]) out << "A[" << i << "]: " << *A[i];
		if (b[i]) out << "b[" << i << "]: " << *b[i];
		if (s[i]) out << "s[" << i << "]: " << *s[i];
	}
}

void SepQcFunc::print(ostream& out, vector<Pointer<char> > var_names) const {
	out << c;
	ios::fmtflags old_flags=out.flags();
	for (int k=0; k<block.size(); k++) {
		if ((!A[k]) && (!b[k]) && (!s[k])) continue;
		out << " + (";
		bool first=true;
		if (b[k]) {
			for (int i=0; i<block[k].size(); i++) {
				double coeff=2* (*b[k])(i);
				if (coeff) {
					if (!first) out << " ";
					if (fabs(coeff)!=1) { out.flags(old_flags | ios::showpos); out << coeff << "*"; out.flags(old_flags); }
					if (!first && coeff==1) out << "+";
					if (coeff==-1) out << "-";
					if (var_names[block[k][i]]) out << var_names[block[k][i]];
					else out << "var" << block[k][i];
					first=false;
				}
			}
		}
		if (A[k]) {
			/*if (!first) out << " + "; */first=false;
			Pointer<ExtUserMatrix> ExtA(dynamic_cast<ExtUserMatrix*>((UserMatrix*)A[k]));
			if (!ExtA) ExtA=new DenseMatrix(*A[k]);
			double coeff;
			for (int i=0; i<ExtA->dim(); i++) {
				for (int j=0; j<i; j++) {
					coeff=(*ExtA)(i,j)+(*ExtA)(j,i);
					if (!coeff) continue;
					out.flags(old_flags | ios::showpos); out << " " << coeff; out.flags(old_flags);
					if (var_names[block[k][i]]) out << "*" << var_names[block[k][i]];
					else out << "*var" << block[k][i];
					if (var_names[block[k][j]]) out << "*" << var_names[block[k][j]];
					else out << "*var" << block[k][j];
				}
				coeff=(*ExtA)(i,i);
				if (!coeff) continue;
				out.flags(old_flags | ios::showpos); out << " " << coeff; out.flags(old_flags);
				if (var_names[block[k][i]]) out << "*sqr(" << var_names[block[k][i]] << ")";
				else out << "*var" << block[k][i] << "^2";
			}
		}
		if (s[k]) {
			if (!first) out << " + "; first=false;
			out << "s" << k << "(x)";
		}
		out << ")";
	}
	out << endl;
	for (int k=0; k<block.size(); k++) {
//		if (A[k] && !dynamic_cast<ExtUserMatrix*>((UserMatrix*)A[k])) out << "  q" << k << " = " << *A[k];
		if (s[k]) out << "  s" << k << " = " << *s[k];
	}
/*	for (int k=0; k<block.size(); k++)
		if (sparsity_block[k] && (A[k] || s[k] || b[k])) out << "  sparsity" << k << " = " << *sparsity_block[k];*/
}
