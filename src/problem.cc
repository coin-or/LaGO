// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

#include "problem.h"
#include "decomp.h"

// -------------------------- MinlpProblem ---------------------

MinlpProblem::MinlpProblem(const MinlpProblem& p, int k)
: dim_(p.block[k].size()), problem_type(p.problem_type), prob_name(p.prob_name),
  out_problem_p(p.out_problem_p), out_problem_log_p(p.out_problem_log_p),
	lower(p.lower, p.block[k]), upper(p.upper, p.block[k]), block(1), discr(p.block[k].size()),
	primal_point(p.primal_point, p.block[k])
{	block[0].resize(dim());
	for (int i=0; i<dim(); i++) {
		block[0][i]=i;
		int i0=p.block[k][i];
		discr[i]=p.discr[i0];
		if (discr[i]) i_discr.push_back(i);
		else i_cont.push_back(i);
		var_names.push_back(p.var_names[i0]);
	}
}

void MinlpProblem::add_var(int i_, int bnum, bool discr_, double lower_, double upper_, char* name) {
  for (int i=0; i<block.size(); i++)
    for (int j=0; j<block[i].size(); j++)
      if (i_==block[i][j]) {
        out_err << "MinlpProblem::add_var: Index " << i_ << " is in block " << i << ". Cannot add him again." << endl;
				exit(-1);
        return;
      }

  // Resizes block-structure, if necessary.
  if (block.size()<=bnum) block.resize(bnum+1);

  // Resizes the block and adds the index.
  block[bnum].resize(block[bnum].size()+1);
  block[bnum][block[bnum].size()-1]=i_;

  if (dim_<=i_) dim_=i_+1;

  if (var_names.size()<dim_) var_names.resize(dim_);

  var_names[i_]=(name ? strdup(name) : NULL);

  if (discr_) i_discr.push_back(i_);
  else i_cont.push_back(i_);

  if (discr.size()<dim_) discr.resize(dim_);
  discr[i_]=discr_;

  if (lower.size()<dim_) lower.resize(dim_);
  lower[i_]=lower_ < -INFINITY ? -INFINITY : lower_;

  if (upper.size()<dim_) upper.resize(dim_);
  upper[i_]=upper_ > INFINITY ? INFINITY : upper_;

  if (primal_point.size()<dim_) primal_point.resize(dim_);

	if (obj) obj->add_var(i_, bnum);
	for (int c=0; c<con.size(); c++) con[c]->add_var(i_, bnum);
};

void MinlpProblem::add_con(Pointer<SepQcFunc> f, bool eq, char* name) {
  con.push_back(f);
  con_names.push_back(name ? strdup(name) : NULL);
  con_eq.push_back(eq);
}

void MinlpProblem::del_con(int connr) {
	assert(connr<con.size());
  vector<Pointer<SepQcFunc> >::iterator it_con(con.begin());
  vector<Pointer<char> >::iterator it_con_names(con_names.begin());
  vector<bool>::iterator it_con_eq(con_eq.begin());
	for (int c=0; c<connr; c++) { it_con++; it_con_names++; it_con_eq++; }
  con.erase(it_con);
  con_names.erase(it_con_names);
  con_eq.erase(it_con_eq);
}

void MinlpProblem::taylor_approx(UserVector<double>& point, const int degree) {
  bool really_needed=false;
  for (int k=0; (!really_needed) && k<block.size(); k++) really_needed=(obj->s[k])||(degree==1 && obj->A[k]);

  if (really_needed) obj=new SepQcFunc(*obj, point, degree);  // get a degree-th order taylor approximation of this function

  for (int c=0; c<con.size(); c++) {
	  really_needed=false;
  	for (int k=0; (!really_needed) && k<block.size(); k++) really_needed=(con[c]->s[k])||(degree==1 && con[c]->A[k]);
    if (really_needed) con[c]=new SepQcFunc(*con[c], point, degree);  // get a degree-th order taylor approximation of this function
  }

  problem_type=QQP;
}

int MinlpProblem::feasible(const UserVector<double>& x, double tol, ostream* out) {
  int ret=0;
	for (int i=0; i<dim(); i++) {
		if ((x(i)!=x(i)) || lower(i)-tol>x(i)) {
			ret++;
			if (out) {
				*out << "! box   : ";
				if (var_names[i]) *out << var_names[i];
				else *out << "var" << i;
				*out << ": " << x(i) << " >= " << lower(i) << endl;
			}
		}
		if ((x(i)!=x(i)) || upper(i)+tol<x(i)) {
			ret++;
			if (out) {
				*out << "! box   : ";
				if (var_names[i]) *out << var_names[i];
				else *out << "var" << i;
				*out << ": " << x(i) << " <= " << upper(i) << endl;
			}
		}
	}

  for (int i=0; i<con.size(); i++) {
    double val=con[i]->eval(x);
    if ((con_eq[i] && ((val!=val) || (fabs(val)>tol))) || ((!con_eq[i]) && ((val!=val) || (val>tol)))) {
      ret++;
			if (!out) continue;
      *out << "! ";
      bool linear=true; for (int k=0; linear && k<block.size(); k++) linear&=(!con[i]->A[k]) && (!con[i]->s[k]);
      *out << i << (linear ? ": linear: " : ": nonlin: ");
			if (con_names[i]) *out << con_names[i];
			else *out << "con" << i;
			*out << ": " << val << (con_eq[i] ? " = 0" : " <= 0") << endl;
    }
  }

  return ret;
}

double MinlpProblem::compute_scale(char* option, int c, double eps, UserVector<double>& x) {
  if (!strcmp(option, "none")) return 1.;

  if (!strcmp(option, "condition Jac")) {
	  Pointer<UserVector<double> > g(x.getemptycopy());
	  con[c]->grad(*g, x);
	  return sqrt(g->sq_norm2())/eps;
  }

  if (!strcmp(option, "equalized")) {
    double val=con[c]->eval(x);
    if (val>rtol || (con_eq[c] && val<-rtol))
      return fabs(val)/eps;
    else
      return 1.;
      }
  if (!strcmp(option, "combo")) {
    // do some scaling in some direction, somehow, somewhere, ...
    return 1.;
  }

  if (!strcmp(option, "W*gradient")) {
    // do new idea of scaling \tilde h(x)=h(x)/\norm{W\nabla h(\hat x)}, \quad  W=\Diag(\ub x-\lb x)
    Pointer<UserVector<double> > g(x.getemptycopy());
    con[c]->grad(*g, x);
    g->diagmult(*g,upper-lower);
    return sqrt(g->sq_norm2())/eps;

  }

  if (!strcmp(option, "adjust by box ends")) {
    double val=fabs(con[c]->eval(lower));
    double val2=fabs(con[c]->eval(upper));
    return MAX(val, MAX(val2, 1.));
  }

  return atof(option);
}

int MinlpProblem::scale(UserVector<double>& x, Pointer<Param> param) {
  dvector con_scale(con.size(),1.);

  int scaled=0; // number of scaled constraints

  double eps=param->get_d("constraint scaling eps", 1.);
  char* default_option=param->get("constraint scaling default", "none");

  if (strcmp(default_option, "none"))
    for (int c=0; c<con.size(); c++) con_scale[c]=compute_scale(default_option, c, eps, x);
	delete default_option;
	
  int max_pattern=param->get_i("constraint scaling max patterns", 0);
  assert(max_pattern<1000);
  for (int i=0; i<=max_pattern; i++) {
    char* pname=new char[strlen("constraint scaling pattern ")+4];
    sprintf(pname, "constraint scaling pattern %i", i);
    char* pattern=param->get(pname);
    if (pattern) {
      pattern=strdup(pattern);
      char* option=pattern;
    	for (; *option && (*option!=' '); option++);
    	if (*option==' ') *(option++)=0;
    	else {
    	  out_err << "Skipping " << pname << ": " << pattern << endl;
    	  delete pname;
    	  delete pattern;
    	  continue;
    	}
	    for (int c=0; c<con.size(); c++)
	      if (strstr(con_names[c], pattern)) con_scale[c]=compute_scale(option, c, eps, x);

	    delete pattern;
	  }
	  delete pname;
  }

  for (int c=0; c<con.size(); c++) {
		if (con_scale(c)>1+rtol) {
		  //		 	double factor=1/MIN(con_scale(c), 1.e+6);
		   	double factor=1/con_scale(c);
//		 	out_problem_log << "Scaling " << con_names[c] << " by " << factor << endl;
		 	for (int k=0; k<con[c]->block.size(); k++) {
		 	  if (con[c]->A[k]) con[c]->A[k]=new SumMatrix(&*con[c]->A[k], NULL, factor);
		  	if (con[c]->b[k]) {
					con[c]->b[k]=con[c]->b[k]->getcopy();
					*con[c]->b[k]*=factor;
				}
			  if (con[c]->s[k]) con[c]->s[k]=new SumFunc(&*con[c]->s[k], NULL, factor);
			}
		  con[c]->c*=factor;
		  scaled++;
		}
		//		else  out_problem_log << "Not Scaling " << con_names[c] << " because " << con_scale(c) << endl;
	}

	return scaled;
}

void MinlpProblem::get_sparsity(vector<Pointer<SparsityInfo> >& si) const {
	si.clear(); si.reserve(block.size());
	for (int k=0; k<block.size(); ++k) si.push_back(get_sparsity(k));
}

Pointer<SparsityInfo> MinlpProblem::get_sparsity(int k) const {
	if (!obj->sparsity_available(k)) return NULL;
	Pointer<SparsityInfo> si(new SparsityInfo(obj->get_sparsity(k)));
	for (int c=0; c<con.size(); ++c) {
		if (!con[c]->sparsity_available(k)) return NULL;
		si->add(con[c]->get_sparsity(k));
	}
	return si;
}


void MinlpProblem::print(ostream& out) const {
  out << "This is a " << problem_type << " problem of the dimension " << dim_ << ": ";
	if (prob_name) out << prob_name << endl;
	else out << "noname" << endl;
  out << block.size() << " blocks:";
  for (int i=0; i<block.size(); i++) {
    out << endl << " block " << i << ": ";
    for (int j=0; j<block[i].size(); j++)
      out << block[i][j] << "(" << (var_names[block[i][j]] ? &*var_names[block[i][j]] : "") << ") ";
  }
  out << endl << "discrete variables: ";
  for (int i=0; i<i_discr.size(); i++)
    if (var_names[i_discr[i]]) out << var_names[i_discr[i]] << " ";
    else out << "var" << i_discr[i] << " ";
  out << endl << "continuous variables: ";
  for (int i=0; i<i_cont.size(); i++)
    if (var_names[i_cont[i]]) out << var_names[i_cont[i]] << " ";
    else out << "var" << i_cont[i] << " ";
  out << endl;
  for (int i=0; i<discr.size(); i++) {
    if (var_names[i]) out << var_names[i];
    else out << "var" << i;
    out << "\tis " << (discr[i] ? "discrete" : "continuous") << " \tlow: " << lower[i] << " \tup: " << upper[i] << " \tstart: " << primal_point[i] << endl;
  }

	if (obj) { out << "objective:\t "; obj->print(out, var_names); }
	else out << "no objective" << endl;

	out << con.size() << " constraints: " << endl;
	for (int c=0; c<con.size(); c++) {
		out << "con " << c;
		if (con_names[c]) out << "(" << con_names[c] << ")";
		out << ":\t 0 " << (con_eq[c] ? "= " : ">= "); con[c]->print(out, var_names);
	}
}

void MinlpProblem::print_as_gams(ostream& out) const {
	out << "Variables ";
	for (int i=0; i<dim(); i++) {
		if (var_names[i]) out << var_names[i] << ',' << endl;
		else out << "var" << i << ',' << endl;
	}
	for (int i=0; i<i_discr.size(); i++)
		if (var_names[i_discr[i]]) out << "my" << var_names[i_discr[i]] << ',' << endl;
		else out << "myvar" << i_discr[i] << ',' << endl;
	out << "myobjvar;" << endl;

	if (i_discr.size()) {
		out << "Binary Variables ";
		for (int i=0; i<i_discr.size(); i++) {
			if (i) out << ',';
			if (var_names[i_discr[i]]) out << "my" << var_names[i_discr[i]] << endl;
			else out << "myvar" << i_discr[i] << endl;
//			assert(lower[i]==0 && upper[i]==1);
		}
		out << ";" << endl;
	}

	out << "Equations ";
	for (int c=0; c<con.size(); c++) {
		if (con_names[c]) out << con_names[c] << ',' << endl;
		else out << "con" << c << ',' << endl;
	}
	for (int i=0; i<i_discr.size(); i++) {
		out << "bin" << i << ',' << endl;
	}
	out << "myobj;" << endl;

	out << "myobj..\t myobjvar =E= "; obj->print(out, var_names);
	out << ";" << endl;

	for (int c=0; c<con.size(); c++) {
		if (con_names[c]) out << con_names[c];
		else out << "con" << c;
		out << "..\t 0 " << (con_eq[c] ? "=E= " : " =G= "); con[c]->print(out, var_names);
		out << ";" << endl;
	}
	for (int i=0; i<i_discr.size(); i++) {
		int i0=i_discr[i];
		out << "bin" << i0 << "..\t my";
		if (var_names[i0]) out << var_names[i0] << " =E= (" << var_names[i0];
		else out << "var" << i0 << " =E= (var" << i0;
		out << " - (" << lower[i0] << ") ) / " << (upper[i0]-lower[i0]) << ";" << endl << endl;
	}

	for (int i=0; i<dim(); i++) {
		if (lower[i]>-INFINITY) {
			if (var_names[i]) out << var_names[i];
			else out << "var" << i;
			out << ".lo = " << lower[i] << ';' << endl;
		}
		if (upper[i]<INFINITY) {
			if (var_names[i]) out << var_names[i];
			else out << "var" << i;
			out << ".up = " << upper[i] << ';' << endl;
		}
		if (var_names[i]) out << var_names[i];
		else out << "var" << i;
		out << ".l = " << primal_point(i) << ';' << endl;
	}

	out << "Model m / all /; " << endl << "Solve m using MINLP minimizing myobjvar;" << endl;
}

//  -------------------------- PenaltyFunc ----------------------------------

double MinlpPenaltyFunc::eval(const UserVector<double>& x) const {
  double val=0;

  for (int i=0; i<minlp->con.size(); i++)
    if (minlp->con_eq[i]) val+=delta(i)*fabs(minlp->con[i]->eval(x));
    else {
      double conval=minlp->con[i]->eval(x);
      val+=delta(i)*MAX(0, conval);
    }

  return val;
}

int MinlpPenaltyFunc::valgrad(double& val, UserVector<double>& y, const UserVector<double>& x) const {
  y=0;

  Pointer<UserVector<double> > congrad(y.getemptycopy(minlp->dim()));
  double conval;

  for (int i=0; i<minlp->con.size(); i++) {
    minlp->con[i]->valgrad(conval, *congrad, x);
    if (minlp->con_eq[i]) {
      val+=delta(i)*fabs(conval);
      if (conval<-rtol) y-=delta(i)**congrad;
      else if (conval>rtol) y+=delta(i)**congrad;
    } else {
      val+=delta(i)*MAX(0, conval);
      if (conval>rtol) y+=delta(i)**congrad;
    }
  }
  return 0;
}

void MinlpPenaltyFunc::grad(UserVector<double>& y, const UserVector<double>& x) const {
  double val;
  valgrad(val, y, x);
}

void MinlpPenaltyFunc::HessMult(UserVector<double>& y, const UserVector<double>& z, const UserVector<double>& x) const {
  y=0;

  Pointer<UserVector<double> > conhessmult(y.getemptycopy(minlp->dim()));
  double conval;

  for (int i=0; i<minlp->con.size(); i++) {
    conval=minlp->con[i]->eval(z);
    minlp->con[i]->HessMult(*conhessmult, z, x);
    if (minlp->con_eq[i]) {
      if (conval<-rtol) y-=delta(i)**conhessmult;
      else if (conval>rtol) y+=delta(i)**conhessmult;
    }
    else if (conval>rtol) y+=delta(i)**conhessmult;
  }
}

// --------------------------------------- MipProblem ---------------------------------------------

MipProblem::MipProblem(int nr_col_, int nr_row_)
: nr_col(nr_col_), nr_row(nr_row_),
	col_lower(nr_col_, -INFINITY), col_upper(nr_col_, INFINITY), col_names(nr_col_), col_types(nr_col_, CONTINUOUS),
	obj(new SparseVector<double>(nr_col_)), obj_const(0.),
	A(nr_row_, nr_col_), row_lower(nr_row_), row_upper(nr_row_), row_names(nr_row_)
{ }

MipProblem::MipProblem(const MinlpProblem& minlp)
: nr_col(minlp.dim()), nr_row(minlp.con.size()),
  col_lower(minlp.lower), col_upper(minlp.upper), col_names(minlp.var_names), col_types(minlp.dim(), CONTINUOUS),
	A(minlp.con.size(), minlp.dim()), row_lower(minlp.con.size(), -INFINITY), row_upper(minlp.con.size()), row_names(minlp.con_names)
{ for (int i=0; i<minlp.i_discr.size(); i++) setColType(minlp.i_discr[i], BINARY);
	for (int c=0; c<minlp.con.size(); c++) {
		if (minlp.con_eq[c]) row_lower[c]=-minlp.con[c]->c;
		row_upper[c]=-minlp.con[c]->c;
		for (int k=0; k<minlp.block.size(); k++) {
			assert((!minlp.con[c]->A[k]) && (!minlp.con[c]->s[k]));
			if (minlp.con[c]->b[k])
				for (int i=0; i<minlp.block[k].size(); i++)
					A.AddElement(c, minlp.block[k][i], 2*(*minlp.con[c]->b[k])(i));
		}
	}
	for (int k=0; k<minlp.block.size(); k++) assert((!minlp.obj->A[k]) && (!minlp.obj->s[k]));
	Pointer<UserVector<double> > b(new SparseVector<double>(minlp.obj->b, minlp.block)); *b*=2;
	setObj(b, minlp.obj->c);
	finish();
}

MipProblem::MipProblem(const dvector& lower, const dvector& upper, const vector<int>& i_discr, const vector<Pointer<char> >& names)
: nr_col(lower.dim()), nr_row(0),
  col_lower(lower), col_upper(upper), col_types(lower.dim(), CONTINUOUS),
	A(0, lower.dim())
{	if (names.size()) col_names=names;
	else col_names.resize(lower.dim());
	for (vector<int>::const_iterator it(i_discr.begin()); it!=i_discr.end(); ++it)
		col_types[*it]=(col_lower[*it]==0 && col_upper[*it]==1) ? BINARY : INTEGER;
}

void MipProblem::reset() {
	nr_row=0;
	row_lower.resize(0);
	row_upper.resize(0);
	row_names.resize(0);
	A.resize(0, dim());
}

void MipProblem::setColBounds(int index, double low, double up) {
	assert(index<=dim());
	col_lower[index]=low;
	col_upper[index]=up;
}

void MipProblem::setColType(int index, VarType type) {
	assert(index<=dim());
	col_types[index]=type;
}

void MipProblem::setColName(int index, Pointer<char> name) {
	assert(index<=dim());
	col_names[index]=name;
}

void MipProblem::setCol(int index, double low, double up, VarType type, Pointer<char> name) {
	assert(index<=dim());
	col_lower[index]=low;
	col_upper[index]=up;
	col_types[index]=type;
	col_names[index]=name;
}

void MipProblem::setObj(Pointer<UserVector<double> > obj_, double obj_const_) {
	assert(obj_->dim()==dim());
	obj=obj_;
	obj_const=obj_const_;
}

void MipProblem::setRow(int index, const UserVector<double>& row, double low, double up, Pointer<char> name) {
	for (int i=0; i<dim(); i++) A.AddElement(index, i, row(i));
	row_lower[index]=low;
	row_upper[index]=up;
	row_names[index]=name;
}

void MipProblem::setRow(int index, const UserVector<double>& row, const ivector& indices, double low, double up, Pointer<char> name) {
	for (int i=0; i<indices.size(); i++) A.AddElement(index, indices(i), row(i));
	row_lower[index]=low;
	row_upper[index]=up;
	row_names[index]=name;
}

void MipProblem::setRowBounds(int index, double low, double up) {
	row_lower[index]=low;
	row_upper[index]=up;
}

void MipProblem::resize_col(int nr_col_) {
	int oldsize=dim();
	nr_col=nr_col_;

	col_lower.resize(dim());
	col_upper.resize(dim());
	col_types.resize(dim(), CONTINUOUS);
	col_names.resize(dim());
	for (int i=oldsize; i<dim(); i++) {
		col_lower[i]=-INFINITY;
		col_upper[i]=INFINITY;
	}

	if (obj) obj=NULL;

	A.resize(nr_row, nr_col);

}

void MipProblem::resize_row(int nr_row_) {
	int oldsize=rows();
	nr_row=nr_row_;

	row_lower.resize(rows());
	row_upper.resize(rows());
	row_names.resize(rows());
	A.resize(nr_row, nr_col);
}

void MipProblem::finish() {
	assert(obj);
	A.finish();
}

ostream& operator<<(ostream& out, const MipProblem& mip) {
	out << "Objective: ";
	if (!mip.obj) out << "not set" << endl;
	else {
		out << mip.obj_const << " ";
		double coeff;
		for (int i=0; i<mip.dim(); i++) {
			if (!(coeff=(*mip.obj)(i))) continue;
			out << setiosflags(ios_base::showpos) << coeff << resetiosflags(ios_base::showpos);
			if (mip.col_names[i]) out << mip.col_names[i] << " ";
			else out << "col" << i << " ";
		}
		out << endl;
	}
	out << "Columns (" << mip.dim() << "): " << endl;
	for (int i=0; i<mip.dim(); i++) {
		out << i << ": " << mip.col_lower(i) << "\t <= ";
		if (mip.col_names[i]) out << mip.col_names[i]; else out << "col" << i;
		out << "\t <= " << mip.col_upper(i) << "\t " << mip.col_types[i] << endl;
	}
	out << "Rows (" << mip.rows() << "): " << endl;
	for (int c=0; c<mip.rows(); c++) {
		out << c << ": ";
		if (mip.row_names[c]) out << mip.row_names[c] << ": \t"; else out << "row " << c << ": \t";
		out << mip.row_lower[c] << "\t <= ";
		double coeff;
		for (int i=0; i<mip.dim(); i++) {
			if (!(coeff=mip.A(c,i))) continue;
			out << setiosflags(ios_base::showpos) << coeff << resetiosflags(ios_base::showpos);
			if (mip.col_names[i]) out << mip.col_names[i] << " ";
			else out << "col" << i << " ";
		}
		out << "\t <= " << mip.row_upper[c] << endl;
	}

	return out;
}
