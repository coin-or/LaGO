// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

#include "usermatrix.h"

#ifdef ARPACK_AVAILABLE
#include "arpack++/arssym.h"
#endif

#ifdef LANCZOS_AVAILABLE
int UserMatrix::eig_lanczos(vector<dvector> &eig_vec, vector<double> &eig_val, Param *param) const {
  int nreig = eig_vec.size();
  int maxiter;
  Matrix Meig_val(nreig, 1, 0.);
  Matrix Meig_vec(dim(), nreig, 0.);
  for (int k=0; k<nreig; k++)
    for (int i=0; i<dim(); i++)
      Meig_vec(k*dim()+i)=eig_vec[k](i);

  Lanczpol lanczos;
  if (! param) maxiter=dim();
  else {
    maxiter=param->get_i("max_eig_iter_factor", 0)*dim();
    if (! maxiter) maxiter=param->get_i("max_eig_iter", dim());
  }
  lanczos.set_maxiter(maxiter);
//  lanczos.set_maxiter(param ? param->get_i("max_eig_iter",10000) : 10000);
//  lanczos.set_maxiter((param ? param->get_i("max_eig_iter_factor", 1) : 1) * dim());
  lanczos.set_relprec(param ? param->get_d("eig_tol", 0.001) : 0.001);
  lanczos.enable_stop_above(1e40);
  if (! (eig_vec[0]==0)) lanczos.set_nchebit(0);//-1);
  lanczos.set_nblockmult(-1);

//  lanczos.set_do_out(1);
//  lanczos.set_myout(cout);

  lanczos.compute((UserMatrix*)this, Meig_val, Meig_vec, nreig, 1);

  for (int i=0; i<nreig; i++) {
    eig_val[i]=-Meig_val[i];
    for (int j=0; j<dim(); j++) eig_vec[i][j]=Meig_vec[i,j];
  }

  return lanczos.get_err();
}

int UserMatrix::eig_lanczos(dvector &eig_vec, double &eig_val, Param *param) const {
  vector<dvector> evec(1); evec[0].resize(dim()); evec[0]=eig_vec;
  vector<double> eval(1);

  int ret=eig_lanczos(evec, eval, param);

  eig_val=eval[0];
  eig_vec=evec[0];

  return ret;
}
#endif // LANCZOS_AVAILABLE

#ifdef ARPACK_AVAILABLE
int UserMatrix::eig_arpack(dvector &eig_vec, double &eig_val, Param *param) {
  vector<double> EigVal(1, eig_val);
  vector<dvector> EigVec(1, eig_vec);

  int ret=eig_arpack(EigVec, EigVal, param);

  eig_val=EigVal[0];
  eig_vec=EigVec[0];

  return ret;
}

int UserMatrix::eig_arpack(vector<dvector>& eig_vec, vector<double>& eig_val, Param *param) {
  int nev=MAX(eig_vec.size(), 1);  // compute as least one eigenvalue, seems to work or do we need 2?

  int maxiter;
  if (! param) maxiter=dim();
  else {
    maxiter=param->get_i("max_eig_iter_factor", 0)*dim();
    if (! maxiter) maxiter=param->get_i("max_eig_iter", 10*dim());
  }
  maxiter=MAX(maxiter, 100*nev);

	ARSymStdEig<double, UserMatrix>* EigProb;

	if (eig_vec[0]==0) EigProb=new ARSymStdEig<double, UserMatrix>(dim(), nev, this, &UserMatrix::arpackmult, "SA", MIN(2*nev+1, dim()-1), param ? param->get_d("eig_tol", 0.001) : 0.001, maxiter);
	else EigProb=new ARSymStdEig<double, UserMatrix>(dim(), nev, this, &UserMatrix::arpackmult, "SA", MIN(2*nev+1, dim()-1), param ? param->get_d("eig_tol", 0.001) : 0.001, maxiter, (double*)(Pointer<double>)(eig_vec[0]));

	int ret=EigProb->FindEigenvectors();

  for (int i=0; i<MIN(ret, eig_vec.size()); i++) {
	  eig_val[i]=EigProb->Eigenvalue(i);
	  for (int j=0; j<dim(); j++) eig_vec[i][j]=EigProb->Eigenvector(i,j);
	}

		delete EigProb;

  return ret;
}
#endif // ARPACK_AVAILABLE

int UserMatrix::eig_ql(vector<dvector>& eig_vec, vector<double>& eig_val) const {
  DenseMatrix D(*this, true);
  return D.eig_ql(eig_vec, eig_val);
}

int UserMatrix::eig(dvector &eig_vec, double &eig_val, Param *param) {
  vector<dvector> EigVec(1, eig_vec);
  vector<double> EigVal(1, eig_val);

  int ret=eig(EigVec, EigVal, param);

  eig_vec=EigVec[0];
  eig_val=EigVal[0];

  return ret;
}

int UserMatrix::eig(vector<dvector> &eig_vec, vector<double> &eig_val, Param *param) {
	if (dim()<(param ? param->get_i("ql_lanczos_switch", 30) : 30)) {
		int ret=eig_ql(eig_vec, eig_val);
		if (!ret) return 0; // all fine, so return; else try arpack
	}

	Pointer<char> ea;
	if (param && (ea=param->get("eig algorithm"))) {
		if (!strcmp("lanczos", ea))
#ifdef LANCZOS_AVAILABLE
			return eig_lanczos(eig_vec, eig_val, param);
#else
			{ out_err << "Sorry, Lanczos not available. Try to set parameter 'eig algorithm' to 'arpack' or 'ql'. Aborting." << endl;
			  exit(-1);
			}
#endif
		if (!strcmp("ql", ea)) return eig_ql(eig_vec, eig_val);
	}
#ifdef ARPACK_AVAILABLE
	int ret = eig_vec.size() - eig_arpack(eig_vec, eig_val, param);
	if (ret) { // some problem in ARPACK, try QL
		out_err << "Problem in Eigenvalue computation with ARPACK, falling back to QL algorithm." << endl;
		ret=eig_ql(eig_vec, eig_val);
		if (!ret) ret=0;
	}
#else // use QL
  int ret = eig_ql(eig_vec, eig_val);
	if (!ret) ret=0;
#endif
  if (ret<0) return 0;
  return ret;  // so its the number of not converged eigenvalues
}

void UserMatrix::print(ostream &out) const {
  DenseMatrix D(*(UserMatrix*)this);
  out << D;
}

//-------------------------------------------------------------

void BlockMatrix::MultV(UserVector<double>& y, const UserVector<double>& x) const {
  y=0;
  for(int i=0; i<block.size(); i++)
    if (A[i]) {
      Pointer<UserVector<double> > yb(y.getemptycopy(block[i].size()));
      Pointer<UserVector<double> > xb(x.getcopy(block[i]));
      A[i]->MultV(*yb, *xb);
      y.set_block(*yb, block[i]);
    }
}

double BlockMatrix::yAx(const UserVector<double>& y, const UserVector<double>& x) const {
  double val=0;

  for(int i=0; i<block.size(); i++) {
    if (A[i]) {
	    Pointer<UserVector<double> > yb(y.getcopy(block[i]));
  	  Pointer<UserVector<double> > xb(x.getcopy(block[i]));
			val+=A[i]->yAx(*yb, *xb);
		}
  }
  return val;
}

#ifdef FILIB_AVAILABLE
void BlockMatrix::MultV(IntervalVector& y, const IntervalVector& x) const {
  y=interval<double>(0.);
  for(int i=0; i<block.size(); i++)
    if (A[i]) {
			assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A[i]));
      IntervalVector yb(block[i].size());
      IntervalVector xb(x, block[i]);
      ((const IntervalCompliantMatrix*)(const UserMatrix*)A[i])->MultV(yb, xb);
      y.set_block(yb, block[i]);
    }
}

interval<double> BlockMatrix::yAx(const IntervalVector& y, const IntervalVector& x) const {
	interval<double> val(0.);
  for(int i=0; i<block.size(); i++)
    if (A[i]) {
			assert(dynamic_cast<const IntervalCompliantMatrix*>((const UserMatrix*)A[i]));
      IntervalVector yb(y, block[i]);
      IntervalVector xb(x, block[i]);
      val+=((const IntervalCompliantMatrix*)(const UserMatrix*)A[i])->yAx(yb, xb);
    }
	return val;
}
#endif

//-------------------------------------------------------------

DenseMatrix::DenseMatrix(const UserMatrix& A_, bool allow_destroy_)
: ExtUserMatrix(A_.dim()), A(A_.dim(), A_.dim()), allow_destroy(allow_destroy_)
{
  dvector v(dim());
  SparseVector<double> e(dim());
  for (int i=0; i<dim(); i++) {
    e[i]=1;
    A_.MultV(v,e);
    for (int j=0; j<dim(); j++) A(i+1,j+1) = v(j);

    e.DelElement(i);
  }
}

void DenseMatrix::MultV(UserVector<double>& y_, const UserVector<double>& x_) const {
  TNT::Vector<double> x(x_.dim(), (const double*)(const Pointer<double>)x_);
//  for (int i=0; i<x_.dim(); i++) x[i]=x_(i);

  TNT::Vector<double> y(A*x);
//  y = A * x;
  for (int i=0; i<y_.dim(); i++) y_.SetElement(i,y[i]);
}

#ifdef FILIB_AVAILABLE
void DenseMatrix::MultV(IntervalVector& y, const IntervalVector& x) const {
	for (int i=0; i<dim(); i++) {
		y[i]=0.;
		for (int j=0; j<dim(); j++)
			y[i]+=A[i][j]*x(j);
	}
}
#endif

void DenseMatrix::MultV(double* y_, const double* x_) const {
  TNT::Vector<double> x(dim(), x_);
  TNT::Vector<double> y(dim());
  y = A * x;
  for (int i=0; i<dim(); i++) y_[i]=y[i];  // improve: memcpy
}

void DenseMatrix::set_random(const dvector& lambda) {
   assert(lambda.dim()==dim());
   A.newsize(dim(), dim());
   double c, s;
   dvector v(dim());

   // random vector with -1 < v[i] < 1
   v.set_random(-1., 1.);
   c = 2 / (v * v);
   // A = (I-c vv^t) * D
   for (int i = 0; i < dim(); i++)
     for (int j = 0; j < dim(); j++) {
			 if (i == j) (*this)(i,i) = lambda[i] * (1 - c * v[i] * v[i]);
			 else (*this)(i,j) = -lambda[j] * c * v[i] * v[j];
     }

   // A = A * (I-c vv^t)
   for (int i = 0; i < dim(); i++) {
     s = 0;
     for (int j = 0; j < dim(); j++)
       s += (*this)(i,j) * v[j];
     s *= c;
     for (int j = 0; j <= i; j++)
		   (*this)(i,j) -= s * v[j];
   }

   // assure exact symmetry
   // (obsolete for compactly stored sym matrices)
   for (int i = 0; i < dim(); i++)
     for (int j = 0; j < i; j++)
       (*this)(j,i) = (*this)(i,j);

}

int DenseMatrix::eig_ql(vector<dvector> &eig_vec, vector<double> &eig_val) const {
  double* diag=new double[dim()];
  double* subdiag=new double[dim()];

  double* content=*(double**)A;

  double* content_copy=NULL;
  if (allow_destroy) content_copy=content;
  else {
    content_copy=new double[dim()*dim()];
    for (int i=0; i<dim()*dim(); i++) content_copy[i]=content[i];
  }

  tred2(dim(), dim(), content, diag, subdiag, content_copy);
  int ret=imtql2(dim(), dim(), diag, subdiag, content_copy);
	
	if (ret) return ret;

  int k=0;
  for (int i=0; i<eig_vec.size(); i++) {
    eig_val[i]=diag[i];
    for (int j=0; j<dim(); j++) eig_vec[i][j]=content_copy[k++];
  }

  delete[] diag;
  delete[] subdiag;
  if (!allow_destroy) delete[] content_copy;

  return 0;
}

int DenseMatrix::eig_ql(dvector &eig_vec, double &eig_val) const {
  vector<dvector> vec(dim());
  for (int i=0; i<dim(); i++) vec[i].resize(dim());
  vector<double> val(dim());

  int ret=eig_ql(vec, val);
	
  if (ret==0) {  // all right, so first eigenvalue is the smallest one
    eig_vec=vec[0];
    eig_val=val[0];
    return 0;
  }
  // else: look for the smallest eigenvalue, which was found
  int min_index=-1;
	double min_value=INFINITY;
  for (int i=0; i<ret; i++) if (val[i]<min_value) { min_index=i; min_value=val[i]; }
	if (min_index>=0) {
	  eig_vec=vec[min_index];
  	eig_val=min_value;
	}
  return ret;
}

// from Helmbergs SBmethod:
/*     this subroutine is a translation of the algol procedure tred2, */
/*     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson. */
/*     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971). */
/*     this subroutine reduces a real symmetric matrix to a */
/*     symmetric tridiagonal matrix using and accumulating */
/*     orthogonal similarity transformations. */

void DenseMatrix::tred2(int nm, int n, double* a, double *d, double *e, double *z) const {
    /* System generated locals */
    Integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static Real f, g, h;
    static Integer i, j, k, l;
    static Real scale, hh;
    static Integer ii;

    /* Parameter adjustments */
    z_dim1 = nm;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --e;
    --d;
    a_dim1 = nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    for (i = 1; i <= n; ++i) {
	for (j = i; j <= n; ++j) {
	    z[j + i * z_dim1] = a[j + i * a_dim1];
	}
	d[i] = a[n + i * a_dim1];
    }

    if (n == 1) {
	goto L510;
    }

/*     .......... for i=n step -1 until 2 do -- .......... */
    i__1 = n;
    for (ii = 2; ii <= i__1; ++ii) {
	i = n + 2 - ii;
	l = i - 1;
	h = 0.;
	scale = 0.;
	if (l < 2) {
	    goto L130;
	}
/*     .......... scale row (algol tol then not needed) .......... */
	for (k = 1; k <= l; ++k) {
	    scale += fabs( d[k]);
	}

	if (scale != 0.) {
	    goto L140;
	}
L130:
	e[i] = d[l];

	for (j = 1; j <= l; ++j) {
	    d[j] = z[l + j * z_dim1];
	    z[i + j * z_dim1] = 0.;
	    z[j + i * z_dim1] = 0.;
	}

	goto L290;

L140:
	for (k = 1; k <= l; ++k) {
	    d[k] /= scale;
	    h += d[k] * d[k];
	}

	f = d[l];
	g = -d_sign(sqrt(h),f);
	e[i] = scale * g;
	h -= f * g;
	d[l] = f - g;

/*     .......... form a*u .......... */
	for (j = 1; j <= l; ++j) {
	    e[j] = 0.;
	}

	for (j = 1; j <= l; ++j) {
	    f = d[j];
	    z[j + i * z_dim1] = f;
	    g = e[j] + z[j + j * z_dim1] * f;

	    for (k = j+1; k <= l; ++k) {
		g += z[k + j * z_dim1] * d[k];
		e[k] += z[k + j * z_dim1] * f;
	    }

	    e[j] = g;
	}

/*     .......... form p .......... */
	f = 0.;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    e[j] /= h;
	    f += e[j] * d[j];
	}

	hh = f / (h + h);

/*     .......... form q .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    e[j] -= hh * d[j];
	}

/*     .......... form reduced a .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d[j];
	    g = e[j];
            Real *zpj=z+j+j * z_dim1;
            Real *ep=e+j;
            Real *dp=d+j;
	    for (k = j; k <= l; ++k) {
		//z[k + j * z_dim1] -= f * e[k] + g * d[k];
                (*zpj++)-= f*(*ep++)+g*(*dp++);
	    }

	    d[j] = z[l + j * z_dim1];
	    z[i + j * z_dim1] = 0.;
	}

L290:
	d[i] = h;
    }

/*     .......... accumulation of transformation matrices .......... */
    i__1 = n;
    for (i = 2; i <= i__1; ++i) {
	l = i - 1;
	z[n + l * z_dim1] = z[l + l * z_dim1];
	z[l + l * z_dim1] = 1.;
	h = d[i];
	if (h == 0.) {
	    goto L380;
	}

	for (k = 1; k <= l; ++k) {
	    d[k] = z[k + i * z_dim1] / h;
	}
        //mat_xeya(l,d+1,z+1+i*z_dim1,1./h);

	for (j = 1; j <= l; ++j) {
	    g = 0.;
	    for (k = 1; k <= l; ++k) {
	        g += z[k + i * z_dim1] * z[k + j * z_dim1];
	    }
            //g=mat_ip(l,z+1+i * z_dim1,z+1+j * z_dim1);

	    for (k = 1; k <= l; ++k) {
	        z[k + j * z_dim1] -= g * d[k];
	    }
            //mat_xpeya(l,z+1+j*z_dim1,d+1,-g);
	}

L380:
	for (k = 1; k <= l; ++k) {
	    z[k + i * z_dim1] = 0.;
	}
        //mat_xea(l,z+1+i * z_dim1,0.);

    }

L510:
    for (i = 1; i <= n; ++i) {
	d[i] = z[n + i * z_dim1];
	z[n + i * z_dim1] = 0.;
    }

    z[n + n * z_dim1] = 1.;
    e[1] = 0.;
}

// from Helmbergs SBmethod:
/*     this subroutine is a translation of the algol procedure imtql2, */
/*     num. math. 12, 377-383(1968) by martin and wilkinson, */
/*     as modified in num. math. 15, 450(1970) by dubrulle. */
/*     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971). */

/*     this subroutine finds the eigenvalues and eigenvectors */
/*     of a symmetric tridiagonal matrix by the implicit ql method. */
/*     the eigenvectors of a full symmetric matrix can also */
/*     be found if  tred2  has been used to reduce this */
/*     full matrix to tridiagonal form. */

int DenseMatrix::imtql2(int nm, int n, double *d, double *e, double *z) const {
    /* System generated locals */
    Integer z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static Real b, c, f, g;
    static Integer i, j, k, l, m;
    static Real p, r, s;
    static Integer ii, mml;
    static Real tst1, tst2;
    static Integer ierr;
    Integer iter_max=max(30,n);

    /* Parameter adjustments */
    z_dim1 = nm;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --e;
    --d;

    /* Function Body */
    ierr = 0;
    if (n == 1) {
	goto L1001;
    }

    i__1 = n;
    for (i = 2; i <= i__1; ++i) {
/* L100: */
	e[i - 1] = e[i];
    }

    e[n] = 0.;

    i__1 = n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
/*     .......... look for small sub-diagonal element .......... */
L105:
	i__2 = n;
	for (m = l; m <= i__2; ++m) {
	    if (m == n) {
		goto L120;
	    }
	    tst1 = fabs(d[m])+fabs(d[m+1]);
	    tst2 = tst1 + fabs(e[m]);
	    if (tst2 == tst1) {
		goto L120;
	    }
/* L110: */
	}

L120:
	p = d[l];
	if (m == l) {
	    goto L240;
	}
	if (j == iter_max) {
	    goto L1000;
	}
	++j;
/*     .......... form shift .......... */
	g = (d[l + 1] - p) / (e[l] * 2.);
	r = sqrt(g * g + 1.);
	g = d[m] - p + e[l] / (g + d_sign(r, g));
	s = 1.;
	c = 1.;
	p = 0.;
	mml = m - l;
/*     .......... for i=m-1 step -1 until l do -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i = m - ii;
	    f = s * e[i];
	    b = c * e[i];
	    r = sqrt(f * f + g * g);
	    e[i + 1] = r;
	    if (r == 0.) {
		goto L210;
	    }
	    s = f / r;
	    c = g / r;
	    g = d[i + 1] - p;
	    r = (d[i] - g) * s + c * 2. * b;
	    p = s * r;
	    d[i + 1] = g + p;
	    g = c * r - b;
/*     .......... form vector .......... */
	    i__3 = n;
            Real* zpi=z+i*z_dim1;
            Real* zpi1=zpi+z_dim1;
	    for (k = 1; k <= i__3; ++k) {
		//f = z[k + (i + 1) * z_dim1];
                f=*(++zpi1);
		//z[k + (i + 1) * z_dim1] = s * z[k + i * z_dim1] + c * f;
                *zpi1=s*(*(++zpi))+c*f;
		//z[k + i * z_dim1] = c * z[k + i * z_dim1] - s * f;
                *zpi*=c; *zpi-=s*f;
/* L180: */
	    }

/* L200: */
	}

	d[l] -= p;
	e[l] = g;
	e[m] = 0.;
	goto L105;
/*     .......... recover from underflow .......... */
L210:
	d[i + 1] -= p;
	e[m] = 0.;
	goto L105;
L240:
	;
    }
/*     .......... order eigenvalues and eigenvectors .......... */
    i__1 = n;
    for (ii = 2; ii <= i__1; ++ii) {
	i = ii - 1;
	k = i;
	p = d[i];

	i__2 = n;
	for (j = ii; j <= i__2; ++j) {
	    if (d[j] >= p) {
		goto L260;
	    }
	    k = j;
	    p = d[j];
L260:
	    ;
	}

	if (k == i) {
	    goto L300;
	}
	d[k] = d[i];
	d[i] = p;

	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    p = z[j + i * z_dim1];
	    z[j + i * z_dim1] = z[j + k * z_dim1];
	    z[j + k * z_dim1] = p;

/* L280: */
	}
        //mat_swap(i__2,z+1+i*z_dim1,z+1+k*z_dim1);

L300:
	;
    }

    goto L1001;
/*     .......... set error -- no convergence to an */
/*                eigenvalue after 30 iterations .......... */
L1000:
    ierr = l;
L1001:
    return ierr;

}

// ------------------------------ SparseMatrix ---------------------------------

SparseMatrix::SparseMatrix(const UserMatrix& A_, bool no_finish)
: rows_(A_.dim()), cols_(A_.dim()), val(0), row_ind(0), col_ptr(0), nz(-1)
{  // could be improved
  SparseVector<double> v(rows_);
  SparseVector<double> e(cols_);
  for (int col=0; col<cols_; col++) {
    e[col]=1;
    A_.MultV(v,e);
    for (int row=0; row<rows_; row++) AddElement(row, col, v(row));

    e.DelElement(col);
  }
 	if (!no_finish) finish();
}

void SparseMatrix::resize(int rows__, int cols__) {
	if (rows__==0 || cols__==0) {
		values.clear();
		if (val) {
			delete[] val; val=NULL;
			delete[] row_ind; row_ind=NULL;
			delete[] col_ptr; col_ptr=NULL;
		}
	} else {
		assert(val==NULL);
		if (rows__<rows() || cols__<cols())
			for (map<pair<int, int>, double>::iterator it(values.begin()); it!=values.end();)
				if (it->first.first>=cols__ || it->first.second>=rows__) {
					map<pair<int, int>, double>::iterator next(it); ++next;
					values.erase(it);
					it=next;
				} else ++it;
	}
	rows_=rows__;
	cols_=cols__;
}
		
SparseMatrix& SparseMatrix::operator+=(const UserMatrix& A_) {
#ifndef NO_SPARSEMATRIX_ASSERTS
	assert(rows_==A_.dim());
	assert(cols_==A_.dim());
#endif
	dvector v(rows_);
	SparseVector<double> e(cols_);
	for (int col=0; col<cols_; col++) {
		e[col]=1.;
		A_.MultV(v,e);
		for (int row=0; row<rows_; row++) AddToElement(row, col, v(row));
		
		e.DelElement(col);
	}

	return *this;
}

SparseMatrix& SparseMatrix::operator=(const double v) {
	if (val==NULL) 
		for (map<pair<int,int>, double>::iterator it(values.begin()); it!=values.end(); it++)
			it->second=v;
	else
		for (int i=0; i<nz; i++) val[i]=v;
		
	return *this;
}

SparseMatrix& SparseMatrix::operator*=(const double v) {
	if (val==NULL)
		for (map<pair<int,int>, double>::iterator it(values.begin()); it!=values.end(); it++)
			it->second*=v;
	else
		for (int i=0; i<nz; i++) val[i]*=v;
		
	return *this;
}

void SparseMatrix::set_block(const SparseMatrix& A, const ivector& indices) {
	assert(val==NULL);
	for (map<pair<int,int>, double>::const_iterator it(A.values.begin()); it!=A.values.end(); ++it)
		values.insert(pair<pair<int,int>, double>(pair<int,int>(indices(it->first.first), indices(it->first.second)), it->second));
}

void SparseMatrix::finish() {
  if (val) return;

  nz=values.size();

  val=new double[nz];
  row_ind=new int[nz];
  col_ptr=new int[cols_+1];

	*col_ptr=0;
	int i=0; int colnr=0;
  for (map<pair<int,int>, double>::iterator it(values.begin()); i<nz; it++, i++) {
  	val[i]=it->second;
    row_ind[i]=it->first.second;
    while (it->first.first > colnr) col_ptr[++colnr]=i;
	}
  while (colnr < cols_) col_ptr[++colnr]=nz;

  values.clear();
}

void SparseMatrix::print(ostream& out) const {
  out << "SparseMatrix: rows=" << rows_ << " cols=" << cols_;
  if (!val) out << " not finished" << endl;
  else {
    out << endl;
   	for (int col=0; col<cols_; col++)
	    for (int j=col_ptr[col]; j<col_ptr[col+1]; j++)
		   	out_log << row_ind[j] << "\t " << col << "\t " << val[j] << endl;
	}
}

void SparseMatrix::plot(char* filename) const {
  assert(filename);
  ofstream out(filename);

  for (int col=0; col<cols_; col++)
    for (int j=col_ptr[col]; j<col_ptr[col+1]; j++)
      out << row_ind[j] << '\t' << -col << '\t' << val[j] << endl;
}

// ----------------------------- SparseMatrix2 ---------------------------------

void SparseMatrix2::make_symmetric() {
	assert(val==NULL);
	for (map<pair<int,int>, double>::iterator it(values.begin()); it!=values.end(); it++) {
		map<pair<int,int>, double>::iterator it2(values.find(pair<int,int>(it->first.second, it->first.first)));
		if (it2!=values.end())
			it->second=it2->second=.5*(it->second+it2->second);
		else {
			it->second*=.5;
			values.insert(pair<pair<int,int>, double>(pair<int,int>(it->first.second, it->first.first), it->second));
		}
	}

}

void SparseMatrix2::set_random(int num_el, double max) {		
  if (val) delete val;
  if (row_ind) delete row_ind;
  if (col_ptr) delete col_ptr;
  values.clear();
  for (int i=0; i<dim(); i++) AddElement(i, i, ::random(-max, max));
	for (int i=0; i<num_el-dim(); i++) {
	  int row=::random(0, dim()-1);
	  int col=::random(0, dim()-1);
		double val=::random(-max, max);
		AddElement(row, col, val);
		AddElement(col, row, val);
	}
	finish();
}

