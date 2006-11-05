// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

#include "dualqqp.h"

void QqpDualFunc::init_ext() {
   // check if partial Lagrange problem is extended
   int i,j,j0;
   bool equal, lin;

   num_ext=0;
   n_q=0;
   for(i=0;i<num_blocks();i++)
   {
      lin=true;
      if (qqp.obj->A[i]) lin=false;
      else for (int l=0; l<qqp.con.size(); l++)
        if (qqp.con[l]->A[i]) {
          lin=false;
          break;
        }
      if (lin) block_type[i]=LIN;
      else {
        // check if bounds are symmetric
        equal=true;
        for(j=0;j<qqp.block[i].size();j++)
        {
  	   j0=qqp.block[i][j];
	   if(fabs(qqp.lower(j0)+qqp.upper(j0))>rtol)
	   {
	    equal=false; block_type[i]=EXTQUAD; break;
	   }
        }
	 // check b not zero
	if (equal)
	{
	  if(qqp.obj->b[i]!=NULL) block_type[i]=EXTQUAD;
	  else
	  {
	    block_type[i]=QUAD;
	    for(j=0;j<qqp.con.size();j++)
	    {
	       //	    cerr << "i: " << i << " j: " << j <<endl;
	       if(qqp.con[j]->b[i]!=NULL)
	       {
		  block_type[i]=EXTQUAD;
		  //	       cerr << "oj\n";
		  break;
	       }
	     }
	  }
	}
      }
      if (block_type[i]==EXTQUAD) i_ext[i]=num_ext++;

      if (block_type[i]!=LIN) {
        i_q[i]=n_q;
        n_q+=qqp.block[i].size();
      }
      else i_q[i]=-1;

     out_func_log << "lag probl " << i << " is ";
     switch(block_type[i]) {
       case LIN: out_func_log << "linear" << endl;
                 break;
       case QUAD: out_func_log << "quadratic" << endl;
                  break;
       case EXTQUAD: out_func_log << "extended quadratic" << endl;
     }
  }
}

// set type of ball-constraints

void QqpDualFunc::set_ball_con()
{
   for(int k=0;k<num_blocks();k++)
     for(int i=0;i<qqp.block[k].size();i++)
       if (ball_con[k]=(!qqp.discr[qqp.block[k][i]])) break;
}


// compute the modified value b and c

void QqpDualFunc::get_bc(double &c_,dvector &b_,SepQcFunc &q)
{
   c_=q.c;

   for (int k=0; k<num_blocks(); k++) {
     dvector b0(q.block[k].size(), 0);
     if (q.b[k]) b0=*q.b[k];
     if (block_type[k]!=LIN) {
       dvector m0(mid_point, q.block[k]);
       if (q.A[k]) c_+=q.A[k]->xAx(m0);
       if (q.b[k]) c_+=*q.b[k] * m0;

       if (q.A[k]) b0+=*q.A[k] * m0;
       b0 = *W[k] * b0;
     }
     b_.set_block(b0, q.block[k]);
   }
}


void QqpDualFunc::set_bc()
{
   int i;
   for(i=0;i<qqp.con.size();i++)
      get_bc(c_con[i],b_con[i],*qqp.con[i]);
   get_bc(c_obj,b_obj,*qqp.obj);
}


void QqpDualFunc::init() {
   int k,i;
   out_func_log << "QqpDualFunc::init" << endl;

   init_ext();

   // set dual bounds for quadratic box constraints (only for nonlinear variables)
   for (k=0; k<num_blocks(); k++)
     if (block_type[k]!=LIN)
       for (i=0; i<qqp.block[k].size(); i++)
         dual_bounds.push_back(qqp.discr[qqp.block[k][i]] ? -INFINITY : 0);

   // set dual bounds for additional extended quadratic box constraints
   for(i=0;i<num_ext; i++) dual_bounds.push_back(-INFINITY);

   // set dual bounds of original constraints
   for(i=0;i<qqp.con.size();i++) dual_bounds.push_back(qqp.con_eq[i] ? -INFINITY : 0);

   dim_=dual_bounds.size();

   subgrad.resize(dim());
   dual_point.resize(dim());

   int nr_eig=param.get_i("nr of eig vec", 5);
   for(k=0; k<num_blocks(); k++) nr_eig=MIN(nr_eig, qqp.block[k].size());
   for(k=0; k<num_blocks(); k++)
   {
      eig_vec[k].resize(nr_eig);
      eig_val[k].resize(nr_eig);
      for (int l=0; l<nr_eig; l++) {
        eig_vec[k][l].resize(qqp.block[k].size()+(block_type[k]==EXTQUAD ? 1 : 0));
//        eig_vec[k][l]=1.0; // start vector for the eigenvalue computation
      }
   }
   for(i=0;i<qqp.con.size();i++) b_con[i].resize(primal_dim());
   for(k=0;k<num_blocks();k++) {
     switch(block_type[k]) {
       case LIN: A_lag[k]=NULL; break;
       case QUAD: A_lag[k]=new QqpMatrix(*this,k); break;
       case EXTQUAD: A_lag[k]=new QqpExtMatrix(*this,k); break;
     }
     b_lag[k].resize(qqp.block[k].size());
   }

   x_lag.resize(primal_dim()+num_ext);

   mid_point=0.5*(qqp.lower+qqp.upper);

   dvector w(0.5*(qqp.lower-qqp.upper));
   for (int k=0; k<num_blocks(); k++) {
     W.A[k]=new DiagMatrix(w.getcopy(qqp.block[k]));
     radius[k]=sqrt(qqp.block[k].size()+(block_type[k]==EXTQUAD ? 1 : 0));
   }

   set_bc();
   set_ball_con();
}


//compute value of the dual function

int QqpDualFunc::set_dual_val() const {
   // compute c_lag
   c_lag=c_obj;
   for (int i=0; i<qqp.con.size(); i++) c_lag+=dual_point[i+num_qbox()]*c_con[i];
   for (int i=0; i<num_qbox(); i++) c_lag-=dual_point[i];

   // compute b_lag
   for (int k=0; k<num_blocks(); k++) {
     if (block_type[k]!=LIN) {
       b_lag[k]=b_obj(qqp.block[k]);
       for (int i=0; i<qqp.con.size(); i++) b_lag[k]+=dual_point[i+num_qbox()] * b_con[i](qqp.block[k]);
     }
     else {
       if (qqp.obj->b[k]) b_lag[k]=*qqp.obj->b[k];
       else b_lag[k]=0;
       for (int i=0; i<qqp.con.size(); i++)
         if (qqp.con[i]->b[k]) b_lag[k]+=dual_point[i+num_qbox()] * *qqp.con[i]->b[k];
     }
   }
   // compute dual_val
   dual_val=c_lag;
   int ret=0;
   for(int k=0;k<num_blocks();k++)
     if (block_type[k]!=LIN) {

			 clock_t e_time=clock();
			 int ret0=A_lag[k]->eig(eig_vec[k], eig_val[k], &param);  // compute eigenvalue
			 eig_time+=(double)(clock()-e_time)/(double)CLOCKS_PER_SEC;
       if (ret0) {
         out_err << " Got an error from eigenvalue computation: " << ret0 << endl;
         out_func_log << endl;
         ret+=ret0;
       }
       else out_func_log << " min_eig_val=" << eig_val[k][0] << endl;

       if (ball_con[k]) dual_val+=MIN(0, radius[k]*radius[k]*eig_val[k][0]);
       else dual_val+=radius[k]*radius[k]*eig_val[k][0];
     }
     else
       for (int i=0; i<qqp.block[k].size(); i++)
         dual_val+=2*MIN(qqp.lower(qqp.block[k][i]) * b_lag[k][i], qqp.upper(qqp.block[k][i]) * b_lag[k][i]);

   return ret;
}

double QqpDualFunc::eval_mod_con(int con_nr) const {
  double val=c_con[con_nr];

  for (int k=0; k<num_blocks(); k++) {
    dvector x0(x_lag, qqp.block[k]);
    if (block_type[k]!=LIN) {
      if (qqp.con[con_nr]->A[k]) val+= qqp.con[con_nr]->A[k]->xAx(*W[k] * x0);
      if (block_type[k]==EXTQUAD) val+= 2 * x_lag[qqp.dim()+i_ext[k]] * x0 * b_con[con_nr](qqp.block[k]);
    }
    else if (qqp.con[con_nr]->b[k]) val+= 2 * *qqp.con[con_nr]->b[k] * x0;
  }
  return val;
}

double QqpDualFunc::eval_mod_obj() const {
  double val=c_obj;
  for (int k=0; k<num_blocks(); k++) {
    dvector x0(x_lag, qqp.block[k]);
    if (block_type[k]!=LIN) {
      if (qqp.obj->A[k]) val+= qqp.obj->A[k]->xAx(*W[k] * x0);
      if (block_type[k]==EXTQUAD) val+= 2 * x_lag[qqp.dim()+i_ext[k]] * x0 * b_obj(qqp.block[k]);
    }
    else if (qqp.obj->b[k]) val+= 2 * *qqp.obj->b[k] * x0;
  }
  return val;
}

// set x_lag and compute subgradient
void QqpDualFunc::set_subgrad() const {
  x_lag=0;
  for (int k=0; k<num_blocks(); k++)
    switch (block_type[k]) {
      case QUAD: if (!(eig_val[k][0]>=0 && ball_con[k])) x_lag.set_block(radius[k]*eig_vec[k][0], qqp.block[k]);
                 break;
      case EXTQUAD: if (!(eig_val[k][0]>=0 && ball_con[k])) {
      		      x_lag.set_block(radius[k]*eig_vec[k][0](0,qqp.block[k].size()-1), qqp.block[k]);
                      x_lag[i_ext[k]+primal_dim()]=radius[k]*eig_vec[k][0](qqp.block[k].size());
              	    }
                    break;
      case LIN: for (int i=0; i<qqp.block[k].size(); i++)
                  x_lag[qqp.block[k][i]] =
                    ( qqp.lower(qqp.block[k][i])*b_lag[k][i] < qqp.upper(qqp.block[k][i])*b_lag[k][i] )
		    ? qqp.lower(qqp.block[k][i])
		    : qqp.upper(qqp.block[k][i]);
    }

  // subgradient components for the quadratic box constraints
  for (int k=0; k<num_blocks(); k++)
    if (block_type[k]!=LIN) {
      for (int i=0; i<qqp.block[k].size(); i++)
        subgrad[i_q[k]+i] = x_lag[qqp.block[k][i]] * x_lag[qqp.block[k][i]] - 1;
      if (block_type[k]==EXTQUAD) subgrad[n_q + i_ext[k]] = x_lag[primal_dim()+i_ext[k]] * x_lag[primal_dim()+i_ext[k]] - 1;
    }

  // subgradient components for the original constraints
  for (int i=0; i<qqp.con.size(); i++) subgrad[num_qbox()+i]=eval_mod_con(i);

// test result:
//   zg=dual_point*subgrad;
//   out_val << "Obj: " << obj_val  << " z*g: " << zg << endl;
//   out_val << "******************* dual_val-(obj_val+zg): " << dual_val-(obj_val+zg) << endl;

}

int QqpDualFunc::valgrad(double& val, UserVector<double>& g, const UserVector<double>& z) const {
   assert(g.dim()==dim());
   assert(z.dim()==dim());

   dual_point=z;

   int ret=set_dual_val();
   set_subgrad();

// check subgrad and dualval
/*   double dual_val2 = eval_mod_obj() + subgrad * dual_point;
   out_dualval << " x_lag: " << x_lag;
   out_dualval << " dual_val2: " << dual_val2 << " = " << eval_mod_obj() << " + " << subgrad * dual_point << endl;
   out_dualval << endl;
*/

   g=subgrad;
   val=dual_val;
   return ret;
}

dvector QqpDualFunc::get_dual_point() const {
  dvector dpnt(qqp.con.size()+primal_dim());
  for (int k=0; k<num_blocks(); k++)
    if (block_type[k]!=LIN) {
			double val=MIN(0, eig_val[k][0]);
      for (int i=0; i<qqp.block[k].size(); i++)
				dpnt[qqp.block[k][i]]=dual_point[i_q[k]+i]-val;
    }
  for (int i=0; i<qqp.con.size(); i++)
		dpnt[primal_dim()+i]=dual_point[num_qbox()+i];

  return dpnt;
};

dvector QqpDualFunc::get_orig_point(int use_eig_vec) const {
  dvector sol(primal_dim());

  for (int k=0; k<num_blocks(); k++)
    switch (block_type[k]) {
      case QUAD: case EXTQUAD:
							if (!(eig_val[k][use_eig_vec]>=0 && ball_con[k]))	sol.set_block(radius[k]*eig_vec[k][use_eig_vec], qqp.block[k]);
                 break;
      case LIN: for (int i=0; i<qqp.block[k].size(); i++)
                  sol[qqp.block[k][i]] =
                    ( qqp.lower(qqp.block[k][i])*b_lag[k][i] < qqp.upper(qqp.block[k][i])*b_lag[k][i] )
		    ? qqp.lower(qqp.block[k][i])
		    : qqp.upper(qqp.block[k][i]);
    }
  return sol;
}

void QqpDualFunc::print(ostream &out) const {
   out  << "dual_dim: " << dim();

   for (int k=0; k<num_blocks(); k++) {
     out << endl;
     out << "block " << k << ": " << endl;

     out << "W: "; if (W.A[k]) out << *W.A[k]; else out << endl;
     out << "radius: " << radius[k] << endl;

     out << "A_lag: "; if (A_lag[k]) out << *A_lag[k]; else out << endl;
     out << "b_lag: " << b_lag[k];
     out << "c_lag: " << c_lag << endl;
     out << endl;

     out << "A in obj: "; if (qqp.obj->A[k]) out << *qqp.obj->A[k]; else out << endl;
     out << "b in obj: "; if (qqp.obj->b[k]) out << *qqp.obj->b[k]; else out << endl;
     out << "b_obj: " << b_obj;
     out << "c_obj: " << c_obj << endl;

     for (int i=0; i<qqp.con.size(); i++) {
       out << endl;
       out << "A in con " << i << ": "; if (qqp.con[i]->A[k]) out << *qqp.con[i]->A[k]; else out << endl;
       out << "b in con " << i << ": "; if (qqp.con[i]->b[k]) out << *qqp.con[i]->b[k]; else out << endl;
       out << "b_con(" << i << "): " << b_con[i];
       out << "c_con(" << i << "): " << c_con[i] << endl;
     }
   }
}


// ------------------------------------ QqpMatrix -------------------------------------

void QqpMatrix::MultV(dvector &y,const dvector &x) const {
   dvector x0(x.dim());
   d.W[b_num]->MultV(x0,x);

   if (d.qqp.obj->A[b_num]) d.qqp.obj->A[b_num]->MultV(y,x0);
   else y=0;

   int con_size=d.qqp.con.size();
   vector<Pointer<SepQcFunc> >::iterator it(d.qqp.con.begin());
   for (int i=0; i<con_size; it++, i++)
     if ((*it)->A[b_num]) (*it)->A[b_num]->AddMult(y, x0, d.dual_point(d.num_qbox()+i));
   y=*d.W[b_num]*y;

   y+=d.dual_point(d.i_q[b_num], d.i_q[b_num]+dim()-1).diagmult(x);

   d.num_matmult++;
}

// ------------------------------------ QqpExtMatrix ----------------------------------

void QqpExtMatrix::MultV(dvector &y,const dvector &x) const {
   dvector x0(x, 0, dim()-2);
   dvector y0(dim()-1);

   y[dim()-1]=d.b_lag[b_num]*x0 + d.dual_point[d.n_q+d.i_ext[b_num]]*x[dim()-1];

   x0=*d.W[b_num]*x0;
   if (d.qqp.obj->A[b_num]) d.qqp.obj->A[b_num]->MultV(y0,x0);
   else y0=0;
   int con_size=d.qqp.con.size();
   vector<Pointer<SepQcFunc> >::iterator it(d.qqp.con.begin());
   for (int i=0; i<con_size; it++, i++)
     if ((*it)->A[b_num]) (*it)->A[b_num]->AddMult(y0, x0, d.dual_point[d.num_qbox()+i]);
   y0=*d.W[b_num]*y0;

   y0+=d.dual_point(d.i_q[b_num], d.i_q[b_num]+dim()-2).diagmult(x(0, dim()-2));

   y0.AddMult(x[dim()-1], d.b_lag[b_num]);

   for (int i=0; i<dim()-1; i++) y[i]=y0[i];

   d.num_matmult++;
}

