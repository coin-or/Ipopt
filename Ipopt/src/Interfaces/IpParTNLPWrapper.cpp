// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-06-17

#include "IpParTNLPWrapper.hpp"

using namespace Ipopt;

void calculate_offsets (Index p, Index p_id, Index n, Index &n_first, Index &n_last)
{
  if (p > 1){
    int np = n/p, np_rem = n - np*p;

    if (p_id < np_rem){
      np ++;
      n_first = p_id * np;
    }
    else
      n_first = p_id * np + np_rem; 

    n_last = n_first + np-1;
  }
  else{
    n_first = 0;
    n_last = n-1;
  }
}

bool
ParTNLPWrapper::get_nlp_info(Index num_proc, Index proc_id,
			     Index& n, Index& n_first, Index& n_last,
			     Index& m, Index& m_first, Index& m_last,
			     Index& nnz_jac_g_part, Index& nnz_h_lag_part,
			     TNLP::IndexStyleEnum& index_style)
{
  Index nnz_jac_g, nnz_h_lag;
  bool rval = true;

  rval = tnlpobj_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  if (rval == false) return rval;

  if (index_style != TNLP::C_STYLE){
    n --;
    m --;
  }

  calculate_offsets (num_proc, proc_id, n, n_first, n_last);
  calculate_offsets (num_proc, proc_id, m, m_first, m_last);

  if (index_style != TNLP::C_STYLE){
    n_first ++;
    n_last ++;
    m_first ++;
    m_last ++;
    n ++;
    m ++;
  }

  Index *iRow = new Index[nnz_jac_g];
  Index *jCol = new Index[nnz_jac_g];

  if (num_proc > 1){
    int i, nz=0;

    tnlpobj_->eval_jac_g(n, NULL, true, m, nnz_jac_g, iRow, jCol, NULL);

    nz = 0;
    for (i=0; i<nnz_jac_g; i++)
      if (iRow[i] >= m_first && iRow[i] <= m_last) nz ++;
    nnz_jac_g_part = nz;

    tnlpobj_->eval_h(n, NULL, true, 0, m, NULL, true, nnz_h_lag, iRow, jCol, NULL);
    nz = 0;
    for (i=0; i<nnz_h_lag; i++)
      if (jCol[i] >= n_first && jCol[i] <= n_last) nz ++;
    nnz_h_lag_part = nz;

  }
  else{
    nnz_jac_g_part = nnz_jac_g;
    nnz_h_lag_part = nnz_h_lag;
  }

  delete [] iRow;
  delete [] jCol;

  return rval;
}

bool
ParTNLPWrapper::get_bounds_info(Index num_proc, Index proc_id,
				Index n, Index n_first, Index n_last,
				Number* x_l_part, Number* x_u_part,
				Index m, Index m_first, Index m_last,
				Number* g_l_part, Number* g_u_part)
{
  int i;
  Number *x_l=NULL, *x_u=NULL, *g_l=NULL, *g_u=NULL;
  bool rval = true;

  x_l = new Number[n];
  x_u = new Number[n];
  g_l = new Number[m];
  g_u = new Number[m];

  rval = tnlpobj_->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
  if (rval == false) return rval;

  for (i=n_first; i<=n_last; i++){
    x_l_part[i-n_first] = x_l[i];
    x_u_part[i-n_first] = x_u[i];
  }

  for (i=m_first; i<=m_last; i++){
    g_l_part[i-m_first] = g_l[i];
    g_u_part[i-m_first] = g_u[i];
  }

  delete [] x_l;
  delete [] x_u;
  delete [] g_l;
  delete [] g_u;

  return rval;
}

bool
ParTNLPWrapper::get_starting_point(Index num_proc, Index proc_id,
				   Index n, Index n_first, Index n_last,
				   bool init_x, Number* x_part,
				   bool init_z, Number* z_L_part, Number* z_U_part,
				   Index m, Index m_first, Index m_last,
				   bool init_lambda, Number* lambda_part)
{
  int i;
  Number *x=NULL, *z_L=NULL, *z_U=NULL, *lam=NULL;
  bool rval = true;

  x = new Number[n];
  z_L = new Number[n];
  z_U = new Number[n];
  lam = new Number[m];

  rval = tnlpobj_->get_starting_point(n, init_x, x, init_z, z_L, z_U, m, init_lambda, lam);
  if (rval == false) return rval;

  if (init_x)
    for (i=n_first; i<=n_last; i++)
      x_part[i-n_first] = x[i];

  if (init_z)
    for (i=n_first; i<=n_last; i++){
      z_L_part[i-n_first] = z_L[i];
      z_U_part[i-n_first] = z_U[i];
    }

  if (init_lambda)
    for (i=m_first; i<=m_last; i++)
      lambda_part[i-m_first] = lam[i];

  delete [] x;
  delete [] z_L;
  delete [] z_U;
  delete [] lam;

  return rval;
}

bool
ParTNLPWrapper::eval_f(Index num_proc, Index proc_id,
		       Index n, Index n_first, Index n_last,
		       const Number* x, bool new_x, Number& obj_value)
{
  if (proc_id == 0)
    return tnlpobj_->eval_f(n, x, new_x, obj_value);
  else
    obj_value = 0.0;

  return true;
}

bool
ParTNLPWrapper::eval_grad_f(Index num_proc, Index proc_id,
			    Index n,  Index n_first, Index n_last,
			    const Number* x, bool new_x,
			    Number* grad_f_part)
{
  int i;
  Number *gf=NULL;
  bool rval = true;

  gf = new Number[n];

  rval = tnlpobj_->eval_grad_f(n, x, new_x, gf);
  if (rval == false) return rval;

  for (i=n_first; i<=n_last; i++)
    grad_f_part[i-n_first] = gf[i];

  delete [] gf;

  return rval;
}

bool
ParTNLPWrapper::eval_g(Index num_proc, Index proc_id,
		       Index n, const Number* x, bool new_x,
		       Index m, Index m_first, Index m_last,
		       Number* g_part)
{
  int i;
  Number *g=NULL;
  bool rval = true;

  g = new Number[m];

  rval = tnlpobj_->eval_g(n, x, new_x, m, g);
  if (rval == false) return rval;

  for (i=m_first; i<=m_last; i++)
    g_part[i-m_first] = g[i];

  delete [] g;

  return rval;
}

bool
ParTNLPWrapper::eval_jac_g(Index num_proc, Index proc_id,
			   Index n, const Number* x, bool new_x,
			   Index m, Index m_first, Index m_last,
			   Index nele_jac_part, Index* iRow_part,
			   Index *jCol_part, Number* values_part)
{
  int i;
  TNLP::IndexStyleEnum index_style;
  Index nnz_jac_g, nnz_h_lag, tn, tm;
  bool rval = true;

  rval = tnlpobj_->get_nlp_info(tn, tm, nnz_jac_g, nnz_h_lag, index_style);
  if (rval == false) return rval;

  Index *iRow = new Index[nnz_jac_g];
  Index *jCol = new Index[nnz_jac_g];

  rval = tnlpobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g, iRow, jCol, NULL);
  if (rval == false) return rval;

  if (values_part == NULL){
    int nz = 0;
    for (i=0; i<nnz_jac_g; i++)
      if (iRow[i] >= m_first && iRow[i] <= m_last){
	iRow_part[nz] = iRow[i];
	jCol_part[nz] = jCol[i];
	nz ++;
      }
    assert (nz == nele_jac_part);
  }
  else{
    Number *values = new Number[nnz_jac_g];

    rval = tnlpobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g, iRow, jCol, values);
    if (rval == false) return rval;

    int nz = 0;
    for (i=0; i<nnz_jac_g; i++)
      if (iRow[i] >= m_first && iRow[i] <= m_last){
	values_part[nz] = values[i];
	nz ++;
      }
    assert (nz == nele_jac_part);
    
    delete [] values;
  }

  delete [] iRow;
  delete [] jCol;

  return rval;
}


bool ParTNLPWrapper::eval_h(Index num_proc, Index proc_id,
			    Index n, Index n_first, Index n_last,
			    const Number* x, bool new_x, Number obj_factor,
			    Index m, Index m_first, Index m_last,
			    const Number* lambda,
			    bool new_lambda, Index nele_hess_part,
			    Index* iRow_part, Index* jCol_part,
			    Number* values_part)
{
  int i;
  TNLP::IndexStyleEnum index_style;
  Index nnz_jac_g, nnz_h_lag, tn, tm;
  bool rval = true;

  rval = tnlpobj_->get_nlp_info(tn, tm, nnz_jac_g, nnz_h_lag, index_style);
  if (rval == false) return rval;

  Index *iRow = new Index[nnz_h_lag];
  Index *jCol = new Index[nnz_h_lag];

  rval = tnlpobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag, iRow, jCol, NULL);
  if (rval == false) return rval;

  if (values_part == NULL){
    int nz = 0;
    for (i=0; i<nnz_h_lag; i++)
      if (jCol[i] >= n_first && jCol[i] <= n_last){
	iRow_part[nz] = iRow[i];
	jCol_part[nz] = jCol[i];
	nz ++;
      }
    assert (nz == nele_hess_part);
  }
  else{
    Number *values = new Number[nnz_h_lag];

    rval = tnlpobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag, iRow, jCol, values);
    if (rval == false) return rval;

    int nz = 0;
    for (i=0; i<nnz_h_lag; i++)
      if (jCol[i] >= n_first && jCol[i] <= n_last){
	values_part[nz] = values[i];
	nz ++;
      }
    assert (nz == nele_hess_part);
    
    delete [] values;
  }

  delete [] iRow;
  delete [] jCol;

  return rval;

}

void ParTNLPWrapper::finalize_solution(SolverReturn status,
			 Index n, const Number* x, const Number* z_L, const Number* z_U,
			 Index m, const Number* g, const Number* lambda,
			 Number obj_value,
			 const IpoptData* ip_data,
			 IpoptCalculatedQuantities* ip_cq)
{
  tnlpobj_->finalize_solution (status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
}
