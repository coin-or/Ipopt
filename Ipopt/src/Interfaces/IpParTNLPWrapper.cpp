// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-06-17

#include "IpParTNLPWrapper.hpp"

using namespace Ipopt;

static void calculate_offsets (Index p, Index p_id, Index n, Index &n_first, Index &n_last)
{
  if (p > 1) {
    int np = n/p, np_rem = n - np*p;

    if (p_id < np_rem) {
      np ++;
      n_first = p_id * np;
    }
    else
      n_first = p_id * np + np_rem;

    n_last = n_first + np-1;
  }
  else {
    n_first = 0;
    n_last = n-1;
  }
}

static void
partition_constraints(Index m, Index jac_nz, Index *iRow,
		      Index num_proc, Index proc_id, Index &m_first, Index &m_last)
{
  Index i;
  Index tweight = jac_nz;
  Index *weight = new Index[m];
  Index *p_first, *p_last;
  Number avgw;

  if (num_proc == 1) {
    m_first = 0;
    m_last = m-1;
    delete [] weight;
    return;
  }

  p_first = new Index[num_proc];
  p_last = new Index[num_proc];

  // calculate number of non-zeros in jacobian per constraint
  // assume C-STYLE indexing
  for (i=0; i<m; i++) weight[i] = 0;
  for (i=0; i<jac_nz; i++) weight[iRow[i]] ++;

  avgw = (Number)tweight / num_proc;

  // Now let each processor have approximately the same weight
  Index cur_p = 0;
  Index cur_start = 0;
  Number cur_weight = 0.0;

  for (i=0; i<m && cur_p < num_proc; i++) {

    tweight -= weight[i];
    cur_weight += weight[i];

    if (cur_weight >= avgw || num_proc-cur_p >= m-i || i == m-1) {
      p_first[cur_p] = cur_start;
      p_last[cur_p] = i;
      cur_start = i+1;
      cur_p ++;
      cur_weight = 0.0;
      if (cur_p < num_proc) avgw = (Number)tweight / (num_proc-cur_p);
    }
  }
  // all constraints should have been allocated
  assert(i == m);

  for (i=cur_p; i<num_proc; i++) {
    p_first[i] = m;
    p_last[i] = m-1;
  }

  m_first = p_first[proc_id];
  m_last = p_last[proc_id];

  delete [] weight;
  delete [] p_first;
  delete [] p_last;
  return;
}

ParTNLPWrapper::ParTNLPWrapper(SmartPtr<TNLP> tnlpobj)
    :
    tnlpobj_(tnlpobj),
    jac_map_(NULL),
    hess_map_(NULL)
{}

ParTNLPWrapper::~ParTNLPWrapper()
{
  delete [] jac_map_;
  delete [] hess_map_;
}

bool
ParTNLPWrapper::get_nlp_info(Index num_proc, Index proc_id,
                             Index& n, Index& n_first, Index& n_last,
                             Index& m, Index& m_first, Index& m_last,
                             Index& nnz_jac_g_part, Index& nnz_h_lag_part,
                             IndexStyleEnum& index_style)
{
  bool rval = true;

  rval = tnlpobj_->get_nlp_info(n, m, nnz_jac_g_, nnz_h_lag_,
                                (TNLP::IndexStyleEnum&)index_style);
  if (rval == false) return rval;

  if (index_style != C_STYLE) {
    printf("index_style = FORTRAN not yet cproperly implemented.\n");
    throw;
  }

  calculate_offsets (num_proc, proc_id, n, n_first, n_last);
  // calculate_offsets (num_proc, proc_id, m, m_first, m_last);

  if (m <= 0){
    m_first = 0;
    m_last = -1;
  }
  else if (num_proc > 1) {

  //printf("n = %d m = %d num_proc = %d proc_id = %d\n", n, m, num_proc, proc_id);
  //printf("n_first = %d n_last = %d\n", n_first, n_last);
  //printf("m_first = %d m_last = %d\n", m_first, m_last);

    delete [] jac_map_;
    jac_map_ = NULL;
    delete [] hess_map_;
    hess_map_ = NULL;

    Index* iRow = new Index[nnz_jac_g_];
    Index* jCol = new Index[nnz_jac_g_];

    rval = tnlpobj_->eval_jac_g(n, NULL, true, m, nnz_jac_g_, iRow, jCol, NULL);
    if (!rval) {
      delete [] iRow;
      delete [] jCol;
      return rval;
    }

    partition_constraints (m, nnz_jac_g_, iRow, num_proc, proc_id, m_first, m_last);

    Index nz = 0;
    // count the nonzeros in the m_first to m_last rows
    for (Index i=0; i<nnz_jac_g_; i++) {
      if (iRow[i] >= m_first && iRow[i] <= m_last) {
        nz++;
      }
    }
    nnz_jac_g_part = nz;
    // get the map
    jac_map_ = new Index[nnz_jac_g_part];
    nz = 0;
    for (Index i=0; i<nnz_jac_g_; i++) {
      if (iRow[i] >= m_first && iRow[i] <= m_last) {
        jac_map_[nz++] = i;
      }
    }

    delete [] iRow;
    iRow = NULL;
    delete [] jCol;
    jCol = NULL;

    iRow = new Index[nnz_h_lag_];
    jCol = new Index[nnz_h_lag_];

    rval = tnlpobj_->eval_h(n, NULL, true, 0, m, NULL, true, nnz_h_lag_, iRow, jCol, NULL);
    if (!rval) {
      delete [] iRow;
      delete [] jCol;
      return rval;
    }

    // count nonzeros
    nz = 0;
    for (Index i=0; i<nnz_h_lag_; i++) {
      if (jCol[i] >= n_first && jCol[i] <= n_last) {
        nz ++;
      }
    }
    nnz_h_lag_part = nz;
    // get the map
    hess_map_ = new Index[nnz_h_lag_part];
    nz = 0;
    for (Index i=0; i<nnz_h_lag_; i++) {
      if (jCol[i] >= n_first && jCol[i] <= n_last) {
        hess_map_[nz] = i;
        nz ++;
      }
    }

    delete [] iRow;
    delete [] jCol;
  }
  else {
    delete [] jac_map_;
    jac_map_ = NULL;
    delete [] hess_map_;
    hess_map_ = NULL;

    nnz_jac_g_part = nnz_jac_g_;
    nnz_h_lag_part = nnz_h_lag_;
  }

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
  if (rval) {
    for (i=n_first; i<=n_last; i++) {
      x_l_part[i-n_first] = x_l[i];
      x_u_part[i-n_first] = x_u[i];
    }

    for (i=m_first; i<=m_last; i++) {
      g_l_part[i-m_first] = g_l[i];
      g_u_part[i-m_first] = g_u[i];
    }
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
  Number *x=NULL, *z_L=NULL, *z_U=NULL, *lam=NULL;
  bool rval = true;

  x = new Number[n];
  z_L = new Number[n];
  z_U = new Number[n];
  lam = new Number[m];

  // TODO: Could use BLAS at many places

  rval = tnlpobj_->get_starting_point(n, init_x, x, init_z, z_L, z_U, m, init_lambda, lam);

  if (rval) {
    if (init_x) {
      for (Index i=n_first; i<=n_last; i++) {
        x_part[i-n_first] = x[i];
      }
    }
    if (init_z) {
      for (Index i=n_first; i<=n_last; i++) {
        z_L_part[i-n_first] = z_L[i];
        z_U_part[i-n_first] = z_U[i];
      }
    }

    if (init_lambda) {
      for (Index i=m_first; i<=m_last; i++) {
        lambda_part[i-m_first] = lam[i];
      }
    }
  }

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
  Number* gf = new Number[n];

  bool rval = tnlpobj_->eval_grad_f(n, x, new_x, gf);
  if (rval) {
    for (Index i=n_first; i<=n_last; i++) {
      grad_f_part[i-n_first] = gf[i];
    }
  }

  delete [] gf;

  return rval;
}

bool
ParTNLPWrapper::eval_g(Index num_proc, Index proc_id,
                       Index n, const Number* x, bool new_x,
                       Index m, Index m_first, Index m_last,
                       Number* g_part)
{
  Number* g = new Number[m];

  bool rval = tnlpobj_->eval_g(n, x, new_x, m, g);
  if (rval) {

    for (Index i=m_first; i<=m_last; i++) {
      g_part[i-m_first] = g[i];
    }
  }

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
  bool rval = true;

  if (values_part == NULL) {
    if (!jac_map_) {
      DBG_ASSERT(nnz_jac_g_ == nele_jac_part);
      rval = tnlpobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g_, iRow_part, jCol_part, NULL);
    }
    else {
      Index *iRow = new Index[nnz_jac_g_];
      Index *jCol = new Index[nnz_jac_g_];

      rval = tnlpobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g_, iRow, jCol, NULL);
      if (!rval) {
        delete [] iRow;
        delete [] jCol;
        return rval;
      }

      int nz = 0;
      for (Index i=0; i<nele_jac_part; i++) {
        iRow_part[i] = iRow[jac_map_[i]] - m_first;
        jCol_part[i] = jCol[jac_map_[i]];
      }

      delete [] iRow;
      delete [] jCol;
    }
  }
  else {
    if (!jac_map_) {
      DBG_ASSERT(nnz_jac_g_ == nele_jac_part);
      rval = tnlpobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g_, NULL, NULL, values_part);
    }
    else {
      Number* values = new Number[nnz_jac_g_];

      rval = tnlpobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g_, NULL, NULL, values);
      if (!rval) {
        delete [] values;
        return rval;
      }

      int nz = 0;
      for (Index i=0; i<nele_jac_part; i++) {
        values_part[i] = values[jac_map_[i]];
      }

      delete [] values;
    }
  }

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
  bool rval = true;

  if (values_part == NULL) {
    if (!hess_map_) {
      DBG_ASSERT(nnz_h_lag_ == nele_hess_part);
      rval = tnlpobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag_, iRow_part, jCol_part, NULL);
    }
    else {
      Index *iRow = new Index[nnz_h_lag_];
      Index *jCol = new Index[nnz_h_lag_];

      rval = tnlpobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag_, iRow, jCol, NULL);
      if (!rval) {
        delete [] iRow;
        delete [] jCol;
        return rval;
      }

      for (Index i=0; i<nele_hess_part; i++) {
        iRow_part[i] = iRow[hess_map_[i]];
        jCol_part[i] = jCol[hess_map_[i]];
      }

      delete [] iRow;
      delete [] jCol;
    }
  }
  else {
    if (!hess_map_) {
      DBG_ASSERT(nnz_h_lag_ == nele_hess_part);
      rval = tnlpobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag_, NULL, NULL, values_part);
    }
    else {
      Number *values = new Number[nnz_h_lag_];

      rval = tnlpobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag_, NULL, NULL, values);
      if (!rval) {
        delete [] values;
        return rval;
      }

      for (Index i=0; i<nele_hess_part; i++) {
        values_part[i] = values[hess_map_[i]];
      }

      delete [] values;
    }
  }

  return rval;
}

bool
ParTNLPWrapper::get_scaling_parameters(Index num_proc, Index proc_id,
                                       Number& obj_scaling,
                                       bool& use_x_scaling,
                                       Index n, Index n_first, Index n_last,
                                       Number* x_scaling_part,
                                       bool& use_g_scaling,
                                       Index m, Index m_first, Index m_last,
                                       Number* g_scaling_part)
{
  Number* x_scaling;
  Number* g_scaling;
  bool rval = true;

  x_scaling = new Number[n];
  g_scaling = new Number[m];

  rval =
    tnlpobj_->get_scaling_parameters(obj_scaling, use_x_scaling, n,
                                     x_scaling, use_g_scaling, m, g_scaling);

  if (rval) {
    if (use_x_scaling) {
      for (Index i=n_first; i<=n_last; i++) {
        x_scaling_part[i-n_first] = x_scaling[i];
      }
    }

    if (use_g_scaling) {
      for (Index i=m_first; i<=m_last; i++) {
        g_scaling_part[i-m_first] = g_scaling[i];
      }
    }
  }

  delete [] x_scaling;
  delete [] g_scaling;

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
