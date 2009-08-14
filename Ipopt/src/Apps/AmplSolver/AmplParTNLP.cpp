// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-07-02


#include "AmplParTNLP.hpp"
#include "IpBlas.hpp"

#include "IpMpi.hpp"

/* AMPL includes */
#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"

#define C_TO_FORT_IX4(a,b,c,d) {a++; b++; c++; d++;}

using namespace Ipopt;


// This function divides an index space uniformly across p processors
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
partition_constraints(ASL_pfgh* asl, Index num_proc, Index proc_id,
                      Index &m_first, Index &m_last)
{
  Index i;
  Number tweight;
  Number *weight = new Number[n_con];
  Index *p_first, *p_last;
  Number avgw;

  DBG_ASSERT(asl);

  if (num_proc == 1) {
    m_first = 0;
    m_last = n_con-1;
    delete [] weight;
    return;
  }

  p_first = new Index[num_proc];
  p_last = new Index[num_proc];

  // calculate number of non-zeros in jacobian per constraint
  // set this number to 1 for linear constraints
  //
#if 1
  // AW: I'm not sure treating linear constraints extra is a good
  //     idea, we should rather go for load-balacing?
  tweight = 0.;
  for (i=0; i<nlc; i++) {
    cgrad *cg;
    Index nz;

    for (nz=0, cg=Cgrad[i]; cg; cg = cg->next) nz ++;
    tweight += nz+1;
    weight[i] = nz+1;
  }

  tweight += n_con - nlc;
  for (i=nlc; i<n_con; i++)
    weight[i] = 1;
#else
  tweight = 0.;
  for (i=0; i<nlc; i++) {
    cgrad *cg;
    Index nz;

    for (nz=0, cg=Cgrad[i]; cg; cg = cg->next) nz ++;
    tweight += nz+1;
    weight[i] = nz+1;
  }

  const Number scaling_factor_for_linear_constraints = 0.5;
  for (i=nlc; i<n_con; i++) {
    cgrad *cg;
    Index nz;

    for (nz=0, cg=Cgrad[i]; cg; cg = cg->next) nz ++;
    tweight += (nz+1)*scaling_factor_for_linear_constraints;
    weight[i] = (nz+1)*scaling_factor_for_linear_constraints;
  }
#endif

  avgw = tweight / num_proc;

  // Now let each processor have approximately the same weight
  Index cur_p = 0;
  Index cur_start = 0;
  Number cur_weight = 0.0;

  for (i=0; i<n_con && cur_p < num_proc; i++) {

    tweight -= weight[i];
    cur_weight += weight[i];

    if (cur_weight >= avgw || num_proc-cur_p >= n_con-i || i == n_con-1) {
      p_first[cur_p] = cur_start;
      p_last[cur_p] = i;
      cur_start = i+1;
      cur_p ++;
      cur_weight = 0.0;
      if (cur_p < num_proc) avgw = tweight / (num_proc-cur_p);
    }
  }
  // all constraints should have been allocated
  assert(i == n_con);

  for (i=cur_p; i<num_proc; i++) {
    p_first[i] = n_con;
    p_last[i] = n_con-1;
  }

  m_first = p_first[proc_id];
  m_last = p_last[proc_id];

  delete [] weight;
  delete [] p_first;
  delete [] p_last;
  return;
}


AmplParTNLP::AmplParTNLP(SmartPtr<AmplTNLP> amplobj)
    :
    amplobj_(amplobj),
    nnz_jac_g_(0),
    nnz_h_lag_(0),
    jac_map_(NULL)
{}

AmplParTNLP::~AmplParTNLP()
{
  if (jac_map_) delete [] jac_map_;
}

bool
AmplParTNLP::get_nlp_info(Index num_proc, Index proc_id,
                          Index& n, Index& n_first, Index& n_last,
                          Index& m, Index& m_first, Index& m_last,
                          Index& nnz_jac_g_part, Index& nnz_h_lag_part,
                          IndexStyleEnum& index_style)
{
  Index i;
  bool rval = true;
  ASL_pfgh* asl = amplobj_->AmplSolverObject();

  DBG_ASSERT(asl);

  if (proc_id != 0) amplobj_->set_obj_in_hessian(0);

  rval = amplobj_->get_nlp_info(n, m, nnz_jac_g_, nnz_h_lag_,
                                (TNLP::IndexStyleEnum&)index_style);
  if (!rval) return rval;

  calculate_offsets (num_proc, proc_id, n, n_first, n_last);
  if (m > 0)
    partition_constraints (asl, num_proc, proc_id, m_first, m_last);
  else {
    m_first = 0;
    m_last = -1;
  }

  if (num_proc == 1) {
    nnz_jac_g_part = nnz_jac_g_;
    nnz_h_lag_part = nnz_h_lag_;
  }
  else {
    // nonzeros in jacobian of constraints from m_first to m_last
    // TODO: use nonzeros in jacobian to figure out nonzeros in hessian
    Index nz = 0;
    for (i=m_first; i<=m_last; i++) {
      cgrad *cg;

      for (cg=Cgrad[i]; cg; cg = cg->next) nz ++;
    }
    nnz_jac_g_part = nz;
    nnz_h_lag_part = nnz_h_lag_;

#if 0 // This stuff is for hessian computation based on cgrad
    if (nnz_h_lag_ > 0 && proc_id > 0) {
      Index nv = 0;
      Index *var_exist = new Index[n];

      for (i=0; i<n; i++) var_exist[i] = 0;
      for (i=m_first; i<=m_last && i<nlc; i++) {
        cgrad *cg;

        for (cg=Cgrad[i]; cg; cg = cg->next) var_exist[cg->varno] = 1;
      }

      for (i=0; i<n; i++) if (var_exist[i]) nv ++;

      // if number of variables in all constraints < total
      if (nv < n) {
        Index *iRow = new Index[nnz_h_lag_];
        Index *jCol = new Index[nnz_h_lag_];

        rval = amplobj_->eval_h(n, NULL, false, 0.0, m, NULL, false, nnz_h_lag_, iRow, jCol, NULL);
        if (!rval) goto CLEANUP;

        for (i=0; i<nnz_h_lag_; i++)
          if (var_exist[iRow[i]-1] && var_exist[jCol[i]-1])
            hess_map_.push_back(i);

        if (hess_map_.size() < nnz_h_lag_) {
          nnz_h_lag_part = hess_map_.size();
        }
        else
          hess_map_.clear();

CLEANUP:
        delete [] iRow;
        delete [] jCol;
      }
    }
#endif
  }

  // AMPL's conval gets values of constraints with: n_conjac[0] <= i < n_conjac[1]
  n_conjac[0] = m_first;
  n_conjac[1] = m_last+1;

#if 1 // Hessian computation sparsity based on evaluating at random perturbation of starting point
  if (num_proc > 1 && nnz_h_lag_ > 0) {
    Index *iRow = new Index[nnz_h_lag_];
    Index *jCol = new Index[nnz_h_lag_];
    Number obj_value, *x = new Number[n];
    Number *g = new Number[m], *lambda = new Number[m];
    Number *val = new Number[nnz_h_lag_];

    for (i=0; i<n; i++) {
      if (havex0[i]) {
        x[i] = X0[i];
      }
      else {
        x[i] = 0.0;
      }
    }
    for (i=0; i<n; i++) {
      // ToDo: We need to make sure perturbed point is within bounds
      x[i] += 1e-3 * (IpRandom01() - 0.5) * Max(1., fabs(x[i]));
    }

    rval = amplobj_->eval_h(n, NULL, false, 0.0, m, NULL, false, nnz_h_lag_, iRow, jCol, NULL);
    if (!rval) goto CLEANUP2;

    for (i=0; i<m; i++) lambda[i] = 1. + 10.*IpRandom01();
    // should not be necessary:
    //for (i=0; i<nnz_h_lag_; i++) val[i] = 0.0;

    rval = eval_h(num_proc, proc_id, n, n_first, n_last, x, true, 1.0,
                  m, m_first, m_last, lambda, true, nnz_h_lag_,
                  NULL, NULL, val);
    if (!rval) goto CLEANUP2;

    hess_map_.clear();
    for (i=0; i<nnz_h_lag_; i++)
      if (val[i] != 0.0)
        hess_map_.push_back(i);

    if (hess_map_.size() < nnz_h_lag_) {
      nnz_h_lag_part = hess_map_.size();
    }
    else
      hess_map_.clear();

CLEANUP2:
    delete [] iRow;
    delete [] jCol;
    delete [] x;
    delete [] g;
    delete [] lambda;
    delete [] val;
  }
#endif

  // As AmplTNLP implements Fortran style, we need to too
  C_TO_FORT_IX4(n_first, n_last, m_first, m_last);

  return rval;
}


bool
AmplParTNLP::get_bounds_info(Index num_proc, Index proc_id,
                             Index n, Index n_first, Index n_last,
                             Number* x_l_part, Number* x_u_part,
                             Index m, Index m_first, Index m_last,
                             Number* g_l_part, Number* g_u_part)
{
  Index i;
  ASL_pfgh* asl = amplobj_->AmplSolverObject();

  DBG_ASSERT(asl);

  for (i=n_first; i<=n_last; i++) {
    x_l_part[i-n_first] = LUv[2*i];
    x_u_part[i-n_first] = LUv[2*i+1];
  }

  for (i=m_first; i<=m_last; i++) {
    g_l_part[i-m_first] = LUrhs[2*i];
    g_u_part[i-m_first] = LUrhs[2*i+1];
  }

  return true;
}

bool
AmplParTNLP::get_starting_point(Index num_proc, Index proc_id,
                                Index n, Index n_first, Index n_last,
                                bool init_x, Number* x_part,
                                bool init_z, Number* z_L_part, Number* z_U_part,
                                Index m, Index m_first, Index m_last,
                                bool init_lambda, Number* lambda_part)
{
  Index i;
  ASL_pfgh* asl = amplobj_->AmplSolverObject();

  DBG_ASSERT(asl);

  if (init_x) {
    for (i=n_first; i<=n_last; i++) {
      if (havex0[i])
        x_part[i-n_first] = X0[i];
      else
        x_part[i-n_first] = 0.0;
    }
  }

  // TODO: Could use BLAS at many places

  if (init_z) {
    // Modified for warm-start from AMPL
    DBG_ASSERT(IsValid(amplobj_->SuffixHandler()));
    const double* zL_init = amplobj_->SuffixHandler()->GetNumberSuffixValues("ipopt_zL_in", AmplSuffixHandler::Variable_Source);
    const double* zU_init = amplobj_->SuffixHandler()->GetNumberSuffixValues("ipopt_zU_in", AmplSuffixHandler::Variable_Source);

    if (zL_init) {
      for (i=n_first; i<=n_last; i++)
        z_L_part[i-n_first]=zL_init[i];
    }
    else {
      for (i=n_first; i<=n_last; i++)
        z_L_part[i-n_first]=1.0;
    }

    if (zU_init) {
      for (i=n_first; i<=n_last; i++)
        z_U_part[i-n_first]=zU_init[i];
    }
    else {
      for (i=n_first; i<=n_last; i++)
        z_U_part[i-n_first]=1.0;
    }
  }

  if (init_lambda) {
    for (i=m_first; i<=m_last; i++) {
      if (havepi0[i])
        lambda_part[i-m_first] = pi0[i];
      else
        lambda_part[i-m_first] = 0.0;
    }
  }

  return true;
}

bool
AmplParTNLP::eval_f(Index num_proc, Index proc_id,
                    Index n, Index n_first, Index n_last,
                    const Number* x, bool new_x, Number& obj_value)
{
  bool rval = amplobj_->eval_f(n, x, new_x, obj_value);
  if (proc_id != 0)
    obj_value = 0.0;

  return rval;
}

bool
AmplParTNLP::eval_grad_f(Index num_proc, Index proc_id,
                         Index n,  Index n_first, Index n_last,
                         const Number* x, bool new_x,
                         Number* grad_f_part)
{
  Number* gf = new Number[n];
  bool rval = true;

  rval = amplobj_->eval_grad_f(n, x, new_x, gf);
  if (rval) {
    for (Index i=n_first; i<=n_last; i++) {
      grad_f_part[i-n_first] = gf[i];
    }
  }

  delete [] gf;

  return rval;
}

bool
AmplParTNLP::eval_g(Index num_proc, Index proc_id,
                    Index n, const Number* x, bool new_x,
                    Index m, Index m_first, Index m_last,
                    Number* g_part)
{
  ASL_pfgh* asl = amplobj_->AmplSolverObject();

  DBG_ASSERT(asl);

  if (m == 0 || m_last < m_first) return true;

  // AMPL's conval gets values of constraints with: n_conjac[0] <= i < n_conjac[1]
  // These are already set in get_nlp_info(); repeat for safety
  n_conjac[0] = m_first;
  n_conjac[1] = m_last+1;

  // TODO: pass <= m
  return amplobj_->eval_g(n, x, new_x, m, g_part);
}

bool
AmplParTNLP::eval_jac_g(Index num_proc, Index proc_id,
                        Index n, const Number* x, bool new_x,
                        Index m, Index m_first, Index m_last,
                        Index nele_jac_part, Index* iRow_part,
                        Index *jCol_part, Number* values_part)
{
  bool rval = true;
  Index i;
  ASL_pfgh* asl = amplobj_->AmplSolverObject();

  if (m == 0 || m_last < m_first || nele_jac_part == 0) return true;

  if (values_part == NULL) {
    if (num_proc == 1) {
      DBG_ASSERT(nnz_jac_g_ == nele_jac_part);
      rval = amplobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g_, iRow_part, jCol_part, NULL);
    }
    else { // num_proc > 1
      cgrad *cg;
      Index nz = 0;

      if (jac_map_) delete jac_map_;
      jac_map_ = new Index[nele_jac_part];

      for (i=m_first; i<=m_last; i++) {
        for (cg=Cgrad[i]; cg; cg = cg->next) {
          iRow_part[nz] = i+1-m_first;
          jCol_part[nz] = cg->varno + 1;
          jac_map_[nz] = cg->goff;
          nz ++;
        }
      }
      DBG_ASSERT (nz == nele_jac_part);
    }
  }
  else {
    if (num_proc == 1) {
      DBG_ASSERT(nnz_jac_g_ == nele_jac_part);
      rval = amplobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g_, NULL, NULL, values_part);
    }
    else {
      DBG_ASSERT (jac_map_);
      Number* values = new Number[nnz_jac_g_];

      // AMPL's jacval gets jacobian of constraints with: n_conjac[0] <= i < n_conjac[1]
      // These are already set in get_nlp_info(); repeat for safety
      n_conjac[0] = m_first;
      n_conjac[1] = m_last+1;

      rval = amplobj_->eval_jac_g(n, x, new_x, m, nnz_jac_g_, NULL, NULL, values);
      if (!rval) {
        delete [] values;
        return rval;
      }

      for (i=0; i<nele_jac_part; i++)
        values_part[i] = values[jac_map_[i]];

      delete [] values;
    }
  }

  return rval;
}


bool AmplParTNLP::eval_h(Index num_proc, Index proc_id,
                         Index n, Index n_first, Index n_last,
                         const Number* x, bool new_x, Number obj_factor,
                         Index m, Index m_first, Index m_last,
                         const Number* lambda,
                         bool new_lambda, Index nele_hess_part,
                         Index* iRow_part, Index* jCol_part,
                         Number* values_part)
{
  bool rval = true;
  Index i;

  if (nele_hess_part == 0) return true;
  if (proc_id > 0 && (m == 0 || m_last < m_first)) return true;

  if (iRow_part && jCol_part && !values_part) {
    if (nele_hess_part == nnz_h_lag_)
      rval = amplobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag_, iRow_part, jCol_part, NULL);

    else {
      Index *iRow = new Index[nnz_h_lag_];
      Index *jCol = new Index[nnz_h_lag_];

      rval = amplobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag_, iRow, jCol, NULL);

      for (i=0; i<nele_hess_part; i++) {
        iRow_part[i] = iRow[hess_map_[i]];
        jCol_part[i] = jCol[hess_map_[i]];
      }

      delete [] iRow;
      delete [] jCol;
    }
  }

  else if (!iRow_part & !jCol_part && values_part) {
    // separate calls needed for one & multi processor case

    if (num_proc == 1) {
      DBG_ASSERT(nnz_h_lag_ == nele_hess_part);
      rval = amplobj_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nnz_h_lag_, NULL, NULL, values_part);
    }
    else {
      Number *lambda_part = new Number[m];
      memcpy(lambda_part, lambda, m*sizeof(Number));

      for (i=0; i<m_first; i++) lambda_part[i] = 0.0;
      for (i=m_last+1; i<m; i++) lambda_part[i] = 0.0;

      if (nele_hess_part == nnz_h_lag_) {
        rval = amplobj_->eval_h(n, x, new_x, obj_factor, m, lambda_part, new_lambda, nnz_h_lag_, NULL, NULL, values_part);
      }
      else {
        Number *values = new Number[nnz_h_lag_];

        rval = amplobj_->eval_h(n, x, new_x, obj_factor, m, lambda_part, new_lambda, nnz_h_lag_, NULL, NULL, values);

        for (i=0; i<nele_hess_part; i++)
          values_part[i] = values[hess_map_[i]];

        delete [] values;
      }

      delete [] lambda_part;
    }
  }
  else {
    DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
  }

  return rval;
}

void
AmplParTNLP::finalize_solution(SolverReturn status,
                               Index n, const Number* x, const Number* z_L, const Number* z_U,
                               Index m, const Number* g, const Number* lambda,
                               Number obj_value,
                               const IpoptData* ip_data,
                               IpoptCalculatedQuantities* ip_cq)
{
  int proc_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

  if (proc_id == 0)
    amplobj_->finalize_solution (status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
}
