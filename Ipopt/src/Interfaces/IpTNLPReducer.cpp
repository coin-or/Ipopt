// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                  IBM    2008-08-10

#include "IpTNLPReducer.hpp"

#include <algorithm>
#include <limits>

namespace Ipopt
{
  TNLPReducer::TNLPReducer(TNLP& tnlp,
                           Index n_g_skip, const Index* index_g_skip,
                           Index n_xL_skip, const Index* index_xL_skip,
                           Index n_xU_skip, const Index* index_xU_skip,
                           Index n_x_fix, const Index* index_x_fix)
      :
      tnlp_(&tnlp),
      n_g_skip_(n_g_skip),
      index_g_skip_(NULL),
      g_keep_map_(NULL),
      m_reduced_(-1),
      jac_g_skipped_(NULL),
      n_xL_skip_(n_xL_skip),
      index_xL_skip_(NULL),
      n_xU_skip_(n_xU_skip),
      index_xU_skip_(NULL),
      n_x_fix_(n_x_fix),
      index_x_fix_(NULL)
  {
    index_g_skip_ = new Index[n_g_skip_+1];
    for (Index i=0; i<n_g_skip_; i++) {
      index_g_skip_[i] = index_g_skip[i];
    }
    std::sort(index_g_skip_, index_g_skip_+n_g_skip_);
    index_g_skip_[n_g_skip_] = -1;

    index_xL_skip_ = new Index[n_xL_skip_+1];
    for (Index i=0; i<n_xL_skip_; i++) {
      index_xL_skip_[i] = index_xL_skip[i];
    }
    std::sort(index_xL_skip_, index_xL_skip_+n_xL_skip_);
    index_xL_skip_[n_xL_skip_] = -1;

    index_xU_skip_ = new Index[n_xU_skip_+1];
    for (Index i=0; i<n_xU_skip_; i++) {
      index_xU_skip_[i] = index_xU_skip[i];
    }
    std::sort(index_xU_skip_, index_xU_skip_+n_xU_skip_);
    index_xU_skip_[n_xU_skip_] = -1;

    index_x_fix_ = new Index[n_x_fix_+1];
    for (Index i=0; i<n_x_fix_; i++) {
      index_x_fix_[i] = index_x_fix[i];
    }
    std::sort(index_x_fix_, index_x_fix_+n_x_fix_);
    index_x_fix_[n_x_fix_] = -1;
  }

  TNLPReducer::~TNLPReducer()
  {
    delete [] index_g_skip_;
    delete [] g_keep_map_;
    delete [] jac_g_skipped_;
    delete [] index_xL_skip_;
    delete [] index_xU_skip_;
    delete [] index_x_fix_;
  }

  bool
  TNLPReducer::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style)
  {
    bool retval = tnlp_->get_nlp_info(n, m_orig_, nnz_jac_g_orig_,
                                      nnz_h_lag, index_style_orig_);
    if (!retval) return false;

    // If we haven't computed the conversion map yet, let's to this now
    if (m_reduced_ == -1) {
      if (index_style_orig_ == FORTRAN_STYLE) {
        for (Index i=0; i<n_g_skip_; i++) {
          index_g_skip_[i]--;
        }
        for (Index i=0; i<n_xL_skip_; i++) {
          index_xL_skip_[i]--;
        }
        for (Index i=0; i<n_xU_skip_; i++) {
          index_xU_skip_[i]--;
        }
        for (Index i=0; i<n_x_fix_; i++) {
          index_x_fix_[i]--;
        }
      }

      g_keep_map_ = new Index[m_orig_];
      m_reduced_ = 0;
      Index count = 0;
      for (Index i=0; i<m_orig_; i++) {
        if (index_g_skip_[count] == i) {
          g_keep_map_[i] = -1;
          count++;
        }
        else {
          g_keep_map_[i] = m_reduced_;
          m_reduced_++;
        }
      }

      // We also need to count the number of nonzeros
      Index* iRow = new Index[nnz_jac_g_orig_];
      Index* jCol = new Index[nnz_jac_g_orig_];
      retval = tnlp_->eval_jac_g(n, NULL, false, m_orig_, nnz_jac_g_orig_,
                                 iRow, jCol, NULL);
      if (!retval) {
        delete [] iRow;
        delete [] jCol;
      }
      nnz_jac_g_reduced_ = 0;
      nnz_jac_g_skipped_ = 0;
      for (Index i=0; i<nnz_jac_g_orig_; i++) {
        if (g_keep_map_[iRow[i]] != -1) {
          nnz_jac_g_reduced_++;
        }
        else {
          nnz_jac_g_skipped_++;
        }
      }

      DBG_ASSERT(nnz_jac_g_reduced_ + nnz_jac_g_skipped_ == nnz_jac_g_orig_);

      delete [] iRow;
      delete [] jCol;
    }

    m = m_reduced_;
    nnz_jac_g = nnz_jac_g_reduced_;
    index_style = index_style_orig_;

    return true;
  }

  bool
  TNLPReducer::get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u)
  {
    Number* g_l_orig = new Number[m_orig_];
    Number* g_u_orig = new Number[m_orig_];

    bool retval = tnlp_->get_bounds_info(n, x_l, x_u,
                                         m_orig_, g_l_orig, g_u_orig);

    if (retval) {

      if (n_x_fix_>0 || n_xL_skip_>0 || n_xU_skip_>0) {
        const Number huge = std::numeric_limits<Number>::max();

        // If there are fixed variables, we need to compute the starting
        // point to know what to fix the variable to
        Number* x_start = NULL;
        if (n_x_fix_>0) {
          x_start = new Number[n];
          retval = tnlp_->get_starting_point(n, true, x_start, false, NULL,
                                             NULL, m_orig_, false, NULL);
          if (!retval) return false;
        }

        Index count_xL = 0;
        Index count_xU = 0;
        Index count_x = 0;
        for (Index i=0; i<n; i++) {
          if (index_xL_skip_[count_xL]==i) {
            x_l[i] = -huge;
            count_xL++;
          }
          if (index_xU_skip_[count_xU]==i) {
            x_u[i] = huge;
            count_xU++;
          }
          if (index_x_fix_[count_x]==i) {
            x_l[i] = x_start[i];
            x_u[i] = x_start[i];
            count_x++;
          }
        }
        delete [] x_start;
      }

      for (Index i=0; i<m_orig_; i++) {
        Index& new_index = g_keep_map_[i];
        if (new_index >= 0) {
          g_l[new_index] = g_l_orig[i];
          g_u[new_index] = g_u_orig[i];
        }
      }
    }

    delete [] g_l_orig;
    delete [] g_u_orig;

    return retval;
  }


  bool
  TNLPReducer::get_scaling_parameters(Number& obj_scaling,
                                      bool& use_x_scaling, Index n,
                                      Number* x_scaling,
                                      bool& use_g_scaling, Index m,
                                      Number* g_scaling)
  {
    Number* g_scaling_orig = new Number[m_orig_];
    bool retval =
      tnlp_->get_scaling_parameters(obj_scaling, use_x_scaling, n, x_scaling,
                                    use_g_scaling, m_orig_, g_scaling_orig);

    if (retval && use_g_scaling) {
      for (Index i=0; i<m_orig_; i++) {
        Index& new_index = g_keep_map_[i];
        if (new_index >= 0) {
          g_scaling[new_index] = g_scaling_orig[i];
        }
      }
    }

    delete [] g_scaling_orig;
    return retval;
  }

  bool
  TNLPReducer::get_variables_linearity(Index n, LinearityType* var_types)
  {
    return tnlp_->get_variables_linearity(n, var_types);
  }

  bool
  TNLPReducer::get_constraints_linearity(Index m, LinearityType* const_types)
  {
    LinearityType* const_types_orig = new LinearityType[m_orig_];

    bool retval = tnlp_->get_constraints_linearity(m_orig_, const_types_orig);
    if (retval) {
      for (Index i=0; i<m_orig_; i++) {
        Index& new_index = g_keep_map_[i];
        if (new_index >= 0) {
          const_types[new_index] = const_types_orig[i];
        }
      }
    }

    delete [] const_types_orig;

    return retval;
  }

  bool
  TNLPReducer::get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda)
  {
    Number* lambda_orig = NULL;
    if (init_lambda) {
      lambda_orig = new Number[m_orig_];
    }

    bool retval = tnlp_->get_starting_point(n, init_x, x, init_z, z_L, z_U,
                                            m_orig_, init_lambda, lambda_orig);

    if (retval && init_lambda) {
      for (Index i=0; i<m_orig_; i++) {
        Index& new_index = g_keep_map_[i];
        if (new_index >= 0) {
          lambda[new_index] = lambda_orig[i];
        }
      }
    }

    delete [] lambda_orig;

    return retval;
  }


  bool
  TNLPReducer::get_warm_start_iterate(IteratesVector& warm_start_iterate)
  {
    return tnlp_->get_warm_start_iterate(warm_start_iterate);
  }

  bool
  TNLPReducer::eval_f(Index n, const Number* x, bool new_x,
                      Number& obj_value)
  {
    return tnlp_->eval_f(n, x, new_x, obj_value);
  }

  bool
  TNLPReducer::eval_grad_f(Index n, const Number* x, bool new_x,
                           Number* grad_f)
  {
    return tnlp_->eval_grad_f(n, x, new_x, grad_f);
  }

  bool
  TNLPReducer::eval_g(Index n, const Number* x, bool new_x,
                      Index m, Number* g)
  {
    Number* g_orig = new Number[m_orig_];

    bool retval = tnlp_->eval_g(n, x, new_x, m_orig_, g_orig);
    if (retval) {
      for (Index i=0; i<m_orig_; i++) {
        Index& new_index = g_keep_map_[i];
        if (new_index >= 0) {
          g[new_index] = g_orig[i];
        }
      }
    }

    delete [] g_orig;
    return retval;
  }

  bool
  TNLPReducer::eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow,
                          Index *jCol, Number* values)
  {
    bool retval;

    if (iRow) {
      delete [] jac_g_skipped_;
      jac_g_skipped_ = NULL;

      Index* iRow_orig = new Index[nnz_jac_g_orig_];
      Index* jCol_orig = new Index[nnz_jac_g_orig_];
      retval = tnlp_->eval_jac_g(n, x, new_x, m_orig_, nnz_jac_g_orig_,
                                 iRow_orig, jCol_orig, values);

      Index offset = 0;
      if (index_style_orig_ == FORTRAN_STYLE) {
        offset = 1;
      }
      if (retval) {
        jac_g_skipped_ = new Index[nnz_jac_g_skipped_+1];
        Index count = 0;
        Index count2 = 0;
        for (Index i=0; i<nnz_jac_g_orig_; i++) {
          Index& irow_red = g_keep_map_[iRow_orig[i]-offset];
          if (irow_red>=0) {
            iRow[count] = irow_red + offset;
            jCol[count] = jCol_orig[i];
            count++;
          }
          else {
            jac_g_skipped_[count2] = i;
            count2++;
          }
        }
        DBG_ASSERT(count == nnz_jac_g_reduced_);
        DBG_ASSERT(count2 == nnz_jac_g_skipped_);
        jac_g_skipped_[nnz_jac_g_skipped_] = -1;
      }

      delete [] iRow_orig;
      delete [] jCol_orig;
    }
    else {
      Number* values_orig = new Number[nnz_jac_g_orig_];
      retval = tnlp_->eval_jac_g(n, x, new_x, m_orig_, nnz_jac_g_orig_,
                                 iRow, jCol, values_orig);
      if (retval) {
        Index count = 0;
        Index count2 = 0;
        for (Index i=0; i<nnz_jac_g_orig_; i++) {
          if (jac_g_skipped_[count] == i) {
            count++;
          }
          else {
            values[count2] = values_orig[i];
            count2++;
          }
        }
        DBG_ASSERT(count == nnz_jac_g_skipped_);
        DBG_ASSERT(count2 == nnz_jac_g_reduced_);
      }

      delete [] values_orig;
    }

    return retval;
  }

  bool
  TNLPReducer::eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess,
                      Index* iRow, Index* jCol, Number* values)
  {
    if (!values) {
      return tnlp_->eval_h(n, x, new_x, obj_factor, m_orig_, lambda,
                           new_lambda, nele_hess, iRow, jCol, values);
    }

    Number* lambda_orig = new Number[m_orig_];
    for (Index i=0; i<m_orig_; i++) {
      Index& new_index = g_keep_map_[i];
      if (new_index >= 0) {
        lambda_orig[i] = lambda[new_index];
      }
      else {
        lambda_orig[i] = 0.;
      }
    }

    bool retval = tnlp_->eval_h(n, x, new_x, obj_factor, m_orig_, lambda_orig,
                                new_lambda, nele_hess, iRow, jCol, values);

    delete [] lambda_orig;

    return retval;
  }

  void
  TNLPReducer::finalize_solution(SolverReturn status,
                                 Index n, const Number* x,
                                 const Number* z_L, const Number* z_U,
                                 Index m, const Number* g,
                                 const Number* lambda,
                                 Number obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq)
  {
    Number* g_orig = new Number[m_orig_];
    Number* lambda_orig = new Number[m_orig_];

    // call evaluation method to get correct constraint values
    tnlp_->eval_g(n, x, true, m_orig_, g_orig);

    // fill unknown multipliers with 0
    for (Index i=0; i<m_orig_; i++) {
      Index& new_index = g_keep_map_[i];
      if (new_index >= 0) {
        lambda_orig[i] = lambda[new_index];
      }
      else {
        lambda_orig[i] = 0.;
      }
    }

    tnlp_->finalize_solution(status, n, x, z_L, z_U, m_orig_, g_orig,
                             lambda_orig, obj_value, ip_data, ip_cq);

    delete [] lambda_orig;
    delete [] g_orig;
  }

  bool
  TNLPReducer::intermediate_callback(AlgorithmMode mode,
                                     Index iter, Number obj_value,
                                     Number inf_pr, Number inf_du,
                                     Number mu, Number d_norm,
                                     Number regularization_size,
                                     Number alpha_du, Number alpha_pr,
                                     Index ls_trials,
                                     const IpoptData* ip_data,
                                     IpoptCalculatedQuantities* ip_cq)
  {
    return tnlp_->intermediate_callback(mode, iter, obj_value, inf_pr, inf_du,
                                        mu, d_norm, regularization_size,
                                        alpha_du, alpha_pr, ls_trials,
                                        ip_data, ip_cq);
  }

  Index
  TNLPReducer::get_number_of_nonlinear_variables()
  {
    return tnlp_->get_number_of_nonlinear_variables();
  }

  bool
  TNLPReducer::get_list_of_nonlinear_variables(Index num_nonlin_vars,
      Index* pos_nonlin_vars)
  {
    return tnlp_->get_list_of_nonlinear_variables(num_nonlin_vars,
           pos_nonlin_vars);
  }

}

