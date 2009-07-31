// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash     IBM    2009-06-11
//              (based on IpTNLPAdapter.cpp, rev 1465)

#include "IpParTNLPAdapter.hpp"
#include "IpBlas.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpParVector.hpp"
#include "IpParExpansionMatrix.hpp"
#include "IpParGenMatrix.hpp"
#include "IpParSymMatrix.hpp"
#include "IpOrigIpoptNLP.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#include "IpMpi.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ParTNLPAdapter::ParTNLPAdapter(const SmartPtr<ParTNLP> partnlp,
                                 const SmartPtr<const Journalist> jnlst /* = NULL */)
      :
      partnlp_(partnlp),
      jnlst_(jnlst),
      full_x_(NULL),
      part_x_(NULL),
      full_lambda_(NULL),
      part_g_(NULL),
      jac_g_part_(NULL),
      c_rhs_part_(NULL),
      x_tag_for_iterates_(0),
      y_c_tag_for_iterates_(0),
      y_d_tag_for_iterates_(0),
      x_tag_for_g_(0),
      x_tag_for_jac_g_(0),
      jac_idx_part_map_(NULL),
      h_idx_part_map_(NULL),
      x_fixed_part_map_(NULL)
  {
    ASSERT_EXCEPTION(IsValid(partnlp_), INVALID_PARTNLP,
                     "The ParTNLP passed to ParTNLPAdapter is NULL. This MUST be a valid ParTNLP!");
  }

  ParTNLPAdapter::~ParTNLPAdapter()
  {
    delete [] full_x_;
    delete [] part_x_;
    delete [] full_lambda_;
    delete [] part_g_;
    delete [] jac_g_part_;
    delete [] c_rhs_part_;
    delete [] jac_idx_part_map_;
    delete [] h_idx_part_map_;
    delete [] x_fixed_part_map_;
  }

  void ParTNLPAdapter::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("NLP");
  }

  bool ParTNLPAdapter::ProcessOptions(const OptionsList& options,
                                      const std::string& prefix)
  {
    DBG_START_METH("ParTNLPAdapter::ProcessOptions", dbg_verbosity);
    options.GetNumericValue("nlp_lower_bound_inf", nlp_lower_bound_inf_, prefix);
    options.GetNumericValue("nlp_upper_bound_inf", nlp_upper_bound_inf_, prefix);

    ASSERT_EXCEPTION(nlp_lower_bound_inf_ < nlp_upper_bound_inf_,
                     OPTION_INVALID,
                     "Option \"nlp_lower_bound_inf\" must be smaller than \"nlp_upper_bound_inf\".");

    // Registered in IpOrigIpoptNLP
    options.GetNumericValue("bound_relax_factor", bound_relax_factor_, prefix);

    Index enum_int;
    options.GetEnumValue("fixed_variable_treatment", enum_int, prefix);
    fixed_variable_treatment_ = FixedVariableTreatmentEnum(enum_int);

    // The option warm_start_same_structure is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);
    options.GetEnumValue("hessian_approximation", enum_int, prefix);
    hessian_approximation_ = HessianApproximationType(enum_int);

    return true;
  }

  bool ParTNLPAdapter::GetSpaces(SmartPtr<const VectorSpace>& x_space,
                                 SmartPtr<const VectorSpace>& c_space,
                                 SmartPtr<const VectorSpace>& d_space,
                                 SmartPtr<const VectorSpace>& x_l_space,
                                 SmartPtr<const MatrixSpace>& px_l_space,
                                 SmartPtr<const VectorSpace>& x_u_space,
                                 SmartPtr<const MatrixSpace>& px_u_space,
                                 SmartPtr<const VectorSpace>& d_l_space,
                                 SmartPtr<const MatrixSpace>& pd_l_space,
                                 SmartPtr<const VectorSpace>& d_u_space,
                                 SmartPtr<const MatrixSpace>& pd_u_space,
                                 SmartPtr<const MatrixSpace>& Jac_c_space,
                                 SmartPtr<const MatrixSpace>& Jac_d_space,
                                 SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space)
  {
    DBG_START_METH("ParTNLPAdapter::GetSpaces", dbg_verbosity);

    // Get basic MPI information
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    num_proc_ = size;
    proc_id_ = rank;

    if (warm_start_same_structure_) {
      ASSERT_EXCEPTION(full_x_, INVALID_WARMSTART,
                       "warm_start_same_structure chosen, but ParTNLPAdapter is called for the first time.");
      if (IsValid(jnlst_)) {
        jnlst_->Printf(J_DETAILED, J_INITIALIZATION,
                       "Reusing previous information for warm start in ParTNLPAdapter.\n");
      }
    }
    else {
      // In case the Adapter has been used before, but this is not a
      // warm start, make sure we delete all previously allocated
      // memory
      delete [] full_x_;
      full_x_ = NULL;
      delete [] part_x_;
      part_x_ = NULL;
      delete [] full_lambda_;
      full_lambda_ = NULL;
      delete [] part_g_;
      part_g_ = NULL;
      delete [] jac_g_part_;
      jac_g_part_ = NULL;
      delete [] c_rhs_part_;
      c_rhs_part_ = NULL;
      delete [] jac_idx_part_map_;
      jac_idx_part_map_ = NULL;
      delete [] h_idx_part_map_;
      h_idx_part_map_ = NULL;
      delete [] x_fixed_part_map_;
      x_fixed_part_map_ = NULL;
    }

    // Get the dimensions and other basic size information about the
    // problem

    Index n_full_x, n_first, n_last, n_full_g, m_first, m_last;
    Index nz_part_jac_g, nz_part_h;
    ParTNLP::IndexStyleEnum index_style;

    int retval1 =
      partnlp_->get_nlp_info(num_proc_, proc_id_, n_full_x, n_first, n_last,
                             n_full_g, m_first, m_last, nz_part_jac_g,
                             nz_part_h, index_style);
    int retval;
    MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    ASSERT_EXCEPTION(retval, INVALID_PARTNLP, "get_nlp_info returned false");
    ASSERT_EXCEPTION(!warm_start_same_structure_ ||
                     (n_full_x == n_full_x_ &&
                      n_first == n_first_ &&
                      n_last == n_last_ &&
                      n_full_g == n_full_g_ &&
                      m_first == m_first_ &&
                      m_last == m_last_ &&
                      nz_part_jac_g == nz_part_jac_g_ &&
                      nz_part_h == nz_part_h_ &&
                      index_style == index_style_),
                     INVALID_WARMSTART,
                     "warm_start_same_structure chosen, but problem dimensions are different.");
    ASSERT_EXCEPTION(n_first<=n_last+1, INVALID_PARTNLP,
                     "Condition n_first<=n_last+1 not satisfied in get_nlp_info");
    ASSERT_EXCEPTION(m_first<=m_last+1, INVALID_PARTNLP,
                     "Condition m_first<=m_last+1 not satisfied in get_nlp_info");
    if (index_style == ParTNLP::FORTRAN_STYLE) {
      n_first--;
      n_last--;
      m_first--;
      m_last--;
    }

    n_full_x_ = n_full_x;
    n_first_ = n_first;
    n_last_ = n_last;
    n_part_x_ = n_last_ - n_first_ + 1;
    n_full_g_ = n_full_g;
    m_first_ = m_first;
    m_last_ = m_last;
    n_part_g_ = m_last_ - m_first_ + 1;
    nz_part_jac_g_ = nz_part_jac_g;
    nz_part_h_ = nz_part_h;
    index_style_ = index_style;

    if (!warm_start_same_structure_) {
      // create space to store vectors that are the full length of x
      // including possibly fixed variables
      full_x_ = new Number[n_full_x_];

      // create space to store all x values that this processor is
      // responsible for, including possibly fixed variables
      part_x_ = new Number[n_part_x_];

      // create space to store vectors that area the full length of lambda
      full_lambda_ = new Number[n_full_g_];

      // create space to store vectors that are all elements of g this
      // processor is responsible for
      part_g_ = new Number[n_part_g_];

      // check if there is any meta data for the variables and constraints
      StringMetaDataMapType var_string_md;
      IntegerMetaDataMapType var_integer_md;
      NumericMetaDataMapType var_numeric_md;
      StringMetaDataMapType con_string_md;
      IntegerMetaDataMapType con_integer_md;
      NumericMetaDataMapType con_numeric_md;
      if (!partnlp_->get_var_con_metadata(num_proc_, proc_id_,
                                          n_full_x_, n_first_, n_last_,
                                          var_string_md, var_integer_md, var_numeric_md,
                                          n_full_g_, m_first_, m_last_,
                                          con_string_md, con_integer_md, con_numeric_md)) {
        var_string_md.clear();
        var_integer_md.clear();
        var_numeric_md.clear();
        con_string_md.clear();
        con_integer_md.clear();
        con_numeric_md.clear();
      }

      // allocate internal space to store the jacobian values for
      // constraints this process is responsible for
      jac_g_part_ = new Number[nz_part_jac_g_];

      /* Spaces for bounds this processor is responsible for. We need
       * to remove the fixed variables and find out which bounds do
       * not exist. */
      Number* x_l_part = new Number[n_part_x_];
      Number* x_u_part = new Number[n_part_x_];
      Number* g_l_part = new Number[n_part_g_];
      Number* g_u_part = new Number[n_part_g_];
      retval1 = partnlp_->get_bounds_info(num_proc_, proc_id_,
                                          n_full_x_, n_first_, n_last_,
                                          x_l_part, x_u_part,
                                          n_full_g_, m_first_, m_last_,
                                          g_l_part, g_u_part);
      MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
      ASSERT_EXCEPTION(retval, INVALID_PARTNLP,
                       "get_bounds_info returned false in GetSpaces");

      //*********************************************************
      // Create the spaces and permutation spaces
      //*********************************************************

      Index n_part_x_var;
      Index n_part_x_l;
      Index n_part_x_u;
      Index* x_not_fixed_part_map = new Index[n_part_x_];
      Index* x_l_part_map = new Index[n_part_x_];
      Index* x_u_part_map = new Index[n_part_x_];

      Index n_part_c;
      Index n_part_d;
      Index n_part_d_l;
      Index n_part_d_u;
      Index* c_part_map = new Index[n_part_g_]; // we do not know n_c yet!
      Index* d_part_map = new Index[n_part_g_]; // we do not know n_d yet!
      Index* d_l_part_map = new Index[n_part_g_]; // "
      Index* d_u_part_map = new Index[n_part_g_]; // "

      bool done=false;
      // We might have to do the following twice: If we detect that we
      // don't have enought degrees of freedom, we simply redo
      // everything with fixed_variable_treatment to set RELAX_BOUNDS
      Index n_part_x_fixed_max;
      while (!done) {
        n_part_x_var = 0;
        n_part_x_l = 0;
        n_part_x_u = 0;
        n_part_x_fixed_ = 0;
        Index* x_fixed_part_map_tmp = new Index[n_part_x_];

        for (Index i=0; i<n_part_x_; i++) {
          Number lower_bound = x_l_part[i];
          Number upper_bound = x_u_part[i];
          if (lower_bound == upper_bound) {
            switch (fixed_variable_treatment_) {
            case MAKE_PARAMETER:
              // Variable is fixed, remove it from the problem
              part_x_[i] = lower_bound;
              x_fixed_part_map_tmp[n_part_x_fixed_] = i;
              n_part_x_fixed_++;
              break;
            case MAKE_CONSTRAINT:
              x_fixed_part_map_tmp[n_part_x_fixed_] = i; // don't really need this
              // array then
              n_part_x_fixed_++;
              x_not_fixed_part_map[n_part_x_var] = i;
              n_part_x_var++;
              break;
            case RELAX_BOUNDS:
              x_l_part_map[n_part_x_l] = n_part_x_var;
              n_part_x_l++;
              x_u_part_map[n_part_x_u] = n_part_x_var;
              n_part_x_u++;
              n_part_x_var++;
              break;
            default:
              DBG_ASSERT(false && "invalid fixed_variable_treatment_");
            }
          }
          else if (lower_bound > upper_bound) {
            char string[128];
            snprintf(string, 127, "There are inconsistent bounds on variable %d: lower = %25.16e and upper = %25.16e.", i, lower_bound, upper_bound);
            delete [] x_l_part;
            delete [] x_u_part;
            delete [] g_l_part;
            delete [] g_u_part;
            delete [] x_not_fixed_part_map;
            delete [] x_fixed_part_map_tmp;
            delete [] x_l_part_map;
            delete [] x_u_part_map;
            delete [] c_part_map;
            delete [] d_part_map;
            delete [] d_l_part_map;
            delete [] d_u_part_map;
            THROW_EXCEPTION(INVALID_PARTNLP, string);
          }
          else {
            x_not_fixed_part_map[n_part_x_var] = i;
            if (lower_bound > nlp_lower_bound_inf_) {
              x_l_part_map[n_part_x_l] = n_part_x_var;
              n_part_x_l++;
            }

            if (upper_bound < nlp_upper_bound_inf_) {
              x_u_part_map[n_part_x_u] = n_part_x_var;
              n_part_x_u++;
            }
            n_part_x_var++;
          }
        }

        // If there are fixed variables, we keep their position around
        // for a possible warm start later or if fixed variables are
        // treated by added equality constraints
        // We need to do this for each process, even if that particular one
        // does not have a fixed variable but some other does
        MPI_Allreduce(&n_part_x_fixed_, &n_part_x_fixed_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (n_part_x_fixed_max>0) {
          delete [] x_fixed_part_map_;
          x_fixed_part_map_ = NULL;
          x_fixed_part_map_ = new Index[n_part_x_fixed_];
          for (Index i=0; i<n_part_x_fixed_; i++) {
            x_fixed_part_map_[i] = x_fixed_part_map_tmp[i];
          }
          delete [] x_fixed_part_map_tmp;
        }
        else {
          x_fixed_part_map_ = x_fixed_part_map_tmp;
        }

        // Create the spaces for c and d
        // - includes the internal permutation matrices for
        //  full_g to c and d
        // - includes the permutation matrices for d_l and d_u
        // c(x) = (P_c)T * g(x)
        // d(x) = (P_d)T * g(x)
        // d_L = (P_d_L)T * (P_d)T * g_l
        // d_U = (P_d_U)T * (P_d)T * g_u
        n_part_c = 0;
        n_part_d = 0;
        n_part_d_l = 0;
        n_part_d_u = 0;

        for (Index i=0; i<n_part_g_; i++) {
          Number lower_bound = g_l_part[i];
          Number upper_bound = g_u_part[i];
          if (lower_bound == upper_bound) {
            // equality constraint
            c_part_map[n_part_c] = i;
            n_part_c++;
          }
          else if (lower_bound > upper_bound) {
            delete [] x_l_part;
            delete [] x_u_part;
            delete [] g_l_part;
            delete [] g_u_part;
            delete [] x_not_fixed_part_map;
            delete [] x_l_part_map;
            delete [] x_u_part_map;
            delete [] c_part_map;
            delete [] d_part_map;
            delete [] d_l_part_map;
            delete [] d_u_part_map;
            char string[128];
            snprintf(string, 127, "There are inconsistent bounds on constraint %d: lower = %25.16e and upper = %25.16e.", i, lower_bound, upper_bound);
            THROW_EXCEPTION(INVALID_PARTNLP, string);
          }
          else {
            // inequality constraint
            d_part_map[n_part_d] = i;
            if (lower_bound > nlp_lower_bound_inf_) {
              d_l_part_map[n_part_d_l] = n_part_d;
              n_part_d_l++;
            }
            if (upper_bound < nlp_upper_bound_inf_) {
              d_u_part_map[n_part_d_u] = n_part_d;
              n_part_d_u++;
            }
            n_part_d++;
          }
        }

        if (fixed_variable_treatment_ == RELAX_BOUNDS ||
            n_part_x_fixed_ == 0 || n_part_x_var >= n_part_c) {
          done = true;
        }
        else {
          fixed_variable_treatment_ = RELAX_BOUNDS;
          jnlst_->Printf(J_WARNING, J_INITIALIZATION,
                         "Too few degrees of freedom (n_x = %d, n_c = %d).\n  Trying fixed_variable_treatment = RELAX_BOUNDS\n\n", n_part_x_var, n_part_c);
          THROW_EXCEPTION(OPTION_INVALID,
                          "This is not implemented yet for parallel");
        }
      } // while (!done)

      Index n_part_c_no_fixed = n_part_c;
      if (n_part_x_fixed_>=0 && fixed_variable_treatment_==MAKE_CONSTRAINT) {
        n_part_c += n_part_x_fixed_;
      }

      // Given the partial counts of free variables etc, find out the
      // global count
      int counts[8];
      counts[0] = n_part_x_var;
      counts[1] = n_part_x_l;
      counts[2] = n_part_x_u;
      counts[3] = n_part_x_fixed_;
      counts[4] = n_part_c;
      counts[5] = n_part_d;
      counts[6] = n_part_d_l;
      counts[7] = n_part_d_u;
      int outcounts[8];
      MPI_Allreduce(counts, outcounts, 8, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      Index n_x_var = outcounts[0];
      Index n_x_l = outcounts[1];
      Index n_x_u = outcounts[2];
      Index n_x_fixed = outcounts[3];
      Index n_c = outcounts[4];
      Index n_d = outcounts[5];
      Index n_d_l = outcounts[6];
      Index n_d_u = outcounts[7];

      // Compute the start positions
      MPI_Scan(counts, outcounts, 8, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      Index sp_x_var = outcounts[0] - counts[0];
      Index sp_x_l = outcounts[1] - counts[1];
      Index sp_x_u = outcounts[2] - counts[2];
      Index sp_c = outcounts[4] - counts[4];
      Index sp_d = outcounts[5] - counts[5];
      Index sp_d_l = outcounts[6] - counts[6];
      Index sp_d_u = outcounts[7] - counts[7];

      if (n_x_var == 0) {
        THROW_EXCEPTION(INVALID_PARTNLP, "All variables are fixed.  Special handling of this case not yet implemented for parallel version.");
      }

      delete [] x_l_part;
      x_l_part = NULL;
      delete [] x_u_part;
      x_u_part = NULL;
      delete [] g_l_part;
      g_l_part = NULL;
      delete [] g_u_part;
      g_u_part = NULL;

      // create x spaces
      SmartPtr<ParVectorSpace> pv_x_space =
        new ParVectorSpace(n_x_var, sp_x_var, n_part_x_var);
      x_space_ = GetRawPtr(pv_x_space);
      SmartPtr<DenseVectorSpace> dv_x_space =
        pv_x_space->LocalSpace();
      SmartPtr<ParVectorSpace> pv_x_l_space =
        new ParVectorSpace(n_x_l, sp_x_l, n_part_x_l);
      x_l_space_ = GetRawPtr(pv_x_l_space);
      SmartPtr<DenseVectorSpace> dv_x_l_space =
        pv_x_l_space->LocalSpace();
      SmartPtr<ParVectorSpace> pv_x_u_space =
        new ParVectorSpace(n_x_u, sp_x_u, n_part_x_u);
      x_u_space_ = GetRawPtr(pv_x_u_space);
      SmartPtr<DenseVectorSpace> dv_x_u_space =
        pv_x_u_space->LocalSpace();

      // This one is only for convenience for MPI communication
      pv_full_x_space_ =
        new ParVectorSpace(n_full_x_, n_first_, n_part_x_);

      if (n_part_x_fixed_max>0 && fixed_variable_treatment_==MAKE_PARAMETER) {
        P_x_part_x_space_ =
          new ParExpansionMatrixSpace(ConstPtr(pv_full_x_space_),
                                      ConstPtr(pv_x_space),
                                      x_not_fixed_part_map);
        P_x_part_x_ = P_x_part_x_space_->MakeNewParExpansionMatrix();
      }
      else {
        P_x_part_x_space_ = NULL;
        P_x_part_x_ = NULL;
      }

      P_x_x_L_space_ = new ParExpansionMatrixSpace(ConstPtr(pv_x_space),
                       ConstPtr(pv_x_l_space),
                       x_l_part_map);
      px_l_space_ = GetRawPtr(P_x_x_L_space_);
      P_x_x_L_ = P_x_x_L_space_->MakeNewParExpansionMatrix();
      P_x_x_U_space_ = new ParExpansionMatrixSpace(ConstPtr(pv_x_space),
                       ConstPtr(pv_x_u_space),
                       x_u_part_map);
      px_u_space_ = GetRawPtr(P_x_x_U_space_);
      P_x_x_U_ = P_x_x_U_space_->MakeNewParExpansionMatrix();

      // setup the variable meta data if present
      if (var_string_md.size() > 0) {
        StringMetaDataMapType::iterator iter;
        for (iter=var_string_md.begin(); iter != var_string_md.end(); iter++) {
          std::vector<std::string> string_md(n_part_x_var);
          const Index* pos_idx = NULL;
          if (IsValid(P_x_part_x_space_)) {
            pos_idx = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
            for (Index i=0; i<n_part_x_var; i++) {
              string_md[i] = iter->second[pos_idx[i]];
            }
          }
          else {
            for (Index i=0; i<n_part_x_var; i++) {
              string_md[i] = iter->second[i];
            }
          }
          dv_x_space->SetStringMetaData(iter->first, string_md);

          string_md.clear();
          string_md.resize(n_part_x_l);
          pos_idx = P_x_x_L_space_->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_x_l; i++) {
            string_md[i] = iter->second[pos_idx[i]];
          }
          dv_x_l_space->SetStringMetaData(iter->first, string_md);

          string_md.clear();
          string_md.resize(n_part_x_u);
          pos_idx = P_x_x_U_space_->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_x_u; i++) {
            string_md[i] = iter->second[pos_idx[i]];
          }
          dv_x_u_space->SetStringMetaData(iter->first, string_md);
        }
      }

      if (var_integer_md.size() > 0) {
        IntegerMetaDataMapType::iterator iter;
        for (iter=var_integer_md.begin(); iter != var_integer_md.end(); iter++) {
          std::vector<Index> integer_md(n_part_x_var);
          const Index* pos_idx = NULL;
          if (IsValid(P_x_part_x_space_)) {
            pos_idx = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
            for (Index i=0; i<n_part_x_var; i++) {
              integer_md[i] = iter->second[pos_idx[i]];
            }
          }
          else {
            for (Index i=0; i<n_part_x_var; i++) {
              integer_md[i] = iter->second[i];
            }
          }
          dv_x_space->SetIntegerMetaData(iter->first, integer_md);

          integer_md.clear();
          integer_md.resize(n_part_x_l);
          pos_idx = P_x_x_L_space_->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_x_l; i++) {
            integer_md[i] = iter->second[pos_idx[i]];
          }
          dv_x_l_space->SetIntegerMetaData(iter->first, integer_md);

          integer_md.clear();
          integer_md.resize(n_part_x_u);
          pos_idx = P_x_x_U_space_->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_x_u; i++) {
            integer_md[i] = iter->second[pos_idx[i]];
          }
          dv_x_u_space->SetIntegerMetaData(iter->first, integer_md);
        }
      }

      if (var_numeric_md.size() > 0) {
        NumericMetaDataMapType::iterator iter;
        for (iter=var_numeric_md.begin(); iter != var_numeric_md.end(); iter++) {
          std::vector<Number> numeric_md(n_part_x_var);
          const Index* pos_idx = NULL;
          if (IsValid(P_x_part_x_space_)) {
            pos_idx = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
            for (Index i=0; i<n_part_x_var; i++) {
              numeric_md[i] = iter->second[pos_idx[i]];
            }
          }
          else {
            for (Index i=0; i<n_part_x_var; i++) {
              numeric_md[i] = iter->second[i];
            }
          }
          dv_x_space->SetNumericMetaData(iter->first, numeric_md);

          numeric_md.clear();
          numeric_md.resize(n_part_x_l);
          pos_idx = P_x_x_L_space_->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_x_l; i++) {
            numeric_md[i] = iter->second[pos_idx[i]];
          }
          dv_x_l_space->SetNumericMetaData(iter->first, numeric_md);

          numeric_md.clear();
          numeric_md.resize(n_part_x_u);
          pos_idx = P_x_x_U_space_->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_x_u; i++) {
            numeric_md[i] = iter->second[pos_idx[i]];
          }
          dv_x_u_space->SetNumericMetaData(iter->first, numeric_md);
        }
      }

      delete [] x_not_fixed_part_map;
      x_not_fixed_part_map = NULL;
      delete [] x_l_part_map;
      x_l_part_map = NULL;
      delete [] x_u_part_map;
      x_u_part_map = NULL;

      // create the required c_space

      // This one is only for convenience for MPI communication
      pv_full_g_space_ =
        new ParVectorSpace(n_full_g_, m_first_, n_part_g_);

      SmartPtr<ParVectorSpace> pv_c_space =
        new ParVectorSpace(n_c, sp_c, n_part_c);
      c_space_ = GetRawPtr(pv_c_space);
      c_rhs_part_ = new Number[n_part_c];
      SmartPtr<DenseVectorSpace> dv_c_space =
        pv_c_space->LocalSpace();
      // create the internal expansion matrix for c to g
      P_c_g_space_ =
        new ExpansionMatrixSpace(n_part_g_, n_part_c_no_fixed, c_part_map);
      P_c_g_ = P_c_g_space_->MakeNewExpansionMatrix();
      delete [] c_part_map;
      c_part_map = NULL;

      // create the required d_space
      SmartPtr<ParVectorSpace> pv_d_space =
        new ParVectorSpace(n_d, sp_d, n_part_d);
      d_space_ = GetRawPtr(pv_d_space);
      SmartPtr<DenseVectorSpace> dv_d_space =
        pv_d_space->LocalSpace();
      // create the internal expansion matrix for d to g
      P_d_g_space_ = new ExpansionMatrixSpace(n_full_g_, n_part_d, d_part_map);
      P_d_g_ = P_d_g_space_->MakeNewExpansionMatrix();
      delete [] d_part_map;
      d_part_map = NULL;

      // create the required d_l space
      SmartPtr<ParVectorSpace> pv_d_l_space =
        new ParVectorSpace(n_d_l, sp_d_l, n_part_d_l);
      d_l_space_ = GetRawPtr(pv_d_l_space);
      SmartPtr<DenseVectorSpace> dv_d_l_space =
        pv_d_l_space->LocalSpace();
      // create the required expansion matrix for d_L to d_L_exp
      SmartPtr<ParExpansionMatrixSpace> P_d_l_space =
        new ParExpansionMatrixSpace(ConstPtr(pv_d_space),
                                    ConstPtr(pv_d_l_space), d_l_part_map);
      pd_l_space_ = GetRawPtr(P_d_l_space);
      delete [] d_l_part_map;
      d_l_part_map = NULL;

      // create the required d_u space
      SmartPtr<ParVectorSpace> pv_d_u_space =
        new ParVectorSpace(n_d_u, sp_d_u, n_part_d_u);
      d_u_space_ = GetRawPtr(pv_d_u_space);
      SmartPtr<DenseVectorSpace> dv_d_u_space =
        pv_d_u_space->LocalSpace();
      // create the required expansion matrix for d_u to d_U_exp
      SmartPtr<ParExpansionMatrixSpace> P_d_u_space =
        new ParExpansionMatrixSpace(ConstPtr(pv_d_space),
                                    ConstPtr(pv_d_u_space), d_u_part_map);
      pd_u_space_ = GetRawPtr(P_d_u_space);
      delete [] d_u_part_map;
      d_u_part_map = NULL;

      // set the constraint meta data if present
      if (con_string_md.size() > 0) {
        StringMetaDataMapType::iterator iter;
        for (iter=con_string_md.begin(); iter != con_string_md.end(); iter++) {
          std::vector<std::string> string_md(n_part_c_no_fixed);
          const Index* pos_idx = P_c_g_space_->ExpandedPosIndices();
          for (Index i=0; i<n_part_c_no_fixed; i++) {
            string_md[i] = iter->second[pos_idx[i]];
          }
          dv_c_space->SetStringMetaData(iter->first, string_md);

          string_md.clear();
          string_md.resize(n_part_d);
          pos_idx = P_d_g_space_->ExpandedPosIndices();
          for (Index i=0; i<n_part_d; i++) {
            string_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_space->SetStringMetaData(iter->first, string_md);

          string_md.clear();
          string_md.resize(n_part_d_l);
          pos_idx = P_d_l_space->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_d_l; i++) {
            string_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_l_space->SetStringMetaData(iter->first, string_md);

          string_md.clear();
          string_md.resize(n_part_d_u);
          pos_idx = P_d_u_space->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_d_u; i++) {
            string_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_u_space->SetStringMetaData(iter->first, string_md);
        }
      }

      if (con_integer_md.size() > 0) {
        IntegerMetaDataMapType::iterator iter;
        for (iter=con_integer_md.begin(); iter != con_integer_md.end(); iter++) {
          std::vector<Index> integer_md(n_part_c_no_fixed);
          const Index* pos_idx = P_c_g_space_->ExpandedPosIndices();
          for (Index i=0; i<n_part_c_no_fixed; i++) {
            integer_md[i] = iter->second[pos_idx[i]];
          }
          dv_c_space->SetIntegerMetaData(iter->first, integer_md);

          integer_md.clear();
          integer_md.resize(n_part_d);
          pos_idx = P_d_g_space_->ExpandedPosIndices();
          for (Index i=0; i<n_part_d; i++) {
            integer_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_space->SetIntegerMetaData(iter->first, integer_md);

          integer_md.clear();
          integer_md.resize(n_part_d_l);
          pos_idx = P_d_l_space->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_d_l; i++) {
            integer_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_l_space->SetIntegerMetaData(iter->first, integer_md);

          integer_md.clear();
          integer_md.resize(n_part_d_u);
          pos_idx = P_d_u_space->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_d_u; i++) {
            integer_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_u_space->SetIntegerMetaData(iter->first, integer_md);
        }
      }

      if (con_numeric_md.size() > 0) {
        NumericMetaDataMapType::iterator iter;
        for (iter=con_numeric_md.begin(); iter != con_numeric_md.end(); iter++) {
          std::vector<Number> numeric_md(n_part_c_no_fixed);
          const Index* pos_idx = P_c_g_space_->ExpandedPosIndices();
          for (Index i=0; i<n_part_c_no_fixed; i++) {
            numeric_md[i] = iter->second[pos_idx[i]];
          }
          dv_c_space->SetNumericMetaData(iter->first, numeric_md);

          numeric_md.clear();
          numeric_md.resize(n_part_d);
          pos_idx = P_d_g_space_->ExpandedPosIndices();
          for (Index i=0; i<n_part_d; i++) {
            numeric_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_space->SetNumericMetaData(iter->first, numeric_md);

          numeric_md.clear();
          numeric_md.resize(n_part_d_l);
          pos_idx = P_d_l_space->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_d_l; i++) {
            numeric_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_l_space->SetNumericMetaData(iter->first, numeric_md);

          numeric_md.clear();
          numeric_md.resize(n_part_d_u);
          pos_idx = P_d_u_space->LocalSpace()->ExpandedPosIndices();
          for (Index i=0; i<n_part_d_u; i++) {
            numeric_md[i] = iter->second[pos_idx[i]];
          }
          dv_d_u_space->SetNumericMetaData(iter->first, numeric_md);
        }
      }

      /** Create the matrix space for the jacobians
       */
      // Get the non zero structure
      Index* g_iRow = new Index[nz_part_jac_g_];
      Index* g_jCol = new Index[nz_part_jac_g_];
      retval1 = partnlp_->eval_jac_g(num_proc_, proc_id_, n_full_x, NULL, false,
                                     n_full_g, m_first, m_last, nz_part_jac_g_,
                                     g_iRow, g_jCol, NULL);
      MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
      if (!retval) {
        delete [] g_iRow;
        delete [] g_jCol;
        THROW_EXCEPTION(INVALID_PARTNLP, "eval_jac_g returned false when asked for structure");
      }

      if (index_style_ != ParTNLP::FORTRAN_STYLE) {
        for (Index i=0; i<nz_part_jac_g_; i++) {
          g_iRow[i] += 1;
          g_jCol[i] += 1;
        }
      }
      DBG_DO(for (Index i=0; i<nz_part_jac_g_; i++) assert(g_iRow[i]>0 && g_iRow[i]<=n_part_g_ && g_jCol[i]>0 && g_jCol[i]<=n_full_x);)

      // ... build the non-zero structure for jac_c
      // ... (the permutation from rows in jac_g to jac_c is
      // ...  the same as P_c_g_)
      Index nz_jac_all;
      if (fixed_variable_treatment_==MAKE_PARAMETER) {
        nz_jac_all = nz_part_jac_g_;
      }
      else {
        nz_jac_all = nz_part_jac_g_ + n_part_x_fixed_;
      }
      jac_idx_part_map_ = new Index[nz_jac_all];
      Index* jac_c_iRow = new Index[nz_jac_all];
      Index* jac_c_jCol = new Index[nz_jac_all];
      Index current_nz = 0;
      const Index* c_row_pos = P_c_g_->CompressedPosIndices();
      if (IsValid(P_x_part_x_)) {
        // there are missing variables x
        const Index* c_col_pos = P_x_part_x_->GlobalCompressedPosIndices();
        for (Index i=0; i<nz_part_jac_g_; i++) {
          const Index& c_row = c_row_pos[g_iRow[i]-1];
          const Index& c_col = c_col_pos[g_jCol[i]-1];
          if (c_col != -1 && c_row != -1) {
            jac_idx_part_map_[current_nz] = i;
            jac_c_iRow[current_nz] = c_row + 1;
            jac_c_jCol[current_nz] = c_col + 1;
            current_nz++;
          }
        }
      }
      else {
        for (Index i=0; i<nz_part_jac_g_; i++) {
          const Index& c_row = c_row_pos[g_iRow[i]-1];
          const Index& c_col = g_jCol[i]-1;
          if (c_row != -1) {
            jac_idx_part_map_[current_nz] = i;
            jac_c_iRow[current_nz] = c_row + 1;
            jac_c_jCol[current_nz] = c_col + 1;
            current_nz++;
          }
        }
      }
      nz_part_jac_c_no_extra_ = current_nz;
      if (fixed_variable_treatment_==MAKE_PARAMETER) {
        nz_part_jac_c_ = nz_part_jac_c_no_extra_;
      }
      else {
        nz_part_jac_c_ = nz_part_jac_c_no_extra_ + n_part_x_fixed_;
        for (Index i=0; i<n_part_x_fixed_; i++) {
          jac_c_iRow[current_nz] = n_part_c_no_fixed + i + 1;
          jac_c_jCol[current_nz] = x_fixed_part_map_[i]+1;
          current_nz++;
        }
      }

      Jac_c_space_ =
        new ParGenMatrixSpace(ConstPtr(pv_c_space), n_x_var, nz_part_jac_c_,
                              jac_c_iRow, jac_c_jCol);
      delete [] jac_c_iRow;
      jac_c_iRow = NULL;
      delete [] jac_c_jCol;
      jac_c_jCol = NULL;

      // ... build the nonzero structure for jac_d
      // ... (the permuation from rows in jac_g to jac_c is the
      // ...  the same as P_d_g_)
      Index* jac_d_iRow = new Index[nz_part_jac_g_];
      Index* jac_d_jCol = new Index[nz_part_jac_g_];
      current_nz = 0;
      const Index* d_row_pos = P_d_g_->CompressedPosIndices();
      if (IsValid(P_x_part_x_)) {
        const Index* d_col_pos = P_x_part_x_->GlobalCompressedPosIndices();
        for (Index i=0; i<nz_part_jac_g_; i++) {
          const Index& d_row = d_row_pos[g_iRow[i]-1];
          const Index& d_col = d_col_pos[g_jCol[i]-1];
          if (d_col != -1 && d_row != -1) {
            jac_idx_part_map_[current_nz + nz_part_jac_c_no_extra_] = i;
            jac_d_iRow[current_nz] = d_row + 1;
            jac_d_jCol[current_nz] = d_col + 1;
            current_nz++;
          }
        }
      }
      else {
        for (Index i=0; i<nz_part_jac_g_; i++) {
          const Index& d_row = d_row_pos[g_iRow[i]-1];
          const Index& d_col = g_jCol[i]-1;
          if (d_row != -1) {
            jac_idx_part_map_[current_nz + nz_part_jac_c_no_extra_] = i;
            jac_d_iRow[current_nz] = d_row + 1;
            jac_d_jCol[current_nz] = d_col + 1;
            current_nz++;
          }
        }
      }
      nz_part_jac_d_ = current_nz;
      Jac_d_space_ =
        new ParGenMatrixSpace(ConstPtr(pv_d_space), n_x_var, nz_part_jac_d_,
                              jac_d_iRow, jac_d_jCol);
      delete [] jac_d_iRow;
      jac_d_iRow = NULL;
      delete [] jac_d_jCol;
      jac_d_jCol = NULL;

      delete [] g_iRow;
      g_iRow = NULL;
      delete [] g_jCol;
      g_jCol = NULL;

      if (hessian_approximation_==EXACT) {
        /** Create the matrix space for the hessian of the lagrangian */
        Index* full_h_iRow = new Index[nz_part_h_];
        Index* full_h_jCol = new Index[nz_part_h_];
        Index* h_iRow = new Index[nz_part_h_];
        Index* h_jCol = new Index[nz_part_h_];
        retval1 = partnlp_->eval_h(num_proc_, proc_id_,
                                   n_full_x_, n_first_, n_last_, NULL, false,
                                   0., n_full_g_, m_first_, m_last_,
                                   NULL, false,
                                   nz_part_h_, full_h_iRow, full_h_jCol, NULL);

        MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        if (!retval) {
          delete [] full_h_iRow;
          delete [] full_h_jCol;
          delete [] h_iRow;
          delete [] h_jCol;
          jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                         "Option \"hessian_approximation\" is not chosen as \"limited-memory\", but eval_h returns false.\n");
          THROW_EXCEPTION(OPTION_INVALID, "eval_h is called but has not been implemented");
        }

        if (index_style_ != ParTNLP::FORTRAN_STYLE) {
          for (Index i=0; i<nz_part_h_; i++) {
            full_h_iRow[i] += 1;
            full_h_jCol[i] += 1;
          }
        }

        current_nz = 0;
        if (IsValid(P_x_part_x_)) {
          h_idx_part_map_ = new Index[nz_part_h_];
          const Index* h_pos = P_x_part_x_->GlobalCompressedPosIndices();
          for (Index i=0; i<nz_part_h_; i++) {
            const Index& h_row = h_pos[full_h_iRow[i]-1];
            const Index& h_col = h_pos[full_h_jCol[i]-1];
            if (h_row != -1 && h_col != -1) {
              h_idx_part_map_[current_nz] = i;
              h_iRow[current_nz] = h_row + 1;
              h_jCol[current_nz] = h_col + 1;
              current_nz++;
            }
          }
        }
        else {
          h_idx_part_map_ = NULL;
          for (Index i=0; i<nz_part_h_; i++) {
            const Index& h_row = full_h_iRow[i]-1;
            const Index& h_col = full_h_jCol[i]-1;
            h_iRow[i] = h_row + 1;
            h_jCol[i] = h_col + 1;
            current_nz++;
          }
          current_nz = nz_part_h_;
        }
        nz_h_ = current_nz;
        Hess_lagrangian_space_ =
          new ParSymMatrixSpace(n_x_var, nz_h_, h_iRow, h_jCol);
        delete [] full_h_iRow;
        full_h_iRow = NULL;
        delete [] full_h_jCol;
        full_h_jCol = NULL;
        delete [] h_iRow;
        h_iRow = NULL;
        delete [] h_jCol;
        h_jCol = NULL;
      }
      else { /* if (hessian_approximation_==EXACT) {*/
        nz_h_ = 0;
        Hess_lagrangian_space_ = NULL;
      }
    } /* if (warm_start_same_structure_) { */

    // Assign the spaces to the returned pointers
    x_space = x_space_;
    c_space = c_space_;
    d_space = d_space_;
    x_l_space = x_l_space_;
    px_l_space = px_l_space_;
    x_u_space = x_u_space_;
    px_u_space = px_u_space_;
    d_l_space = d_l_space_;
    pd_l_space = pd_l_space_;
    d_u_space = d_u_space_;
    pd_u_space = pd_u_space_;
    Jac_c_space = Jac_c_space_;
    Jac_d_space = Jac_d_space_;
    Hess_lagrangian_space = Hess_lagrangian_space_;

    if (IsValid(jnlst_)) {
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS,
                     "Number of nonzeros in equality constraint Jacobian...:%9d\n", nz_part_jac_c_);
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS,
                     "Number of nonzeros in inequality constraint Jacobian.:%9d\n", nz_part_jac_d_);
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS,
                     "Number of nonzeros in Lagrangian Hessian.............:%9d\n\n", nz_h_);
    }

    return true;
  }

  bool ParTNLPAdapter::GetBoundsInformation(const Matrix& Px_L,
      Vector& x_L,
      const Matrix& Px_U,
      Vector& x_U,
      const Matrix& Pd_L,
      Vector& d_L,
      const Matrix& Pd_U,
      Vector& d_U)
  {
    // This could be done more efficiently, I have already called this method
    // once to setup the structure for the problem, I could store the values
    // and use them here ?
    // Actually, this is better for a warm start
    Number* x_l_part = new Number[n_part_x_];
    Number* x_u_part = new Number[n_part_x_];
    Number* g_l_part = new Number[n_part_g_];
    Number* g_u_part = new Number[n_part_g_];
    int retval1 = partnlp_->get_bounds_info(num_proc_, proc_id_,
                                            n_full_x_, n_first_, n_last_,
                                            x_l_part, x_u_part,
                                            n_full_g_, m_first_, m_last_,
                                            g_l_part, g_u_part);
    int retval;
    MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    ASSERT_EXCEPTION(retval, INVALID_PARTNLP,
                     "get_bounds_info returned false in GetBoundsInformation");

    if (fixed_variable_treatment_==MAKE_PARAMETER) {
      // Set the values of fixed variables
      for (Index i=0; i<n_part_x_fixed_; i++) {
        DBG_ASSERT(x_l_part[x_fixed_part_map_[i]] == x_u_part[x_fixed_part_map_[i]]);
        part_x_[x_fixed_part_map_[i]] = x_l_part[x_fixed_part_map_[i]];
      }
    }
    else if (fixed_variable_treatment_==RELAX_BOUNDS) {
      // Relax the bounds for fixed variables
      const Number bound_relax = Max(1e-8, bound_relax_factor_);
      for (Index i=0; i<n_part_x_fixed_; i++) {
        if (x_l_part[i] == x_u_part[i]) {
          x_l_part[i] -= bound_relax*Max(1.,fabs(x_l_part[i]));
          x_u_part[i] += bound_relax*Max(1.,fabs(x_u_part[i]));
        }
      }
    }

    // Set the bounds values for x
    ParVector* px_L = static_cast<ParVector*>(&x_L);
    DBG_ASSERT(dynamic_cast<ParVector*>(&x_L));
    DenseVector* dx_L = px_L->LocalVector();
    Number* values = dx_L->Values();

    const ParExpansionMatrix* pem_Px_L =
      static_cast<const ParExpansionMatrix*>(&Px_L);
    DBG_ASSERT(dynamic_cast<const ParExpansionMatrix*>(&Px_L));
    const ExpansionMatrix* em_Px_L = pem_Px_L->LocalMatrix();

    if (IsValid(P_x_part_x_)) {
      for (Index i=0; i<dx_L->Dim(); i++) {
        const Index& ipopt_idx = em_Px_L->ExpandedPosIndices()[i];
        const Index& full_idx = P_x_part_x_->LocalMatrix()->ExpandedPosIndices()[ipopt_idx];
        const Number& lower_bound = x_l_part[full_idx];
        values[i] = lower_bound;
      }
    }
    else {
      for (Index i=0; i<dx_L->Dim(); i++) {
        const Index& ipopt_idx = em_Px_L->ExpandedPosIndices()[i];
        const Number& lower_bound = x_l_part[ipopt_idx];
        values[i] = lower_bound;
      }
    }

    ParVector* px_U = static_cast<ParVector*>(&x_U);
    DBG_ASSERT(dynamic_cast<ParVector*>(&x_U));
    DenseVector* dx_U = px_U->LocalVector();
    values = dx_U->Values();

    const ParExpansionMatrix* pem_Px_U =
      static_cast<const ParExpansionMatrix*>(&Px_U);
    DBG_ASSERT(dynamic_cast<const ParExpansionMatrix*>(&Px_U));
    const ExpansionMatrix* em_Px_U = pem_Px_U->LocalMatrix();

    if (IsValid(P_x_part_x_)) {
      for (Index i=0; i<dx_U->Dim(); i++) {
        const Index& ipopt_idx = em_Px_U->ExpandedPosIndices()[i];
        const Index& full_idx = P_x_part_x_->LocalMatrix()->ExpandedPosIndices()[ipopt_idx];
        const Number& upper_bound = x_u_part[full_idx];
        values[i] = upper_bound;
      }
    }
    else {
      for (Index i=0; i<dx_U->Dim(); i++) {
        const Index& ipopt_idx = em_Px_U->ExpandedPosIndices()[i];
        const Number& upper_bound = x_u_part[ipopt_idx];
        values[i] = upper_bound;
      }
    }

    // get the bounds values (rhs values to subtract) for c
    // i.e. if gL == gU, then we actually have g(x) = gL = gU,
    // since we solve c(x) = 0, we actually need c(x) - gL = 0
    for (Index i=0; i<P_c_g_->NCols(); i++) {
      Index full_idx = P_c_g_->ExpandedPosIndices()[i];
      Number rhs = g_l_part[full_idx];
      c_rhs_part_[i] = rhs;
    }
    // similarly, if we have fixed variables, consider them here
    if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
      Index n_part_c_no_fixed = P_c_g_->NCols();
      for (Index i=0; i<n_part_x_fixed_; i++) {
        DBG_ASSERT(x_l_part[x_fixed_part_map_[i]]==x_u_part[x_fixed_part_map_[i]]);
        c_rhs_part_[n_part_c_no_fixed+i] = x_l_part[x_fixed_part_map_[i]];
      }
    }

    // get the bounds values for d
    ParVector* pd_L = static_cast<ParVector*>(&d_L);
    DBG_ASSERT(dynamic_cast<ParVector*>(&d_L));
    DenseVector* dd_L = pd_L->LocalVector();
    values = dd_L->Values();

    const ParExpansionMatrix* pem_Pd_L =
      static_cast<const ParExpansionMatrix*>(&Pd_L);
    DBG_ASSERT(dynamic_cast<const ParExpansionMatrix*>(&Pd_L));
    const ExpansionMatrix* em_Pd_L = pem_Pd_L->LocalMatrix();
    for (Index i=0; i<em_Pd_L->NCols(); i++) {
      Index d_exp_idx = em_Pd_L->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number lower_bound = g_l_part[full_idx];
      values[i] = lower_bound;
    }

    ParVector* pd_U = static_cast<ParVector*>(&d_U);
    DBG_ASSERT(dynamic_cast<ParVector*>(&d_U));
    DenseVector* dd_U = pd_U->LocalVector();
    values = dd_U->Values();

    const ParExpansionMatrix* pem_Pd_U =
      static_cast<const ParExpansionMatrix*>(&Pd_U);
    DBG_ASSERT(dynamic_cast<const ParExpansionMatrix*>(&Pd_U));
    const ExpansionMatrix* em_Pd_U = pem_Pd_U->LocalMatrix();
    for (Index i=0; i<em_Pd_U->NCols(); i++) {
      Index d_exp_idx = em_Pd_U->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number upper_bound = g_u_part[full_idx];
      values[i] = upper_bound;
    }

    delete [] x_l_part;
    x_l_part = NULL;
    delete [] x_u_part;
    x_u_part = NULL;
    delete [] g_l_part;
    g_l_part = NULL;
    delete [] g_u_part;
    g_u_part = NULL;

    return true;
  }

  bool ParTNLPAdapter::GetStartingPoint(SmartPtr<Vector> x,
                                        bool need_x,
                                        SmartPtr<Vector> y_c,
                                        bool need_y_c,
                                        SmartPtr<Vector> y_d,
                                        bool need_y_d,
                                        SmartPtr<Vector> z_L,
                                        bool need_z_L,
                                        SmartPtr<Vector> z_U,
                                        bool need_z_U)
  {
    Number* x_part = new Number[n_part_x_];
    Number* z_l_part = new Number[n_part_x_];
    Number* z_u_part = new Number[n_part_x_];
    Number* lambda_part = new Number[n_part_g_];
    bool init_x = need_x;
    bool init_z = need_z_L || need_z_U;
    bool init_lambda = need_y_c || need_y_d;

    int retval1 =
      partnlp_->get_starting_point(num_proc_, proc_id_,
                                   n_full_x_, n_first_, n_last_,
                                   init_x, x_part,
                                   init_z, z_l_part, z_u_part,
                                   n_full_g_, m_first_, m_last_,
                                   init_lambda, lambda_part);
    int retval;
    MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    if (!retval) {
      delete [] x_part;
      delete [] z_l_part;
      delete [] z_u_part;
      delete [] lambda_part;
      return false;
    }

    if (need_x) {
      ParVector* px = static_cast<ParVector*>(GetRawPtr(x));
      DBG_ASSERT(dynamic_cast<ParVector*>(GetRawPtr(x)));
      DenseVector* dx = px->LocalVector();
      Number* values = dx->Values();
      const Index& n_x_var = dx->Dim();
      if (IsValid(P_x_part_x_)) {
        const Index* x_pos = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
        for (Index i=0; i<n_x_var; i++) {
          values[i] = x_part[x_pos[i]];
        }
      }
      else {
        IpBlasDcopy(n_x_var, x_part, 1, values, 1);
      }
    }

    if (need_y_c) {
      ParVector* py_c = static_cast<ParVector*>(GetRawPtr(y_c));
      DBG_ASSERT(dynamic_cast<ParVector*>(GetRawPtr(y_c)));
      DenseVector* dy_c = py_c->LocalVector();
      Number* values = dy_c->Values();
      const Index* y_c_pos = P_c_g_->ExpandedPosIndices();
      for (Index i=0; i<P_c_g_->NCols(); i++) {
        values[i] = lambda_part[y_c_pos[i]];
      }
      if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
        // ToDo maybe use info from z_L and Z_U here?
        const Number zero = 0.;
        IpBlasDcopy(n_part_x_fixed_, &zero, 0, &values[P_c_g_->NCols()], 1);
      }
    }

    if (need_y_d) {
      ParVector* py_d = static_cast<ParVector*>(GetRawPtr(y_d));
      DBG_ASSERT(dynamic_cast<ParVector*>(GetRawPtr(y_d)));
      DenseVector* dy_d = py_d->LocalVector();
      Number* values = dy_d->Values();
      const Index* y_d_pos = P_d_g_->ExpandedPosIndices();
      for (Index i=0; i<P_d_g_->NCols(); i++) {
        values[i] = lambda_part[y_d_pos[i]];
      }
    }

    if (need_z_L) {
      ParVector* pz_l = static_cast<ParVector*>(GetRawPtr(z_L));
      DBG_ASSERT(dynamic_cast<ParVector*>(GetRawPtr(z_L)));
      DenseVector* dz_L = pz_l->LocalVector();
      Number* values = dz_L->Values();
      const Index& n_part_z_l = dz_L->Dim();
      const Index* z_l_pos = P_x_x_L_->LocalMatrix()->ExpandedPosIndices();
      if (IsValid(P_x_part_x_)) {
        const Index* x_pos = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
        for (Index i=0; i<n_part_z_l; i++) {
          Index idx = z_l_pos[i]; // convert from x_L to x (ipopt)
          idx = x_pos[idx]; // convert from x (ipopt) to x_full
          values[i] = z_l_part[idx];
        }
      }
      else {
        for (Index i=0; i<n_part_z_l; i++) {
          Index idx = z_l_pos[i]; // convert from x_L to x (ipopt)
          values[i] = z_l_part[idx];
        }
      }
    }

    if (need_z_U) {
      ParVector* pz_u = static_cast<ParVector*>(GetRawPtr(z_U));
      DBG_ASSERT(dynamic_cast<ParVector*>(GetRawPtr(z_U)));
      DenseVector* dz_U = pz_u->LocalVector();
      Number* values = dz_U->Values();
      const Index& n_part_z_u = dz_U->Dim();
      const Index* z_u_pos = P_x_x_U_->LocalMatrix()->ExpandedPosIndices();
      if (IsValid(P_x_part_x_)) {
        const Index* x_pos = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
        for (Index i=0; i<n_part_z_u; i++) {
          Index idx = z_u_pos[i]; // convert from x_u to x (ipopt)
          idx = x_pos[idx]; // convert from x (ipopt) to x_full
          values[i] = z_u_part[idx];
        }
      }
      else {
        for (Index i=0; i<n_part_z_u; i++) {
          Index idx = z_u_pos[i]; // convert from x_u to x (ipopt)
          values[i] = z_u_part[idx];
        }
      }
    }

    delete [] x_part;
    x_part = NULL;
    delete [] z_l_part;
    z_l_part = NULL;
    delete [] z_u_part;
    z_u_part = NULL;
    delete [] lambda_part;
    lambda_part = NULL;

    return true;
  }

  bool ParTNLPAdapter::Eval_f(const Vector& x, Number& f)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }
    Number f1;
    int retval1 = partnlp_->eval_f(num_proc_, proc_id_,
                                   n_full_x_, n_first_, n_last_,
                                   full_x_, new_x, f1);

    // synchonize return values (TODO: What operation?)
    int retval;
    MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    if (retval==0) return false;

    // sum up all returned objective function values
    MPI_Allreduce(&f1, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return true;
  }

  bool ParTNLPAdapter::Eval_grad_f(const Vector& x, Vector& g_f)
  {
    int retvalue;
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    ParVector* pg_f = static_cast<ParVector*>(&g_f);
    DBG_ASSERT(dynamic_cast<ParVector*>(&g_f));
    DenseVector* dg_f = pg_f->LocalVector();
    Number* values = dg_f->Values();
    if (IsValid(P_x_part_x_)) {
      Number* grad_f_part = new Number[n_part_x_];
      if (partnlp_->eval_grad_f(num_proc_, proc_id_,
                                n_full_x_, n_first_, n_last_,
                                full_x_, new_x, grad_f_part)) {
        const Index* x_pos = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
        for (Index i=0; i<dg_f->Dim(); i++) {
          values[i] = grad_f_part[x_pos[i]];
        }
        retvalue = 1;
      }
      delete [] grad_f_part;
    }
    else {
      retvalue = partnlp_->eval_grad_f(num_proc_, proc_id_,
                                       n_full_x_, n_first_, n_last_,
                                       full_x_, new_x, values);
    }

    // synchonize return values (TODO: What operation?)
    int retvalue_all;
    MPI_Allreduce(&retvalue, &retvalue_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    return retvalue_all;
  }

  bool ParTNLPAdapter::Eval_c(const Vector& x, Vector& c)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_g(new_x)) {
      ParVector* pc = static_cast<ParVector*>(&c);
      DBG_ASSERT(dynamic_cast<ParVector*>(&c));
      DenseVector* dc = pc->LocalVector();
      Number* values = dc->Values();
      const Index* c_pos = P_c_g_->ExpandedPosIndices();
      Index n_part_c_no_fixed = P_c_g_->NCols();
      for (Index i=0; i<n_part_c_no_fixed; i++) {
        values[i] = part_g_[c_pos[i]];
        values[i] -= c_rhs_part_[i];
      }
      if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
        for (Index i=0; i<n_part_x_fixed_; i++) {
          values[n_part_c_no_fixed+i] =
            part_x_[x_fixed_part_map_[i]] - c_rhs_part_[n_part_c_no_fixed+i];
        }
      }
      return true;
    }

    return false;
  }

  bool ParTNLPAdapter::Eval_jac_c(const Vector& x, Matrix& jac_c)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_jac_g(new_x)) {
      ParGenMatrix* p_jac_c = static_cast<ParGenMatrix*>(&jac_c);
      DBG_ASSERT(dynamic_cast<ParGenMatrix*>(&jac_c));
      GenTMatrix* gt_jac_c = p_jac_c->LocalMatrix();
      Number* values = gt_jac_c->Values();

      for (Index i=0; i<nz_part_jac_c_no_extra_; i++) {
        // Assume the same structure as initially given
        values[i] = jac_g_part_[jac_idx_part_map_[i]];
      }
      if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
        const Number one = 1.;
        IpBlasDcopy(n_part_x_fixed_, &one, 0, &values[nz_part_jac_c_no_extra_], 1);
      }
      return true;
    }
    return false;
  }

  bool ParTNLPAdapter::Eval_d(const Vector& x, Vector& d)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_g(new_x)) {
      ParVector* pd = static_cast<ParVector*>(&d);
      DBG_ASSERT(dynamic_cast<ParVector*>(&d));
      DenseVector* dd = pd->LocalVector();
      Number* values = dd->Values();
      const Index* d_pos = P_d_g_->ExpandedPosIndices();
      for (Index i=0; i<dd->Dim(); i++) {
        values[i] = part_g_[d_pos[i]];
      }
      return true;
    }

    return false;
  }

  bool ParTNLPAdapter::Eval_jac_d(const Vector& x, Matrix& jac_d)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_jac_g(new_x)) {
      ParGenMatrix* p_jac_d = static_cast<ParGenMatrix*>(&jac_d);
      DBG_ASSERT(dynamic_cast<ParGenMatrix*>(&jac_d));
      GenTMatrix* gt_jac_d = p_jac_d->LocalMatrix();
      Number* values = gt_jac_d->Values();

      for (Index i=0; i<nz_part_jac_d_; i++) {
        // Assume the same structure as initially given
        values[i] = jac_g_part_[jac_idx_part_map_[nz_part_jac_c_no_extra_ + i]];
      }
      return true;
    }
    return false;
  }

  bool ParTNLPAdapter::Eval_h(const Vector& x,
                              Number obj_factor,
                              const Vector& yc,
                              const Vector& yd,
                              SymMatrix& h)
  {
    ParSymMatrix* p_h = static_cast<ParSymMatrix*>(&h);
    DBG_ASSERT(dynamic_cast<ParSymMatrix*>(&h));
    SymTMatrix* st_h = p_h->LocalMatrix();
    Number* values = st_h->Values();

    // First see if all weights are set to zero (for example, when
    // computing the least square multiplier estimates, this is what
    // we do).  In that case, there is no need to compute values, just
    // set them to zero.
    if (obj_factor==0. && yc.Asum()==0. && yd.Asum()==0.) {
      for (Index i=0; i<nz_h_; i++) {
        values[i] = 0.;
      }
      return true;
    }

    int retval = 0;
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }
    bool new_y = false;
    if (update_local_lambda(yc, yd)) {
      new_y = true;
    }

    if (h_idx_part_map_) {
      Number* part_h = new Number[nz_part_h_];

      if (partnlp_->eval_h(num_proc_, proc_id_, n_full_x_, n_first_, n_last_,
                           full_x_, new_x, obj_factor,
                           n_full_g_, m_first_, m_last_, full_lambda_, new_y,
                           nz_part_h_, NULL, NULL, part_h)) {
        for (Index i=0; i<nz_h_; i++) {
          values[i] = part_h[h_idx_part_map_[i]];
        }
        retval = 1;
      }
      delete [] part_h;
    }
    else {
      retval =
        partnlp_->eval_h(num_proc_, proc_id_, n_full_x_, n_first_, n_last_,
                         full_x_, new_x, obj_factor,
                         n_full_g_, m_first_, m_last_, full_lambda_, new_y,
                         nz_part_h_, NULL, NULL, values);
    }

    int retval_all;
    MPI_Allreduce(&retval, &retval_all, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    return (bool)retval_all;
  }


  void ParTNLPAdapter::FinalizeSolution(SolverReturn status,
                                        const Vector& x, const Vector& z_L, const Vector& z_U,
                                        const Vector& c, const Vector& d,
                                        const Vector& y_c, const Vector& y_d,
                                        Number obj_value,
                                        const IpoptData* ip_data,
                                        IpoptCalculatedQuantities* ip_cq)
  {
    DBG_START_METH("ParTNLPAdapter::FinalizeSolution", dbg_verbosity);

    update_local_x(x);
    update_local_lambda(y_c, y_d);

    ResortX(x, full_x_);
    ResortG(y_c, y_d, full_lambda_);

    Number* full_g = new Number[n_full_g_];
    // TODO:
    if (c.Dim() + d.Dim() < n_full_g_) {
      const Number zero = 0.;
      IpBlasDcopy(n_full_g_, &zero, 0, full_g, 1);
    }
    ResortG(c, d, full_g);
    // To Ipopt, the equality constraints are presented with right
    // hand side zero, so we correct for the original right hand side.
    const Index* c_pos = P_c_g_->ExpandedPosIndices();
    Index n_c_no_fixed = P_c_g_->NCols();
    for (Index i=0; i<n_c_no_fixed; i++) {
      full_g[c_pos[i]] += c_rhs_part_[i];
    }

    Number* full_z_L = new Number[n_full_x_];
    Number* full_z_U = new Number[n_full_x_];
    ResortBnds(z_L, full_z_L, z_U, full_z_U);

    // Hopefully the following is correct to recover the bound
    // multipliers for fixed variables (sign ok?)
    if (fixed_variable_treatment_==MAKE_CONSTRAINT && n_part_x_fixed_>0) {
      assert(false && "Not yet implemented");
#if 0
      const DenseVector* dy_c = static_cast<const DenseVector*>(&y_c);
      DBG_ASSERT(dynamic_cast<const DenseVector*>(&y_c));
      DBG_ASSERT(!dy_c->IsHomogeneous());
      const Number* values = dy_c->Values();
      Index n_c_no_fixed = y_c.Dim() - n_x_fixed_;
      for (Index i=0; i<n_x_fixed_; i++) {
        full_z_L[x_fixed_part_map_[i]] = Max(0., -values[n_c_no_fixed+i]);
        full_z_U[x_fixed_part_map_[i]] = Max(0., values[n_c_no_fixed+i]);
      }
#endif
    }

    partnlp_->finalize_solution(status,
                                n_full_x_, full_x_, full_z_L, full_z_U,
                                n_full_g_, full_g, full_lambda_,
                                obj_value, ip_data, ip_cq);

    delete [] full_z_L;
    full_z_L = NULL;
    delete [] full_z_U;
    full_z_U = NULL;
    delete [] full_g;
    full_g = NULL;

  }

  bool ParTNLPAdapter::
  IntermediateCallBack(AlgorithmMode mode,
                       Index iter, Number obj_value,
                       Number inf_pr, Number inf_du,
                       Number mu, Number d_norm,
                       Number regularization_size,
                       Number alpha_du, Number alpha_pr,
                       Index ls_trials,
                       const IpoptData* ip_data,
                       IpoptCalculatedQuantities* ip_cq)
  {
    return
      partnlp_->intermediate_callback(mode, iter, obj_value, inf_pr, inf_du,
                                      mu, d_norm, regularization_size,
                                      alpha_du, alpha_pr, ls_trials,
                                      ip_data, ip_cq);
  }

  void ParTNLPAdapter::ResortX(const Vector& x, Number* x_orig)
  {
    const ParVector* px = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));
    const DenseVector* dx = px->LocalVector();

    Number* x_orig_part = new Number[n_part_x_];

    if (IsValid(P_x_part_x_)) {
      const Index* x_pos = P_x_part_x_->LocalMatrix()->CompressedPosIndices();

      if (dx->IsHomogeneous()) {
        const Number& scalar = dx->Scalar();
        for (Index i=0; i<n_part_x_; i++) {
          Index idx = x_pos[i];
          if (idx != -1) {
            x_orig_part[i] = scalar;
          }
          else {
            x_orig_part[i] = part_x_[i];
          }
        }
      }
      else {
        const Number* x_values = dx->Values();
        for (Index i=0; i<n_part_x_; i++) {
          Index idx = x_pos[i];
          if (idx != -1) {
            x_orig_part[i] = x_values[idx];
          }
          else {
            x_orig_part[i] = part_x_[i];
          }
        }
      }
    }
    else {
      if (dx->IsHomogeneous()) {
        const Number& scalar = dx->Scalar();
        IpBlasDcopy(n_part_x_, &scalar, 0, x_orig_part, 1);
      }
      else {
        IpBlasDcopy(n_part_x_, dx->Values(), 1, x_orig_part, 1);
      }
    }

    MPI_Allgatherv(x_orig_part, n_part_x_, MPI_DOUBLE, x_orig,
                   const_cast<int*>(&(pv_full_x_space_->RecvCounts()[0])),
                   const_cast<int*>(&(pv_full_x_space_->Displs()[0])),
                   MPI_DOUBLE, MPI_COMM_WORLD);

    delete [] x_orig_part;
  }

  void ParTNLPAdapter::ResortG(const Vector& c, const Vector& d, Number* g_orig)
  {
    Number* g_orig_part = new Number[n_part_g_];

    const ParVector* pc = static_cast<const ParVector*>(&c);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&c));
    const DenseVector* dc = pc->LocalVector();

    const Index* c_pos = P_c_g_->ExpandedPosIndices();
    if (dc->IsHomogeneous()) {
      Number scalar = dc->Scalar();
      for (Index i=0; i<P_c_g_->NCols(); i++) {
        g_orig_part[c_pos[i]] = scalar;
      }
    }
    else {
      const Number* c_values = dc->Values();
      for (Index i=0; i<P_c_g_->NCols(); i++) {
        g_orig_part[c_pos[i]] = c_values[i];
      }
    }

    const ParVector* pd = static_cast<const ParVector*>(&d);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&d));
    const DenseVector* dd = pd->LocalVector();

    const Index* d_pos = P_d_g_->ExpandedPosIndices();
    if (dd->IsHomogeneous()) {
      Number scalar = dd->Scalar();
      for (Index i=0; i<dd->Dim(); i++) {
        g_orig_part[d_pos[i]] = scalar;
      }
    }
    else {
      const Number* d_values = dd->Values();
      for (Index i=0; i<dd->Dim(); i++) {
        g_orig_part[d_pos[i]] = d_values[i];
      }
    }

    MPI_Allgatherv(g_orig_part, n_part_g_, MPI_DOUBLE, g_orig,
                   const_cast<int*>(&(pv_full_g_space_->RecvCounts()[0])),
                   const_cast<int*>(&(pv_full_g_space_->Displs()[0])),
                   MPI_DOUBLE, MPI_COMM_WORLD);

    delete [] g_orig_part;
  }

  void ParTNLPAdapter::ResortBnds(const Vector& x_L, Number* x_L_orig,
                                  const Vector& x_U, Number* x_U_orig)
  {
    if (x_L_orig) {
      Number* x_L_orig_part = new Number[n_part_x_];
      const Number zero = 0.;
      IpBlasDcopy(n_part_x_, &zero, 0, x_L_orig_part, 1);

      const ParVector* px_L = static_cast<const ParVector*>(&x_L);
      DBG_ASSERT(dynamic_cast<const ParVector*>(&x_L));
      const DenseVector* dx_L = px_L->LocalVector();

      const Index* bnds_pos_not_fixed =
        P_x_x_L_->LocalMatrix()->ExpandedPosIndices();
      const Index& n_xL = dx_L->Dim();

      if (IsValid(P_x_part_x_)) {
        const Index* bnds_pos_full = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
        if (dx_L->IsHomogeneous()) {
          Number scalar = dx_L->Scalar();
          for (Index i=0; i<n_xL; i++) {
            int idx = bnds_pos_not_fixed[i];
            int idx2 = bnds_pos_full[idx];
            x_L_orig_part[idx2] = scalar;
          }
        }
        else {
          const Number* x_L_values = dx_L->Values();
          for (Index i=0; i<n_xL; i++) {
            int idx = bnds_pos_not_fixed[i];
            int idx2 = bnds_pos_full[idx];
            x_L_orig_part[idx2] = x_L_values[i];
          }
        }
      }
      else {
        if (dx_L->IsHomogeneous()) {
          Number scalar = dx_L->Scalar();
          for (Index i=0; i<n_xL; i++) {
            int idx = bnds_pos_not_fixed[i];
            x_L_orig_part[idx] = scalar;
          }
        }
        else {
          const Number* x_L_values = dx_L->Values();
          for (Index i=0; i<n_xL; i++) {
            int idx = bnds_pos_not_fixed[i];
            x_L_orig_part[idx] = x_L_values[i];
          }
        }
      }

      MPI_Allgatherv(x_L_orig_part, n_part_x_, MPI_DOUBLE, x_L_orig,
                     const_cast<int*>(&(pv_full_x_space_->RecvCounts()[0])),
                     const_cast<int*>(&(pv_full_x_space_->Displs()[0])),
                     MPI_DOUBLE, MPI_COMM_WORLD);

      delete [] x_L_orig_part;
    }

    if (x_U_orig) {
      Number* x_U_orig_part = new Number[n_part_x_];
      const Number zero = 0.;
      IpBlasDcopy(n_part_x_, &zero, 0, x_U_orig_part, 1);

      const ParVector* px_U = static_cast<const ParVector*>(&x_U);
      DBG_ASSERT(dynamic_cast<const ParVector*>(&x_U));
      const DenseVector* dx_U = px_U->LocalVector();

      const Index* bnds_pos_not_fixed =
        P_x_x_U_->LocalMatrix()->ExpandedPosIndices();
      const Index& n_xU = dx_U->Dim();

      if (IsValid(P_x_part_x_)) {
        const Index* bnds_pos_full = P_x_part_x_->LocalMatrix()->ExpandedPosIndices();
        if (dx_U->IsHomogeneous()) {
          Number scalar = dx_U->Scalar();
          for (Index i=0; i<n_xU; i++) {
            int idx = bnds_pos_not_fixed[i];
            int idx2 = bnds_pos_full[idx];
            x_U_orig_part[idx2] = scalar;
          }
        }
        else {
          const Number* x_U_values = dx_U->Values();
          for (Index i=0; i<n_xU; i++) {
            int idx = bnds_pos_not_fixed[i];
            int idx2 = bnds_pos_full[idx];
            x_U_orig_part[idx2] = x_U_values[i];
          }
        }
      }
      else {
        if (dx_U->IsHomogeneous()) {
          Number scalar = dx_U->Scalar();
          for (Index i=0; i<n_xU; i++) {
            int idx = bnds_pos_not_fixed[i];
            x_U_orig_part[idx] = scalar;
          }
        }
        else {
          const Number* x_U_values = dx_U->Values();
          for (Index i=0; i<n_xU; i++) {
            int idx = bnds_pos_not_fixed[i];
            x_U_orig_part[idx] = x_U_values[i];
          }
        }
      }

      MPI_Allgatherv(x_U_orig_part, n_part_x_, MPI_DOUBLE, x_U_orig,
                     const_cast<int*>(&(pv_full_x_space_->RecvCounts()[0])),
                     const_cast<int*>(&(pv_full_x_space_->Displs()[0])),
                     MPI_DOUBLE, MPI_COMM_WORLD);

      delete [] x_U_orig_part;
    }

  }

  bool ParTNLPAdapter::update_local_x(const Vector& x)
  {
    if (x.GetTag() == x_tag_for_iterates_) {
      return false;
    }

    ResortX(x, full_x_);

    x_tag_for_iterates_ = x.GetTag();

    return true;
  }

  bool ParTNLPAdapter::update_local_lambda(const Vector& y_c, const Vector& y_d)
  {
    if (y_c.GetTag() == y_c_tag_for_iterates_
        && y_d.GetTag() == y_d_tag_for_iterates_) {
      return false;
    }

    ResortG(y_c, y_d, full_lambda_);

    y_c_tag_for_iterates_ = y_c.GetTag();
    y_d_tag_for_iterates_ = y_d.GetTag();

    return true;
  }

  bool ParTNLPAdapter::internal_eval_g(bool new_x)
  {
    if (x_tag_for_g_ == x_tag_for_iterates_) {
      // already calculated!
      return true;
    }

    x_tag_for_g_ = x_tag_for_iterates_;

    int retval1 = partnlp_->eval_g(num_proc_, proc_id_,
                                   n_full_x_, full_x_, new_x,
                                   n_full_g_, m_first_, m_last_, part_g_);

    // synchonize return values (TODO: What operation?)
    int retval;
    MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    if (!retval) {
      x_tag_for_jac_g_ = 0;
    }

    return (bool)retval;
  }

  bool ParTNLPAdapter::internal_eval_jac_g(bool new_x)
  {
    if (x_tag_for_jac_g_ == x_tag_for_iterates_) {
      // already calculated!
      return true;
    }

    x_tag_for_jac_g_ = x_tag_for_iterates_;

    int retval1 = partnlp_->eval_jac_g(num_proc_, proc_id_,
                                       n_full_x_, full_x_, new_x,
                                       n_full_g_, m_first_, m_last_,
                                       nz_part_jac_g_, NULL, NULL, jac_g_part_);

    // synchonize return values (TODO: What operation?)
    int retval;
    MPI_Allreduce(&retval1, &retval, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    if (!retval) {
      x_tag_for_jac_g_ = 0;
    }

    return (bool)retval;
  }

} // namespace Ipopt
