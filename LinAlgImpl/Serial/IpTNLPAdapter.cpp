// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpTNLPAdapter.hpp"

namespace Ipopt
{

  TNLPAdapter::TNLPAdapter(const SmartPtr<TNLP> tnlp)
      :
      tnlp_(tnlp),
      nlp_lower_bound_inf_(-1e19),
      nlp_upper_bound_inf_(1e19),
      n_full_x_(-1),
      n_full_g_(-1),
      nz_jac_c_(-1),
      nz_jac_d_(-1),
      nz_full_jac_g_(-1),
      nz_full_h_(-1),
      nz_h_(-1),
      full_x_(NULL),
      full_lambda_(NULL),
      full_g_(NULL),
      jac_g_(NULL),
      c_rhs_(NULL),
      jac_idx_map_(NULL),
      h_idx_map_(NULL),
      x_tag_for_iterates_(0),
      y_c_tag_for_iterates_(0),
      y_d_tag_for_iterates_(0),
      x_tag_for_g_(0),
      x_tag_for_jac_g_(0)
  {
    ASSERT_EXCEPTION(IsValid(tnlp_), INVALID_TNLP,
                     "The TNLP passed to TNLPAdapter is NULL. This MUST be a valid TNLP!");
  }

  TNLPAdapter::~TNLPAdapter()
  {
    delete [] full_x_;
    full_x_ = NULL;
    delete [] full_lambda_;
    full_lambda_ = NULL;
    delete [] full_g_;
    full_g_ = NULL;
    delete [] jac_g_;
    jac_g_ = NULL;
    delete [] c_rhs_;
    c_rhs_ = NULL;
    delete [] jac_idx_map_;
    jac_idx_map_ = NULL;
    delete [] h_idx_map_;
    h_idx_map_ = NULL;
  }

  bool TNLPAdapter::ProcessOptions(const OptionsList& options,
                                   const std::string& prefix)
  {
    Number value = 0.0;

    // Check for the algorithm options
    if (options.GetNumericValue("nlp_lower_bound_inf", value, prefix)) {
      //      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
      //		       "Option \"nlp_lower_bound_inf\": This value must be larger than 0.");
      nlp_lower_bound_inf_ = value;
    }
    else {
      nlp_lower_bound_inf_ = -1e19;
    }

    if (options.GetNumericValue("nlp_upper_bound_inf", value, prefix)) {
      //      ASSERT_EXCEPTION(value > nlp_lower_bound_inf_, OptionsList::OPTION_OUT_OF_RANGE,
      //                       "Option \"theta_max_fact\": This value must be larger than nlp_lower_bound_inf_.");
      nlp_upper_bound_inf_ = value;
    }
    else {
      nlp_upper_bound_inf_ = 1e19;
    }

    // allow TNLP to process some options
    //tnlp_->ProcessOptions(....);

    return true;
  }

  bool TNLPAdapter::GetSpaces(SmartPtr<VectorSpace>& x_space,
                              SmartPtr<VectorSpace>& c_space,
                              SmartPtr<VectorSpace>& d_space,
                              SmartPtr<VectorSpace>& x_l_space,
                              SmartPtr<MatrixSpace>& px_l_space,
                              SmartPtr<VectorSpace>& x_u_space,
                              SmartPtr<MatrixSpace>& px_u_space,
                              SmartPtr<VectorSpace>& d_l_space,
                              SmartPtr<MatrixSpace>& pd_l_space,
                              SmartPtr<VectorSpace>& d_u_space,
                              SmartPtr<MatrixSpace>& pd_u_space,
                              SmartPtr<MatrixSpace>& Jac_c_space,
                              SmartPtr<MatrixSpace>& Jac_d_space,
                              SmartPtr<SymMatrixSpace>& Hess_lagrangian_space)
  {
    // Get the full dimensions of the problem
    tnlp_->get_nlp_info(n_full_x_, n_full_g_, nz_full_jac_g_, nz_full_h_);

    // create space to store vectors that are the full length of x
    full_x_ = new Number[n_full_x_];
    //    full_grad_x_ = new Number[n_full_x_];

    // create space to store vectors that area the full length of lambda
    full_lambda_ = new Number[n_full_g_];

    // create space to store vectors that are the full length of g
    full_g_ = new Number[n_full_g_];

    // allocate internal space to store the full jacobian
    jac_g_ = new Number[nz_full_jac_g_];

    // allocate internal space to store the full hessian
    //    h_lag_ = new Number[nz_full_h];
    //    Index* full_h_lag_iRow = new Number[nz_full_h];
    //    Index* full_h_lag_jCol = new Number[nz_full_h];


    //*********************************************************
    // Create the spaces and permutation spaces
    //*********************************************************
    /* Spaces for x, x_L, and x_U. We need to remove the fixed variables
     * and find out which bounds do not exist. */
    Number* x_l = new Number[n_full_x_];
    Number* x_u = new Number[n_full_x_];
    Number* g_l = new Number[n_full_g_];
    Number* g_u = new Number[n_full_g_];
    tnlp_->get_bounds_info(n_full_x_, x_l, x_u, n_full_g_, g_l, g_u);

    Index n_x_not_fixed = 0;
    Index n_x_fixed = 0;
    Index n_x_l = 0;
    Index n_x_u = 0;
    Index* x_not_fixed_map = new Index[n_full_x_];
    Index* x_l_map = new Index[n_full_x_];
    Index* x_u_map = new Index[n_full_x_];

    for (Index i=0; i<n_full_x_; i++) {
      Number lower_bound = x_l[i];
      Number upper_bound = x_u[i];
      if (lower_bound == upper_bound) {
        // Variable is fixed, remove it from the problem
        //	DBG_ASSERT(false && "ToDo: write code to handle fixed variables. For now, remove the fixed variable from the problem or ensure an interior.\n");
        n_x_fixed++;
        full_x_[i] = lower_bound;
      }
      else if (lower_bound > upper_bound) {
	char string[128];
	sprintf(string, "There are inconsistent bounds on variable %d: lower = %25.16e and upper = %25.16e.", i, lower_bound, upper_bound);
	THROW_EXCEPTION(INVALID_TNLP, string);
      }
      else {
        x_not_fixed_map[n_x_not_fixed] = i;
        if (lower_bound > nlp_lower_bound_inf_) {
          x_l_map[n_x_l] = n_x_not_fixed;
          n_x_l++;
        }

        if (upper_bound < nlp_upper_bound_inf_) {
          x_u_map[n_x_u] = n_x_not_fixed;
          n_x_u++;
        }
        n_x_not_fixed++;
      }
    }


    x_space = new DenseVectorSpace(n_x_not_fixed);
    x_l_space = new DenseVectorSpace(n_x_l);
    x_u_space = new DenseVectorSpace(n_x_u);

    P_x_full_x_space_ = new ExpansionMatrixSpace(n_full_x_, n_x_not_fixed, x_not_fixed_map);
    P_x_full_x_ = P_x_full_x_space_->MakeNewExpansionMatrix();

    P_x_x_L_space_ = new ExpansionMatrixSpace(n_x_not_fixed, n_x_l, x_l_map);
    px_l_space = GetRawPtr(P_x_x_L_space_);
    P_x_x_L_ = P_x_x_L_space_->MakeNewExpansionMatrix();
    P_x_x_U_space_ = new ExpansionMatrixSpace(n_x_not_fixed, n_x_u, x_u_map);
    px_u_space = GetRawPtr(P_x_x_U_space_);
    P_x_x_U_ = P_x_x_U_space_->MakeNewExpansionMatrix();

    delete [] x_l;
    x_l = NULL;
    delete [] x_u;
    x_u = NULL;
    delete [] x_not_fixed_map;
    x_not_fixed_map = NULL;
    delete [] x_l_map;
    x_l_map = NULL;
    delete [] x_u_map;
    x_u_map = NULL;

    // Create the spaces for c and d
    // - includes the internal permutation matrices for
    //  full_g to c and d
    // - includes the permutation matrices for d_l and d_u
    // c(x) = (P_c)T * g(x)
    // d(x) = (P_d)T * g(x)
    // d_L = (P_d_L)T * (P_d)T * g_l
    // d_U = (P_d_U)T * (P_d)T * g_u
    Index n_c = 0;
    Index n_d = 0;
    Index n_d_l = 0;
    Index n_d_u = 0;

    Index* c_map = new Index[n_full_g_]; // we do not know n_c yet!
    Index* d_map = new Index[n_full_g_]; // we do not know n_d yet!
    Index* d_l_map = new Index[n_full_g_]; // "
    Index* d_u_map = new Index[n_full_g_]; // "
    for (Index i=0; i<n_full_g_; i++) {
      Number lower_bound = g_l[i];
      Number upper_bound = g_u[i];
      if (lower_bound == upper_bound) {
        // equality constraint
        c_map[n_c] = i;
        n_c++;
      }
      else {
        //	printf("pre_d_U[%d] = %lf\n",i,upper_bound);
        // inequality constraint
        //printf("%lf != %lf\n", lower_bound, upper_bound);
        d_map[n_d] = i;
        if (lower_bound > nlp_lower_bound_inf_) {
          d_l_map[n_d_l] = n_d;
          n_d_l++;
        }
        if (upper_bound < nlp_upper_bound_inf_) {
          //	  printf("Setting the upper bound\n");
          d_u_map[n_d_u] = n_d;
          n_d_u++;
        }
        n_d++;
      }
    }

    // create the required c_space
    SmartPtr<DenseVectorSpace> dc_space = new DenseVectorSpace(n_c);
    c_rhs_ = new Number[n_c];
    c_space = GetRawPtr(dc_space);
    // create the internal expansion matrix for c to g
    P_c_g_space_ = new ExpansionMatrixSpace(n_full_g_, n_c, c_map);
    P_c_g_ = P_c_g_space_->MakeNewExpansionMatrix();
    delete [] c_map;
    c_map = NULL;

    // create the required d_space
    d_space = new DenseVectorSpace(n_d);
    // create the internal expansion matrix for d to g
    P_d_g_space_ = new ExpansionMatrixSpace(n_full_g_, n_d, d_map);
    P_d_g_ = P_d_g_space_->MakeNewExpansionMatrix();
    delete [] d_map;
    d_map = NULL;

    // create the required d_l space
    d_l_space = new DenseVectorSpace(n_d_l);
    // create the required expansion matrix for d_L to d_L_exp
    pd_l_space = new ExpansionMatrixSpace(n_d, n_d_l, d_l_map);
    delete [] d_l_map;
    d_l_map = NULL;

    // create the required d_u space
    d_u_space = new DenseVectorSpace(n_d_u);
    // create the required expansion matrix for d_U to d_U_exp
    pd_u_space = new ExpansionMatrixSpace(n_d, n_d_u, d_u_map);
    delete [] d_u_map;
    d_u_map = NULL;

    delete [] g_l;
    g_l = NULL;
    delete [] g_u;
    g_u = NULL;

    /** Create the matrix space for the jacobians
     */
    // Get the non zero structure
    Index* g_iRow = new Index[nz_full_jac_g_];
    Index* g_jCol = new Index[nz_full_jac_g_];
    tnlp_->eval_jac_g(n_full_x_, NULL, false, n_full_g_, nz_full_jac_g_,
                      g_iRow, g_jCol, NULL);

    const Index* c_col_pos = P_x_full_x_->CompressedPosIndices();
    const Index* c_row_pos = P_c_g_->CompressedPosIndices();

    // ... build the non-zero structure for jac_c
    // ... (the permutation from rows in jac_g to jac_c is
    // ...  the same as P_c_g_)
    jac_idx_map_ = new Index[nz_full_jac_g_];
    Index* jac_c_iRow = new Index[nz_full_jac_g_];
    Index* jac_c_jCol = new Index[nz_full_jac_g_];
    Index current_nz = 0;
    for (Index i=0; i<nz_full_jac_g_; i++) {
      const Index c_row = c_row_pos[g_iRow[i]-1];
      const Index c_col = c_col_pos[g_jCol[i]-1];
      if (c_col != -1 && c_row != -1) {
        jac_idx_map_[current_nz] = i;
        jac_c_iRow[current_nz] = c_row + 1;
        jac_c_jCol[current_nz] = c_col + 1;
        current_nz++;
      }
    }
    nz_jac_c_ = current_nz;
    Jac_c_space = new GenTMatrixSpace(n_c, n_x_not_fixed, nz_jac_c_, jac_c_iRow, jac_c_jCol);
    delete [] jac_c_iRow;
    jac_c_iRow = NULL;
    delete [] jac_c_jCol;
    jac_c_jCol = NULL;

    const Index* d_col_pos = P_x_full_x_->CompressedPosIndices();
    const Index* d_row_pos = P_d_g_->CompressedPosIndices();
    // ... build the nonzero structure for jac_d
    // ... (the permuation from rows in jac_g to jac_c is the
    // ...  the same as P_d_g_)
    Index* jac_d_iRow = new Index[nz_full_jac_g_];
    Index* jac_d_jCol = new Index[nz_full_jac_g_];
    current_nz = 0;
    for (Index i=0; i<nz_full_jac_g_; i++) {
      const Index d_row = d_row_pos[g_iRow[i]-1];
      const Index d_col = d_col_pos[g_jCol[i]-1];
      if (d_col != -1 && d_row != -1) {
        jac_idx_map_[current_nz + nz_jac_c_] = i;
        jac_d_iRow[current_nz] = d_row + 1;
        jac_d_jCol[current_nz] = d_col + 1;
        current_nz++;
      }
    }
    nz_jac_d_ = current_nz;
    Jac_d_space = new GenTMatrixSpace(n_d, n_x_not_fixed, nz_jac_d_, jac_d_iRow, jac_d_jCol);
    delete [] jac_d_iRow;
    jac_d_iRow = NULL;
    delete [] jac_d_jCol;
    jac_d_jCol = NULL;

    delete [] g_iRow;
    g_iRow = NULL;
    delete [] g_jCol;
    g_jCol = NULL;

    /** Create the matrix space for the hessian of the lagrangian */
    Index* full_h_iRow = new Index[nz_full_h_];
    Index* full_h_jCol = new Index[nz_full_h_];
    Index* h_iRow = new Index[nz_full_h_];
    Index* h_jCol = new Index[nz_full_h_];
    tnlp_->eval_h(n_full_x_, NULL, false, 0, n_full_g_, NULL, false,
                  nz_full_h_, full_h_iRow, full_h_jCol, NULL);

    h_idx_map_ = new Index[nz_full_h_];
    const Index* h_pos = P_x_full_x_->CompressedPosIndices();
    current_nz = 0;
    for (Index i=0; i<nz_full_h_; i++) {
      const Index h_row = h_pos[full_h_iRow[i]-1];
      const Index h_col = h_pos[full_h_jCol[i]-1];
      if (h_row != -1 && h_col != -1) {
        h_idx_map_[current_nz] = i;
        h_iRow[current_nz] = h_row + 1;
        h_jCol[current_nz] = h_col + 1;
        current_nz++;
      }
    }
    nz_h_ = current_nz;
    Hess_lagrangian_space = new SymTMatrixSpace(n_x_not_fixed, nz_h_, h_iRow, h_jCol);
    delete [] full_h_iRow;
    full_h_iRow = NULL;
    delete [] full_h_jCol;
    full_h_jCol = NULL;
    delete [] h_iRow;
    h_iRow = NULL;
    delete [] h_jCol;
    h_jCol = NULL;

    return true;
  }

  bool TNLPAdapter::GetBoundsInformation(Matrix& Px_L,
                                         Vector& x_L,
                                         Matrix& Px_U,
                                         Vector& x_U,
                                         Matrix& Pd_L,
                                         Vector& d_L,
                                         Matrix& Pd_U,
                                         Vector& d_U)
  {
    // This could be done more efficiently, I have already called this method
    // once to setup the structure for the problem, I could store the values
    // and use them here ?
    Number* x_l = new Number[n_full_x_];
    Number* x_u = new Number[n_full_x_];
    Number* g_l = new Number[n_full_g_];
    Number* g_u = new Number[n_full_g_];
    tnlp_->get_bounds_info(n_full_x_, x_l, x_u, n_full_g_, g_l, g_u);

    // Set the bounds values for x
    DenseVector* dx_L = dynamic_cast<DenseVector*>(&x_L);
    DBG_ASSERT(dx_L);
    Number* values = dx_L->Values();
    ExpansionMatrix* em_Px_L = dynamic_cast<ExpansionMatrix*>(&Px_L);
    DBG_ASSERT(em_Px_L);
    for (Index i=0; i<Px_L.NCols(); i++) {
      Index ipopt_idx = em_Px_L->ExpandedPosIndices()[i];
      Index full_idx = P_x_full_x_->ExpandedPosIndices()[ipopt_idx];
      Number lower_bound = x_l[full_idx];
      values[i] = lower_bound;
    }

    DenseVector* dx_U = dynamic_cast<DenseVector*>(&x_U);
    DBG_ASSERT(dx_U);
    values = dx_U->Values();
    ExpansionMatrix* em_Px_U = dynamic_cast<ExpansionMatrix*>(&Px_U);
    DBG_ASSERT(em_Px_U);
    for (Index i=0; i<Px_U.NCols(); i++) {
      Index ipopt_idx = em_Px_U->ExpandedPosIndices()[i];
      Index full_idx = P_x_full_x_->ExpandedPosIndices()[ipopt_idx];
      Number upper_bound = x_u[full_idx];
      values[i] = upper_bound;
    }

    // get the bounds values (rhs values to subtract) for c
    // i.e. if gL == gU, then we actually have g(x) = gL = gU,
    // since we solve c(x) = 0, we actually need c(x) - gL = 0
    for (Index i=0; i<P_c_g_->NCols(); i++) {
      Index full_idx = P_c_g_->ExpandedPosIndices()[i];
      Number rhs = g_l[full_idx];
      c_rhs_[i] = rhs;
    }

    // get the bounds values for d
    DenseVector* dd_L = dynamic_cast<DenseVector*>(&d_L);
    DBG_ASSERT(dd_L);
    values = dd_L->Values();
    ExpansionMatrix* em_Pd_L = dynamic_cast<ExpansionMatrix*>(&Pd_L);
    DBG_ASSERT(em_Pd_L);
    for (Index i=0; i<Pd_L.NCols(); i++) {
      Index d_exp_idx = em_Pd_L->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number lower_bound = g_l[full_idx];
      values[i] = lower_bound;
    }

    DenseVector* dd_U = dynamic_cast<DenseVector*>(&d_U);
    DBG_ASSERT(dd_U);
    values = dd_U->Values();
    ExpansionMatrix* em_Pd_U = dynamic_cast<ExpansionMatrix*>(&Pd_U);
    DBG_ASSERT(em_Pd_U);
    for (Index i=0; i<Pd_U.NCols(); i++) {
      Index d_exp_idx = em_Pd_U->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number upper_bound = g_u[full_idx];
      values[i] = upper_bound;
      //      printf("d_U[%d] = %lf\n", ampl_idx, upper_bound);
    }

    delete [] x_l;
    x_l = NULL;
    delete [] x_u;
    x_u = NULL;
    delete [] g_l;
    g_l = NULL;
    delete [] g_u;
    g_u = NULL;

    return true;
  }

  bool TNLPAdapter::GetStartingPoint(Vector& x,
                                     bool need_x,
                                     Vector& y_c,
                                     bool need_y_c,
                                     Vector& y_d,
                                     bool need_y_d,
                                     Vector& z_L,
                                     bool need_z_L,
                                     Vector& z_U,
                                     bool need_z_U,
                                     Vector& v_L,
                                     bool need_v_L,
                                     Vector& v_U,
                                     bool need_v_U
                                    )
  {
    Number* full_x = new Number[n_full_x_];
    Number* full_z_l = new Number[n_full_x_];
    Number* full_z_u = new Number[n_full_x_];
    Number* full_lambda = new Number[n_full_g_];
    bool init_x = need_x;
    bool init_z = need_z_L && need_z_U;
    bool init_lambda = need_y_c && need_y_d;

    bool retvalue =
      tnlp_->get_starting_point(n_full_x_, init_x, full_x, init_z, full_z_l, full_z_u, n_full_g_, init_lambda, full_lambda);

    if (!retvalue) {
      return false;
    }

    DenseVector* dx = dynamic_cast<DenseVector*>(&x);
    DBG_ASSERT(dx);
    Number* values = dx->Values();

    if (need_x) {
      const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
      for (Index i=0; i<x.Dim(); i++) {
        values[i] = full_x[x_pos[i]];
      }
    }

    if (need_y_c) {
      DenseVector* dy_c = dynamic_cast<DenseVector*>(&y_c);
      DBG_ASSERT(dy_c);
      values = dy_c->Values();
      const Index* y_c_pos = P_c_g_->ExpandedPosIndices();
      for (Index i=0; i<y_c.Dim(); i++) {
        values[i] = full_lambda[y_c_pos[i]];
      }
    }

    if (need_y_d) {
      DenseVector* dy_d = dynamic_cast<DenseVector*>(&y_d);
      DBG_ASSERT(dy_d);
      values = dy_d->Values();
      const Index* y_d_pos = P_d_g_->ExpandedPosIndices();
      for (Index i=0; i<y_d.Dim(); i++) {
        values[i] = full_lambda[y_d_pos[i]];
      }
    }

    if (need_z_L) {
      DenseVector* dz_l = dynamic_cast<DenseVector*>(&z_L);
      DBG_ASSERT(dz_l);
      values = dz_l->Values();
      const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
      const Index* z_l_pos = P_x_x_L_->ExpandedPosIndices();
      for (Index i=0; i<z_L.Dim(); i++) {
        Index idx = z_l_pos[i]; // convert from x_L to x (ipopt)
        idx = x_pos[idx]; // convert from x (ipopt) to x_full
        values[i] = full_z_l[idx];
      }
    }

    if (need_z_U) {
      DenseVector* dz_u = dynamic_cast<DenseVector*>(&z_U);
      DBG_ASSERT(dz_u);
      values = dz_u->Values();
      const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
      const Index* z_u_pos = P_x_x_L_->ExpandedPosIndices();
      for (Index i=0; i<z_U.Dim(); i++) {
        Index idx = z_u_pos[i]; // convert from x_u to x (ipopt)
        idx = x_pos[idx]; // convert from x (ipopt) to x_full
        values[i] = full_z_u[idx];
      }
    }

    DBG_ASSERT(!need_v_L && !need_v_U
               && "Need to think about initialization of these.");

    delete [] full_x;
    full_x = NULL;
    delete [] full_z_l;
    full_z_l = NULL;
    delete [] full_z_u;
    full_z_u = NULL;
    delete [] full_lambda;
    full_lambda = NULL;

    return true;
  }

  bool TNLPAdapter::Eval_f(const Vector& x, Number& f)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }
    return tnlp_->eval_f(n_full_x_, full_x_, new_x, f);
  }

  bool TNLPAdapter::Eval_grad_f(const Vector& x, Vector& g_f)
  {
    bool retvalue = false;
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    Number* full_grad_f = new Number[n_full_x_];
    if (tnlp_->eval_grad_f(n_full_x_, full_x_, new_x, full_grad_f)) {
      const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
      DenseVector* dg_f = dynamic_cast<DenseVector*>(&g_f);
      DBG_ASSERT(dg_f);
      Number* values = dg_f->Values();

      for (Index i=0; i<g_f.Dim(); i++) {
        values[i] = full_grad_f[x_pos[i]];
      }
      retvalue = true;
    }

    delete [] full_grad_f;
    return retvalue;
  }

  bool TNLPAdapter::Eval_c(const Vector& x, Vector& c)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    DenseVector* dc = dynamic_cast<DenseVector*>(&c);
    DBG_ASSERT(dc);
    Number* values = dc->Values();
    if (internal_eval_g(new_x)) {
      const Index* c_pos = P_c_g_->ExpandedPosIndices();
      for (Index i=0; i<c.Dim(); i++) {
        values[i] = full_g_[c_pos[i]];
        values[i] -= c_rhs_[i];
      }
      return true;
    }

    return false;
  }

  bool TNLPAdapter::Eval_jac_c(const Vector& x, Matrix& jac_c)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_jac_g(new_x)) {
      GenTMatrix* gt_jac_c = dynamic_cast<GenTMatrix*>(&jac_c);
      DBG_ASSERT(gt_jac_c);
      Number* values = gt_jac_c->Values();

      for (Index i=0; i<nz_jac_c_; i++) {
        // Assume the same structure as initially given
        values[i] = jac_g_[jac_idx_map_[i]];
      }
      return true;
    }
    return false;
  }

  bool TNLPAdapter::Eval_d(const Vector& x, Vector& d)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    DenseVector* dd = dynamic_cast<DenseVector*>(&d);
    DBG_ASSERT(dd);
    Number* values = dd->Values();
    if (internal_eval_g(new_x)) {
      const Index* d_pos = P_d_g_->ExpandedPosIndices();
      for (Index i=0; i<d.Dim(); i++) {
        values[i] = full_g_[d_pos[i]];
      }
      return true;
    }

    return false;
  }

  bool TNLPAdapter::Eval_jac_d(const Vector& x, Matrix& jac_d)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_jac_g(new_x)) {
      GenTMatrix* gt_jac_d = dynamic_cast<GenTMatrix*>(&jac_d);
      DBG_ASSERT(gt_jac_d);
      Number* values = gt_jac_d->Values();

      for (Index i=0; i<nz_jac_d_; i++) {
        // Assume the same structure as initially given
        values[i] = jac_g_[jac_idx_map_[nz_jac_c_ + i]];
      }
      return true;
    }
    return false;
  }

  bool TNLPAdapter::Eval_h(const Vector& x,
                           Number obj_factor,
                           const Vector& yc,
                           const Vector& yd,
                           SymMatrix& h)
  {
    // First see if all weights are set to zero (for example, when
    // computing the least square multiplier estimates, this is what
    // we do).  In that case, there is no need to compute values, just
    // set them to zero.
    if (obj_factor==0. && yc.Asum()==0. && yd.Asum()==0.) {
      SymTMatrix* st_h = dynamic_cast<SymTMatrix*>(&h);
      DBG_ASSERT(st_h);
      Number* values = st_h->Values();
      for (Index i=0; i<nz_h_; i++) {
        values[i] = 0.;
      }
      return true;
    }

    bool retval = false;
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }
    bool new_y = false;
    if (update_local_lambda(yc, yd)) {
      new_y = true;
    }

    Number* full_h = new Number[nz_full_h_];

    if (tnlp_->eval_h(n_full_x_, full_x_, new_x, obj_factor, n_full_g_,
                      full_lambda_, new_y, nz_full_h_, NULL, NULL, full_h)) {
      SymTMatrix* st_h = dynamic_cast<SymTMatrix*>(&h);
      DBG_ASSERT(st_h);
      Number* values = st_h->Values();
      for (Index i=0; i<nz_h_; i++) {
        values[i] = full_h[h_idx_map_[i]];
      }
      retval = true;
    }

    delete [] full_h;
    return retval;
  }

  void TNLPAdapter::FinalizeSolution(ApplicationReturnStatus status,
                                     const Vector& x, const Vector& z_L, const Vector& z_U,
                                     const Vector& c, const Vector& d,
                                     const Vector& y_c, const Vector& y_d,
                                     Number obj_value)
  {
    update_local_x(x);
    update_local_lambda(y_c, y_d);

    ResortX(x, full_x_);
    ResortG(y_c, y_d, full_lambda_);

    const DenseVector* dc = dynamic_cast<const DenseVector*>(&c);
    const DenseVector* dd = dynamic_cast<const DenseVector*>(&d);
    DBG_ASSERT(dc && dd);
    Number* full_g = new Number[n_full_g_];
    ResortG(c, d, full_g);

    Number* full_z_L = new Number[n_full_x_];
    Number* full_z_U = new Number[n_full_x_];
    for (int i=0; i<n_full_x_; i++) {
      full_z_L[i] = nlp_lower_bound_inf_;
      full_z_U[i] = nlp_upper_bound_inf_;
    }
    ResortBnds(z_L, full_z_L, z_U, full_z_U);

    tnlp_->finalize_solution(status,
                             n_full_x_, full_x_, full_z_L, full_z_U,
                             n_full_g_, full_g, full_lambda_,
                             obj_value);

    delete [] full_z_L;
    full_z_L = NULL;
    delete [] full_z_U;
    full_z_U = NULL;
    delete [] full_g;
    full_g = NULL;
  }

  void TNLPAdapter::ResortX(const Vector& x, Number* x_orig)
  {
    const DenseVector* dx = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dx);
    const Number* x_values = dx->Values();

    const Index* x_pos = P_x_full_x_->CompressedPosIndices();
    for (Index i=0; i<n_full_x_; i++) {
      Index idx = x_pos[i];
      if (idx != -1) {
        x_orig[i] = x_values[idx];
      }
      else {
        x_orig[i] = full_x_[i];
      }
    }
  }

  void TNLPAdapter::ResortG(const Vector& c, const Vector& d, Number* g_orig)
  {
    const DenseVector* dc = dynamic_cast<const DenseVector*>(&c);
    DBG_ASSERT(dc);
    const Number* c_values = dc->Values();

    const Index* c_pos = P_c_g_->ExpandedPosIndices();
    for (Index i=0; i<c.Dim(); i++) {
      g_orig[c_pos[i]] = c_values[i];
    }

    const DenseVector* dd = dynamic_cast<const DenseVector*>(&d);
    DBG_ASSERT(dd);
    const Number* d_values = dd->Values();

    const Index* d_pos = P_d_g_->ExpandedPosIndices();
    for (Index i=0; i<d.Dim(); i++) {
      g_orig[d_pos[i]] = d_values[i];
    }
  }

  void TNLPAdapter::ResortBnds(const Vector& x_L, Number* x_L_orig,
                               const Vector& x_U, Number* x_U_orig)
  {
    if (x_L_orig) {
      const DenseVector* dx_L = dynamic_cast<const DenseVector*>(&x_L);
      DBG_ASSERT(dx_L);
      const Number* x_L_values = dx_L->Values();

      const Index* bnds_pos_not_fixed = P_x_x_L_->ExpandedPosIndices();
      const Index* bnds_pos_full = P_x_full_x_->ExpandedPosIndices();
      for (Index i=0; i<x_L.Dim(); i++) {
        int idx = bnds_pos_not_fixed[i];
        idx = bnds_pos_full[idx];
        x_L_orig[idx] = x_L_values[i];
      }
    }

    if (x_U_orig) {
      const DenseVector* dx_U = dynamic_cast<const DenseVector*>(&x_U);
      DBG_ASSERT(dx_U);
      const Number* x_U_values = dx_U->Values();

      const Index* bnds_pos_not_fixed = P_x_x_U_->ExpandedPosIndices();
      const Index* bnds_pos_full = P_x_full_x_->ExpandedPosIndices();
      for (Index i=0; i<x_U.Dim(); i++) {
        int idx = bnds_pos_not_fixed[i];
        idx = bnds_pos_full[idx];
        x_U_orig[idx] = x_U_values[i];
      }
    }
  }

  bool TNLPAdapter::update_local_x(const Vector& x)
  {
    if (x.GetTag() == x_tag_for_iterates_) {
      return false;
    }

    ResortX(x, full_x_);

    x_tag_for_iterates_ = x.GetTag();

    return true;
  }

  bool TNLPAdapter::update_local_lambda(const Vector& y_c, const Vector& y_d)
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

  bool TNLPAdapter::internal_eval_g(bool new_x)
  {
    if (x_tag_for_g_ == x_tag_for_iterates_) {
      // already calculated!
      return true;
    }

    x_tag_for_g_ = x_tag_for_iterates_;
    return tnlp_->eval_g(n_full_x_, full_x_, new_x, n_full_g_, full_g_);
  }

  bool TNLPAdapter::internal_eval_jac_g(bool new_x)
  {
    if (x_tag_for_jac_g_ == x_tag_for_iterates_) {
      // already calculated!
      return true;
    }

    x_tag_for_jac_g_ = x_tag_for_iterates_;
    return tnlp_->eval_jac_g(n_full_x_, full_x_, new_x, n_full_g_,                             nz_full_jac_g_, NULL, NULL, jac_g_);
  }

} // namespace Ipopt
