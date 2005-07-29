// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpOrigIpoptNLP.hpp"
#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  DefineIpoptType(OrigIpoptNLP);

  OrigIpoptNLP::OrigIpoptNLP(const SmartPtr<const Journalist>& jnlst,
                             const SmartPtr<NLP>& nlp,
                             const SmartPtr<NLPScalingObject>& nlp_scaling)
      :
      IpoptNLP(nlp_scaling),
      jnlst_(jnlst),
      nlp_(nlp),
      f_cache_(1),
      grad_f_cache_(1),
      c_cache_(1),
      jac_c_cache_(1),
      d_cache_(1),
      jac_d_cache_(1),
      h_cache_(1),
      f_evals_(0),
      grad_f_evals_(0),
      c_evals_(0),
      jac_c_evals_(0),
      d_evals_(0),
      jac_d_evals_(0),
      h_evals_(0),
      initialized_(false)
  {}

  OrigIpoptNLP::~OrigIpoptNLP()
  {}

  void OrigIpoptNLP::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "bound_relax_factor",
      "Factor for initial relaxation of the bounds.",
      0, false,
      1e-8,
      "Before start of the optimization, the bounds given by the user are "
      "relaxed.  This option determines the factor by how much.  If it "
      "is set to zero, then this option is disabled.  (See Eqn.(35) in "
      "implmentation paper.)");
  }

  bool OrigIpoptNLP::Initialize(const Journalist& jnlst,
                                const OptionsList& options,
                                const std::string& prefix)
  {
    options.GetNumericValue("bound_relax_factor", bound_relax_factor_, prefix);

    if (!nlp_->ProcessOptions(options, prefix)) {
      return false;
    }

    initialized_ = true;
    return IpoptNLP::Initialize(jnlst, options, prefix);
  }

  bool OrigIpoptNLP::InitializeStructures(SmartPtr<Vector>& x,
                                          bool init_x,
                                          SmartPtr<Vector>& y_c,
                                          bool init_y_c,
                                          SmartPtr<Vector>& y_d,
                                          bool init_y_d,
                                          SmartPtr<Vector>& z_L,
                                          bool init_z_L,
                                          SmartPtr<Vector>& z_U,
                                          bool init_z_U,
                                          SmartPtr<Vector>& v_L,
                                          SmartPtr<Vector>& v_U
                                         )
  {
    DBG_ASSERT(initialized_);

    bool retValue = nlp_->GetSpaces(x_space_, c_space_, d_space_,
                                    x_l_space_, px_l_space_,
                                    x_u_space_, px_u_space_,
                                    d_l_space_, pd_l_space_,
                                    d_u_space_, pd_u_space_,
                                    jac_c_space_, jac_d_space_,
                                    h_space_);

    if (!retValue) {
      return false;
    }

    NLP_scaling()->DetermineScaling(x_space_,
                                    c_space_, d_space_,
                                    jac_c_space_, jac_d_space_,
                                    h_space_,
                                    scaled_jac_c_space_, scaled_jac_d_space_,
                                    scaled_h_space_);

    ASSERT_EXCEPTION(x_space_->Dim() >= c_space_->Dim(), TOO_FEW_DOF,
                     "Too few degrees of freedom!");

    // cannot have any null pointers, want zero length vectors
    // instead of null - this will later need to be changed for _h;
    retValue = (IsValid(x_space_) && IsValid(c_space_) && IsValid(d_space_)
                && IsValid(x_l_space_) && IsValid(px_l_space_)
                && IsValid(x_u_space_) && IsValid(px_u_space_)
                && IsValid(d_u_space_) && IsValid(pd_u_space_)
                && IsValid(d_l_space_) && IsValid(pd_l_space_)
                && IsValid(jac_c_space_) && IsValid(jac_d_space_)
                && IsValid(h_space_)
                && IsValid(scaled_jac_c_space_) && IsValid(scaled_jac_d_space_)
                && IsValid(scaled_h_space_));

    DBG_ASSERT(retValue && "Model cannot return null vector or matrix prototypes or spaces,"
               " please return zero length vectors instead");

    // Create the bounds structures
    SmartPtr<Vector> x_L = x_l_space_->MakeNew();
    SmartPtr<Matrix> Px_L = px_l_space_->MakeNew();
    SmartPtr<Vector> x_U = x_u_space_->MakeNew();
    SmartPtr<Matrix> Px_U = px_u_space_->MakeNew();
    SmartPtr<Vector> d_L = d_l_space_->MakeNew();
    SmartPtr<Matrix> Pd_L = pd_l_space_->MakeNew();
    SmartPtr<Vector> d_U = d_u_space_->MakeNew();
    SmartPtr<Matrix> Pd_U = pd_u_space_->MakeNew();

    retValue = nlp_->GetBoundsInformation(*Px_L, *x_L, *Px_U, *x_U,
                                          *Pd_L, *d_L, *Pd_U, *d_U);

    if (!retValue) {
      return false;
    }

    relax_bounds(-bound_relax_factor_, *x_L);
    relax_bounds( bound_relax_factor_, *x_U);
    relax_bounds(-bound_relax_factor_, *d_L);
    relax_bounds( bound_relax_factor_, *d_U);

    x_L_ = ConstPtr(x_L);
    Px_L_ = ConstPtr(Px_L);
    x_U_ = ConstPtr(x_U);
    Px_U_ = ConstPtr(Px_U);
    d_L_ = ConstPtr(d_L);
    Pd_L_ = ConstPtr(Pd_L);
    d_U_ = ConstPtr(d_U);
    Pd_U_ = ConstPtr(Pd_U);

    // now create and store the scaled bounds
    x_L_ = NLP_scaling()->apply_vector_scaling_x_LU(*Px_L_, x_L_, *x_space_);
    x_U_ = NLP_scaling()->apply_vector_scaling_x_LU(*Px_U_, x_U_, *x_space_);
    d_L_ = NLP_scaling()->apply_vector_scaling_d_LU(*Pd_L_, d_L_, *d_space_);
    d_U_ = NLP_scaling()->apply_vector_scaling_d_LU(*Pd_U_, d_U_, *d_space_);

    // Create the iterates structures
    x = x_space_->MakeNew();
    y_c = c_space_->MakeNew();
    y_d = d_space_->MakeNew();
    z_L = x_l_space_->MakeNew();
    z_U = x_u_space_->MakeNew();
    v_L = d_l_space_->MakeNew();
    v_U = d_u_space_->MakeNew();

    retValue = nlp_->GetStartingPoint(GetRawPtr(x), init_x,
                                      GetRawPtr(y_c), init_y_c,
                                      GetRawPtr(y_d), init_y_d,
                                      GetRawPtr(z_L), init_z_L,
                                      GetRawPtr(z_U), init_z_U);

    if (!retValue) {
      return false;
    }

    Number obj_scal = NLP_scaling()->apply_obj_scaling(1.);
    if (init_x) {
      if (NLP_scaling()->have_x_scaling()) {
        x = NLP_scaling()->apply_vector_scaling_x_NonConst(ConstPtr(x));
      }
    }
    if (init_y_c) {
      if (NLP_scaling()->have_c_scaling()) {
        y_c = NLP_scaling()->unapply_vector_scaling_c_NonConst(ConstPtr(y_c));
      }
      if (obj_scal!=1.) {
        y_c->Scal(obj_scal);
      }
    }
    if (init_y_d) {
      if (NLP_scaling()->have_d_scaling()) {
        y_d = NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(y_d));
      }
      if (obj_scal!=1.) {
        y_d->Scal(obj_scal);
      }
    }
    if (init_z_L) {
      if (NLP_scaling()->have_x_scaling()) {
        z_L = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_L_, ConstPtr(z_L), *x_space_);
      }
      if (obj_scal!=1.) {
        z_L->Scal(obj_scal);
      }
    }
    if (init_z_U) {
      if (NLP_scaling()->have_x_scaling()) {
        z_U = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_U_, ConstPtr(z_U), *x_space_);
      }
      if (obj_scal!=1.) {
        z_U->Scal(obj_scal);
      }
    }

    return true;
  }

  void
  OrigIpoptNLP::relax_bounds(Number bound_relax_factor, Vector& bounds)
  {
    DBG_START_METH("OrigIpoptNLP::relax_bounds", dbg_verbosity);
    if (bound_relax_factor!=0.) {
      SmartPtr<Vector> tmp = bounds.MakeNew();
      tmp->Copy(bounds);
      tmp->ElementWiseAbs();
      SmartPtr<Vector> ones = bounds.MakeNew();
      ones->Set(1.);
      tmp->ElementWiseMax(*ones);
      DBG_PRINT((1, "bound_relax_factor = %e", bound_relax_factor));
      DBG_PRINT_VECTOR(2, "tmp", *tmp);
      DBG_PRINT_VECTOR(2, "bounds before", bounds);
      bounds.Axpy(bound_relax_factor, *tmp);
      DBG_PRINT_VECTOR(2, "bounds after", bounds);
    }
  }

  Number OrigIpoptNLP::f(const Vector& x)
  {
    DBG_START_METH("OrigIpoptNLP::f", dbg_verbosity);
    Number ret = 0.0;
    DBG_PRINT((2, "x.Tag = %d\n", x.GetTag()));
    if (!f_cache_.GetCachedResult1Dep(ret, &x)) {
      f_evals_++;
      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      bool success = nlp_->Eval_f(*unscaled_x, ret);
      DBG_PRINT((1, "success = %d ret = %e\n", success, ret));
      ASSERT_EXCEPTION(success && FiniteNumber(ret), Eval_Error,
                       "Error evaluating the objective function");
      ret = NLP_scaling()->apply_obj_scaling(ret);
      f_cache_.AddCachedResult1Dep(ret, &x);
    }

    return ret;
  }

  SmartPtr<const Vector> OrigIpoptNLP::grad_f(const Vector& x)
  {
    SmartPtr<Vector> unscaled_grad_f;
    SmartPtr<const Vector> retValue;
    if (!grad_f_cache_.GetCachedResult1Dep(retValue, &x)) {
      grad_f_evals_++;
      unscaled_grad_f = x_space_->MakeNew();

      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      bool success = nlp_->Eval_grad_f(*unscaled_x, *unscaled_grad_f);
      ASSERT_EXCEPTION(success && FiniteNumber(unscaled_grad_f->Nrm2()),
                       Eval_Error, "Error evaluating the gradient of the objective function");
      retValue = NLP_scaling()->apply_grad_obj_scaling(ConstPtr(unscaled_grad_f));
      grad_f_cache_.AddCachedResult1Dep(retValue, &x);
    }

    return retValue;
  }

  /** Equality constraint residual */
  SmartPtr<const Vector> OrigIpoptNLP::c(const Vector& x)
  {
    SmartPtr<const Vector> retValue;
    SmartPtr<const Vector> dep;
    if (c_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Vector has always the same tag (this might make a difference
      // in cases where only the constraints are supposed to change...
      dep = NULL;
      if (!c_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        retValue = c_space_->MakeNew();
        c_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      dep = &x;
    }
    if (!c_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
      SmartPtr<Vector> unscaled_c = c_space_->MakeNew();
      c_evals_++;
      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      bool success = nlp_->Eval_c(*unscaled_x, *unscaled_c);
      ASSERT_EXCEPTION(success && FiniteNumber(unscaled_c->Nrm2()),
                       Eval_Error, "Error evaluating the equality constraints");
      retValue = NLP_scaling()->apply_vector_scaling_c(ConstPtr(unscaled_c));
      c_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
    }

    return retValue;
  }

  SmartPtr<const Vector> OrigIpoptNLP::d(const Vector& x)
  {
    DBG_START_METH("OrigIpoptNLP::d", dbg_verbosity);
    SmartPtr<const Vector> retValue;
    SmartPtr<const Vector> dep;
    if (d_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Vector has always the same tag (this might make a difference
      // in cases where only the constraints are supposed to change...
      dep = NULL;
      if (!d_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        retValue = d_space_->MakeNew();
        d_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      dep = &x;
    }
    if (!d_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
      d_evals_++;
      SmartPtr<Vector> unscaled_d = d_space_->MakeNew();

      DBG_PRINT_VECTOR(2, "scaled_x", x);
      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      bool success = nlp_->Eval_d(*unscaled_x, *unscaled_d);
      DBG_PRINT_VECTOR(2, "unscaled_d", *unscaled_d);
      ASSERT_EXCEPTION(success && FiniteNumber(unscaled_d->Nrm2()),
                       Eval_Error, "Error evaluating the inequality constraints");
      retValue = NLP_scaling()->apply_vector_scaling_d(ConstPtr(unscaled_d));
      d_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
    }

    return retValue;
  }

  SmartPtr<const Matrix> OrigIpoptNLP::jac_c(const Vector& x)
  {
    SmartPtr<const Matrix> retValue;
    SmartPtr<const Vector> dep;
    if (c_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Matrix has always the same tag
      dep = NULL;
      if (!jac_c_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        retValue = jac_c_space_->MakeNew();
        jac_c_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      dep = &x;
    }
    if (!jac_c_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
      jac_c_evals_++;
      SmartPtr<Matrix> unscaled_jac_c = jac_c_space_->MakeNew();

      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      bool success = nlp_->Eval_jac_c(*unscaled_x, *unscaled_jac_c);
      ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the jacobian of the equality constraints");
      retValue = NLP_scaling()->apply_jac_c_scaling(ConstPtr(unscaled_jac_c));
      jac_c_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
    }

    return retValue;
  }

  SmartPtr<const Matrix> OrigIpoptNLP::jac_d(const Vector& x)
  {
    SmartPtr<const Matrix> retValue;
    SmartPtr<const Vector> dep;
    if (d_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Matrix has always the same tag
      dep = NULL;
      if (!jac_d_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        retValue = jac_d_space_->MakeNew();
        jac_d_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      dep = &x;
    }
    if (!jac_d_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
      jac_d_evals_++;
      SmartPtr<Matrix> unscaled_jac_d = jac_d_space_->MakeNew();

      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      bool success = nlp_->Eval_jac_d(*unscaled_x, *unscaled_jac_d);
      ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the jacobian of the inequality constraints");
      retValue = NLP_scaling()->apply_jac_d_scaling(ConstPtr(unscaled_jac_d));
      jac_d_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
    }

    return retValue;
  }

  SmartPtr<const SymMatrix> OrigIpoptNLP::h(const Vector& x,
      Number obj_factor,
      const Vector& yc,
      const Vector& yd)
  {
    std::vector<const TaggedObject*> deps(3);
    deps[0] = &x;
    deps[1] = &yc;
    deps[2] = &yd;
    std::vector<Number> scalar_deps(1);
    scalar_deps[0] = obj_factor;

    SmartPtr<SymMatrix> unscaled_h;
    SmartPtr<const SymMatrix> retValue;
    if (!h_cache_.GetCachedResult(retValue, deps, scalar_deps)) {
      h_evals_++;
      unscaled_h = h_space_->MakeNewSymMatrix();

      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      SmartPtr<const Vector> unscaled_yc = NLP_scaling()->apply_vector_scaling_c(&yc);
      SmartPtr<const Vector> unscaled_yd = NLP_scaling()->apply_vector_scaling_d(&yd);
      Number scaled_obj_factor = NLP_scaling()->apply_obj_scaling(obj_factor);
      bool success = nlp_->Eval_h(*unscaled_x, scaled_obj_factor, *unscaled_yc, *unscaled_yd, *unscaled_h);
      ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the hessian of the lagrangian");
      retValue = NLP_scaling()->apply_hessian_scaling(ConstPtr(unscaled_h));
      h_cache_.AddCachedResult(retValue, deps, scalar_deps);
    }

    return retValue;
  }

  void OrigIpoptNLP::GetSpaces(SmartPtr<const VectorSpace>& x_space,
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
    // Make sure that we already have all the pointers
    DBG_ASSERT(IsValid(x_space_) &&
               IsValid(c_space_) &&
               IsValid(d_space_) &&
               IsValid(x_l_space_) &&
               IsValid(px_l_space_) &&
               IsValid(x_u_space_) &&
               IsValid(px_u_space_) &&
               IsValid(d_l_space_) &&
               IsValid(pd_l_space_) &&
               IsValid(d_u_space_) &&
               IsValid(pd_u_space_) &&
               IsValid(scaled_jac_c_space_) &&
               IsValid(scaled_jac_d_space_) &&
               IsValid(scaled_h_space_));

    DBG_ASSERT(IsValid(NLP_scaling()));

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
    Jac_c_space = scaled_jac_c_space_;
    Jac_d_space = scaled_jac_d_space_;
    Hess_lagrangian_space = scaled_h_space_;
  }

  void OrigIpoptNLP::FinalizeSolution(SolverReturn status,
                                      const Vector& x, const Vector& z_L, const Vector& z_U,
                                      const Vector& c, const Vector& d,
                                      const Vector& y_c, const Vector& y_d,
                                      Number obj_value)
  {
    // need to submit the unscaled solution back to the nlp
    SmartPtr<const Vector> unscaled_x =
      NLP_scaling()->unapply_vector_scaling_x(&x);
    SmartPtr<const Vector> unscaled_c =
      NLP_scaling()->unapply_vector_scaling_c(&c);
    SmartPtr<const Vector> unscaled_d =
      NLP_scaling()->unapply_vector_scaling_d(&d);
    const Number unscaled_obj = NLP_scaling()->unapply_obj_scaling(obj_value);

    SmartPtr<const Vector> unscaled_z_L;
    SmartPtr<const Vector> unscaled_z_U;
    SmartPtr<const Vector> unscaled_y_c;
    SmartPtr<const Vector> unscaled_y_d;

    // The objective function scaling factor also appears in the constraints
    Number obj_unscale_factor = NLP_scaling()->unapply_obj_scaling(1.);
    if (obj_unscale_factor!=1.) {
      SmartPtr<Vector> tmp = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_L_, &z_L, *x_space_);
      tmp->Scal(obj_unscale_factor);
      unscaled_z_L = ConstPtr(tmp);

      tmp = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_U_, &z_U, *x_space_);
      tmp->Scal(obj_unscale_factor);
      unscaled_z_U = ConstPtr(tmp);

      tmp = NLP_scaling()->apply_vector_scaling_c_NonConst(&y_c);
      tmp->Scal(obj_unscale_factor);
      unscaled_y_c = ConstPtr(tmp);

      tmp = NLP_scaling()->apply_vector_scaling_d_NonConst(&y_d);
      tmp->Scal(obj_unscale_factor);
      unscaled_y_d = ConstPtr(tmp);
    }
    else {
      unscaled_z_L = NLP_scaling()->apply_vector_scaling_x_LU(*Px_L_, &z_L, *x_space_);
      unscaled_z_U = NLP_scaling()->apply_vector_scaling_x_LU(*Px_U_, &z_U, *x_space_);
      unscaled_y_c = NLP_scaling()->apply_vector_scaling_c(&y_c);
      unscaled_y_d = NLP_scaling()->apply_vector_scaling_d(&y_d);
    }

    nlp_->FinalizeSolution(status, *unscaled_x,
                           *unscaled_z_L, *unscaled_z_U,
                           *unscaled_c, *unscaled_d,
                           *unscaled_y_c, *unscaled_y_d,
                           unscaled_obj);
  }

  void OrigIpoptNLP::AdjustVariableBounds(const Vector& new_x_L, const Vector& new_x_U,
                                          const Vector& new_d_L, const Vector& new_d_U)
  {
    x_L_ = new_x_L.MakeNewCopy();
    x_U_ = new_x_U.MakeNewCopy();
    d_L_ = new_d_L.MakeNewCopy();
    d_U_ = new_d_U.MakeNewCopy();
  }

} // namespace Ipopt
