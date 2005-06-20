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
                             const SmartPtr<NLP>& nlp)
      :
      IpoptNLP(),
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
    roptions->AddLowerBoundedNumberOption("bound_relax_factor","factor for initial relaxation of the bounds",
                                          0, false, 1e-8);
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
    return true;
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
                                          bool init_v_L,
                                          SmartPtr<Vector>& v_U,
                                          bool init_v_U
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

    ASSERT_EXCEPTION(x_space_->Dim() >= c_space_->Dim(), TOO_FEW_DOF,
                     "Too few degrees of freedom!");

    //    ASSERT_EXCEPTION(x_space_->Dim() != c_space_->Dim(), IpoptException,
    //                     "Currently cannot solve a square problem!");

    // cannot have any null pointers, want zero length vectors
    // instead of null - this will later need to be changed for _h;
    retValue = (IsValid(x_space_) && IsValid(c_space_) && IsValid(d_space_)
                && IsValid(x_l_space_) && IsValid(px_l_space_)
                && IsValid(x_u_space_) && IsValid(px_u_space_)
                && IsValid(d_u_space_) && IsValid(pd_u_space_)
                && IsValid(d_l_space_) && IsValid(pd_l_space_)
                && IsValid(jac_c_space_) && IsValid(jac_d_space_)
                && IsValid(h_space_));

    DBG_ASSERT(retValue && "Model cannot return null vector or matrix prototypes or spaces,"
               " please return zero length vectors instead");

    // Create the bounds structures
    x_L_ = x_l_space_->MakeNew();
    Px_L_ = px_l_space_->MakeNew();
    x_U_ = x_u_space_->MakeNew();
    Px_U_ = px_u_space_->MakeNew();
    d_L_ = d_l_space_->MakeNew();
    Pd_L_ = pd_l_space_->MakeNew();
    d_U_ = d_u_space_->MakeNew();
    Pd_U_ = pd_u_space_->MakeNew();

    retValue = nlp_->GetBoundsInformation(*Px_L_, *x_L_, *Px_U_, *x_U_,
                                          *Pd_L_, *d_L_, *Pd_U_, *d_U_);

    if (!retValue) {
      return false;
    }

    relax_bounds(-bound_relax_factor_, *x_L_);
    relax_bounds( bound_relax_factor_, *x_U_);
    relax_bounds(-bound_relax_factor_, *d_L_);
    relax_bounds( bound_relax_factor_, *d_U_);

    // Create the iterates structures
    x = x_space_->MakeNew();
    y_c = c_space_->MakeNew();
    y_d = d_space_->MakeNew();
    z_L = x_l_space_->MakeNew();
    z_U = x_u_space_->MakeNew();
    v_L = d_l_space_->MakeNew();
    v_U = d_u_space_->MakeNew();

    retValue = nlp_->GetStartingPoint(*x, init_x,
                                      *y_c, init_y_c,
                                      *y_d, init_y_d,
                                      *z_L, init_z_L,
                                      *z_U, init_z_U,
                                      *v_L, init_v_L,
                                      *v_U, init_v_U);

    if (!retValue) {
      return false;
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
    DBG_PRINT((1, "x.Tag = %d\n", x.GetTag()));
    if (!f_cache_.GetCachedResult1Dep(ret, &x)) {
      f_evals_++;
      bool success = nlp_->Eval_f(x, ret);
      ASSERT_EXCEPTION(success && FiniteNumber(ret), Eval_Error,
                       "Error evaluating the objective function");
      f_cache_.AddCachedResult1Dep(ret, &x);
    }

    return ret;
  }

  SmartPtr<const Vector> OrigIpoptNLP::grad_f(const Vector& x)
  {
    SmartPtr<Vector> retValue;
    if (!grad_f_cache_.GetCachedResult1Dep(retValue, &x)) {
      grad_f_evals_++;
      retValue = x_space_->MakeNew();

      bool success = nlp_->Eval_grad_f(x, *retValue);
      ASSERT_EXCEPTION(success && FiniteNumber(retValue->Nrm2()),
                       Eval_Error, "Error evaluating the gradient of the objective function");

      grad_f_cache_.AddCachedResult1Dep(retValue, &x);
    }

    return ConstPtr(retValue);
  }


  /** Equality constraint residual */
  SmartPtr<const Vector> OrigIpoptNLP::c(const Vector& x)
  {
    SmartPtr<Vector> retValue;
    if (!c_cache_.GetCachedResult1Dep(retValue, &x)) {
      c_evals_++;
      retValue = c_space_->MakeNew();
      bool success = nlp_->Eval_c(x, *retValue);
      ASSERT_EXCEPTION(success && FiniteNumber(retValue->Nrm2()),
                       Eval_Error, "Error evaluating the equality constraints");
      c_cache_.AddCachedResult1Dep(retValue, &x);
    }

    return ConstPtr(retValue);
  }


  SmartPtr<const Vector> OrigIpoptNLP::d(const Vector& x)
  {
    DBG_START_METH("OrigIpoptNLP::d", dbg_verbosity);
    SmartPtr<Vector> retValue;
    if (!d_cache_.GetCachedResult1Dep(retValue, &x)) {
      d_evals_++;
      retValue = d_space_->MakeNew();

      DBG_PRINT_VECTOR(2, "x", x);
      bool success = nlp_->Eval_d(x, *retValue);
      DBG_PRINT_VECTOR(2, "retValue", *retValue);
      ASSERT_EXCEPTION(success && FiniteNumber(retValue->Nrm2()),
                       Eval_Error, "Error evaluating the inequality constraints");
      d_cache_.AddCachedResult1Dep(retValue, &x);
    }

    return ConstPtr(retValue);
  }

  SmartPtr<const Matrix> OrigIpoptNLP::jac_c(const Vector& x)
  {
    SmartPtr<Matrix> retValue;
    if (!jac_c_cache_.GetCachedResult1Dep(retValue, &x)) {
      jac_c_evals_++;
      retValue = jac_c_space_->MakeNew();

      bool success = nlp_->Eval_jac_c(x, *retValue);
      ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the jacobian of the equality constraints");
      jac_c_cache_.AddCachedResult1Dep(retValue, &x);
    }

    return ConstPtr(retValue);
  }

  SmartPtr<const Matrix> OrigIpoptNLP::jac_d(const Vector& x)
  {
    SmartPtr<Matrix> retValue;
    if (!jac_d_cache_.GetCachedResult1Dep(retValue, &x)) {
      jac_d_evals_++;
      retValue = jac_d_space_->MakeNew();

      bool success = nlp_->Eval_jac_d(x, *retValue);
      ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the jacobian of the inequality constraints");
      jac_d_cache_.AddCachedResult1Dep(retValue, &x);
    }

    return ConstPtr(retValue);
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

    SmartPtr<SymMatrix> retValue;
    if (!h_cache_.GetCachedResult(retValue, deps, scalar_deps)) {
      h_evals_++;
      retValue = h_space_->MakeNewSymMatrix();

      bool success = nlp_->Eval_h(x, obj_factor, yc, yd, *retValue);
      ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the hessian of the lagrangian");
      h_cache_.AddCachedResult(retValue, deps, scalar_deps);
    }

    return ConstPtr(retValue);
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
               IsValid(jac_c_space_) &&
               IsValid(jac_d_space_) &&
               IsValid(h_space_));

    x_space = ConstPtr(x_space_);
    c_space = ConstPtr(c_space_);
    d_space = ConstPtr(d_space_);
    x_l_space = ConstPtr(x_l_space_);
    px_l_space = ConstPtr(px_l_space_);
    x_u_space = ConstPtr(x_u_space_);
    px_u_space = ConstPtr(px_u_space_);
    d_l_space = ConstPtr(d_l_space_);
    pd_l_space = ConstPtr(pd_l_space_);
    d_u_space = ConstPtr(d_u_space_);
    pd_u_space = ConstPtr(pd_u_space_);
    Jac_c_space = ConstPtr(jac_c_space_);
    Jac_d_space = ConstPtr(jac_d_space_);
    Hess_lagrangian_space = ConstPtr(h_space_);
  }

  void OrigIpoptNLP::AdjustVariableBounds(const Vector& new_x_L, const Vector& new_x_U,
                                          const Vector& new_d_L, const Vector& new_d_U)
  {
    x_L_->Copy(new_x_L);
    x_U_->Copy(new_x_U);
    d_L_->Copy(new_d_L);
    d_U_->Copy(new_d_U);
  }

} // namespace Ipopt
