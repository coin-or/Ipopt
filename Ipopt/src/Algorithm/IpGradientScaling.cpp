// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-07-13

#include "IpGradientScaling.hpp"

namespace Ipopt
{

  void GradientScaling::RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "nlp_scaling_max_gradient", "Maximum gradient after NLP scaling.",
      0, true, 100.0,
      "This is the gradient scaling cut-off. If the maximum"
      " gradient is above this value, then gradient based scaling"
      " will be performed. Scaling parameters are calculated to"
      " scale the maximum gradient back to this value. (This is g_max in "
      "Section 3.8 of the implementation paper.) Note: This"
      " option is only used if \"nlp_scaling_method\" is chosen as"
      " \"gradient-based\".");
    roptions->AddLowerBoundedNumberOption(
      "nlp_scaling_obj_target_gradient",
      "Target value for objective function gradient size.",
      0, false, 0.,
      "If a positive number is chosen, the scaling factor the objective "
      "function is computed so that the gradient has the max norm of the given "
      "size at the starting point.  This overrides nlp_scaling_max_gradient "
      "for the objective function.");
    roptions->AddLowerBoundedNumberOption(
      "nlp_scaling_constr_target_gradient",
      "Target value for constraint function gradient size.",
      0, false, 0.,
      "If a positive number is chosen, the scaling factor the constraint "
      "functions is computed so that the gradient has the max norm of the given "
      "size at the starting point.  This overrides nlp_scaling_max_gradient "
      "for the constraint functions.");
    roptions->AddLowerBoundedNumberOption(
      "nlp_scaling_min_value",
      "Minimum value of gradient-based scaling values.",
      0, false, 1e-8,
      "This is the lower bound for the scaling factors computed by "
      "gradient-based scaling method.  If some derivatives of some functions "
      "are huge, the scaling factors will otherwise become very small, and "
      "the (unscaled) final constraint violation, for example, might then be "
      "significant.  Note: This option is only used if \"nlp_scaling_method\" "
      "is chosen as \"gradient-based\".");
  }

  bool GradientScaling::InitializeImpl(const OptionsList& options,
                                       const std::string& prefix)
  {
    options.GetNumericValue("nlp_scaling_max_gradient",
                            scaling_max_gradient_, prefix);
    options.GetNumericValue("nlp_scaling_obj_target_gradient",
                            scaling_obj_target_gradient_, prefix);
    options.GetNumericValue("nlp_scaling_constr_target_gradient",
                            scaling_constr_target_gradient_, prefix);
    options.GetNumericValue("nlp_scaling_min_value",
                            scaling_min_value_, prefix);
    return StandardScalingBase::InitializeImpl(options, prefix);
  }


  void GradientScaling::DetermineScalingParametersImpl(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    const SmartPtr<const MatrixSpace> jac_c_space,
    const SmartPtr<const MatrixSpace> jac_d_space,
    const SmartPtr<const SymMatrixSpace> h_space,
    const Matrix& Px_L, const Vector& x_L,
    const Matrix& Px_U, const Vector& x_U,
    Number& df,
    SmartPtr<Vector>& dx,
    SmartPtr<Vector>& dc,
    SmartPtr<Vector>& dd)
  {
    DBG_ASSERT(IsValid(nlp_));

    SmartPtr<Vector> x = x_space->MakeNew();
    if (!nlp_->GetStartingPoint(GetRawPtr(x), true,
                                NULL, false,
                                NULL, false,
                                NULL, false,
                                NULL, false)) {
      THROW_EXCEPTION(FAILED_INITIALIZATION,
                      "Error getting initial point from NLP in GradientScaling.\n");
    }

    //
    // Calculate grad_f scaling
    //
    SmartPtr<Vector> grad_f = x_space->MakeNew();
    if (nlp_->Eval_grad_f(*x, *grad_f)) {
      double max_grad_f = grad_f->Amax();
      df = 1.;
      if (scaling_obj_target_gradient_ == 0.) {
        if (max_grad_f > scaling_max_gradient_) {
          df = scaling_max_gradient_ / max_grad_f;
        }
      }
      else {
        if (max_grad_f == 0.) {
          Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                         "Gradient of objective function is zero at starting point.  Cannot determine scaling factor based on scaling_obj_target_gradient option.\n");
        }
        else {
          df = scaling_obj_target_gradient_ / max_grad_f;
        }
      }
      df = Max(df, scaling_min_value_);
      Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                     "Scaling parameter for objective function = %e\n", df);
    }
    else {
      Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                     "Error evaluating objective gradient at user provided starting point.\n  No scaling factor for objective function computed!\n");
      df = 1.;
    }
    //
    // No x scaling
    //
    dx = NULL;

    dc = NULL;
    if (c_space->Dim()>0) {
      //
      // Calculate c scaling
      //
      SmartPtr<Matrix> jac_c = jac_c_space->MakeNew();
      if (nlp_->Eval_jac_c(*x, *jac_c)) {
        dc = c_space->MakeNew();
        const double dbl_min = std::numeric_limits<double>::min();
        dc->Set(dbl_min);
        jac_c->ComputeRowAMax(*dc, false);
        Number arow_max = dc->Amax();
        if (scaling_constr_target_gradient_<=0.) {
          if (arow_max > scaling_max_gradient_) {
            dc->ElementWiseReciprocal();
            dc->Scal(scaling_max_gradient_);
            SmartPtr<Vector> dummy = dc->MakeNew();
            dummy->Set(1.);
            dc->ElementWiseMin(*dummy);
          }
          else {
            dc = NULL;
          }
        }
        else {
          dc->Set(scaling_constr_target_gradient_/arow_max);
        }
        if (IsValid(dc) && scaling_min_value_ > 0.) {
          SmartPtr<Vector> tmp = dc->MakeNew();
          tmp->Set(scaling_min_value_);
          dc->ElementWiseMax(*tmp);
        }
      }
      else {
        Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                       "Error evaluating Jacobian of equality constraints at user provided starting point.\n  No scaling factors for equality constraints computed!\n");
      }
    }

    dd = NULL;
    if (d_space->Dim()>0) {
      //
      // Calculate d scaling
      //
      SmartPtr<Matrix> jac_d = jac_d_space->MakeNew();
      if (nlp_->Eval_jac_d(*x, *jac_d)) {
        dd = d_space->MakeNew();
        const double dbl_min = std::numeric_limits<double>::min();
        dd->Set(dbl_min);
        jac_d->ComputeRowAMax(*dd, false);
        Number arow_max = dd->Amax();
        if (scaling_constr_target_gradient_<=0.) {
          if (arow_max > scaling_max_gradient_) {
            dd->ElementWiseReciprocal();
            dd->Scal(scaling_max_gradient_);
            SmartPtr<Vector> dummy = dd->MakeNew();
            dummy->Set(1.);
            dd->ElementWiseMin(*dummy);
          }
          else {
            dd = NULL;
          }
        }
        else {
          dd->Set(scaling_constr_target_gradient_/arow_max);
        }
        if (IsValid(dd) && scaling_min_value_ > 0.) {
          SmartPtr<Vector> tmp = dd->MakeNew();
          tmp->Set(scaling_min_value_);
          dd->ElementWiseMax(*tmp);
        }
      }
      else {
        Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                       "Error evaluating Jacobian of inequality constraints at user provided starting point.\n  No scaling factors for inequality constraints computed!\n");
      }
    }
  }

} // namespace Ipopt
