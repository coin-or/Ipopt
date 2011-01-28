// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-08-31

#include "IpInexactDoglegNormal.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  InexactDoglegNormalStep::InexactDoglegNormalStep(SmartPtr<InexactNewtonNormalStep> newton_step,
      SmartPtr<InexactNormalTerminationTester> normal_tester /* = NULL */)
      :
      InexactNormalStepCalculator(),
      newton_step_(newton_step),
      normal_tester_(normal_tester)
  {}

  InexactDoglegNormalStep::~InexactDoglegNormalStep()
  {}

  void InexactDoglegNormalStep::RegisterOptions(SmartPtr<RegisteredOptions> reg_options)
  {
    reg_options->AddLowerBoundedNumberOption(
      "omega_init",
      "Initial trust region factor for normal problem.",
      0.0, true, 100.);
    reg_options->AddLowerBoundedNumberOption(
      "omega_max",
      "Maximal trust region factor for normal problem.",
      0.0, true, 1e20);
  }

  bool InexactDoglegNormalStep::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("omega_init", curr_omega_, prefix);
    options.GetNumericValue("omega_max", omega_max_, prefix);

    // We do not want to trigger an increase of the trust region
    // factor in the first iteration, so we initialize this flag to
    // flase
    last_tr_inactive_ = true;

    return newton_step_->Initialize(Jnlst(), IpNLP(), IpData(),
                                    IpCq(), options, prefix);
  }

  bool
  InexactDoglegNormalStep::ComputeNormalStep(SmartPtr<Vector>& normal_x,
      SmartPtr<Vector>& normal_s)

  {
    DBG_START_METH("InexactDoglegNormalStep::ComputeNormalStep",
                   dbg_verbosity);

    // test if we should increase the trust region factor
    if (!last_tr_inactive_ && InexData().full_step_accepted()) {
      if (curr_omega_ >= omega_max_) {
        Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                       "Trust region radius factor would be increased, but it is already at its upper limit %e.\n", curr_omega_);
        IpData().Append_info_string("O");
      }
      else {
        Number omega_old = curr_omega_;
        curr_omega_ = Min(omega_max_, 10.*curr_omega_);
        Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                       "Increasing trust region factor from %e to %e\n.", omega_old, curr_omega_);
        IpData().Append_info_string("o");
      }
    }
    last_tr_inactive_ = false;

    // TODO if (IpCq().curr_primal_infeasibility(NORM_2) == 0.) {
    if (IpCq().curr_primal_infeasibility(NORM_2) <= 1e-12 ) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM, "Dogleg step:  We are at a feasible point, the normal step is set to zero.\n");
      normal_x = IpData().curr()->x()->MakeNew();
      normal_s = IpData().curr()->s()->MakeNew();
      normal_x->Set(0.);
      normal_s->Set(0.);
      last_tr_inactive_ = true;
      return true;
    }

    /////////////////// Cauchy Step

    // Compute the negative of the steepest descent direction.
    // (scaled constraint Jacobian transpose times constraint values)
    SmartPtr<const Vector> curr_jac_cdT_times_curr_cdminuss =
      InexCq().curr_jac_cdT_times_curr_cdminuss();
    SmartPtr<const Vector> curr_slack_scaled_d_minus_s =
      InexCq().curr_slack_scaled_d_minus_s();

    DBG_PRINT_VECTOR(1, "curr_jac_cdT_times_curr_cdminuss", *curr_jac_cdT_times_curr_cdminuss);
    DBG_PRINT_VECTOR(1, "curr_slack_scaled_d_minus_s", *curr_slack_scaled_d_minus_s);

    // Compute the norm of the (scaled) gradient of the objective
    // function (A^T c)
    Number v_ATc_norm = InexCq().curr_scaled_Ac_norm();

    // Compute A * A^T * c
    SmartPtr<const Vector> vec_AATc_c =
      IpCq().curr_jac_c_times_vec(*curr_jac_cdT_times_curr_cdminuss);
    SmartPtr<Vector> vec_AATc_d = curr_slack_scaled_d_minus_s->MakeNewCopy();
    vec_AATc_d->ElementWiseMultiply(*InexCq().curr_scaling_slacks());
    DBG_PRINT_VECTOR(1, "curr_scaling_slacks", *InexCq().curr_scaling_slacks());
    DBG_PRINT_VECTOR(1, "vec_AATc_d", *vec_AATc_d);
    vec_AATc_d->AddOneVector(1., *IpCq().curr_jac_d_times_vec(*curr_jac_cdT_times_curr_cdminuss) , 1.);
    DBG_PRINT_VECTOR(1, "IpCq().curr_jac_d_times_vec(*curr_jac_cdT_times_curr_cdminuss)", *IpCq().curr_jac_d_times_vec(*curr_jac_cdT_times_curr_cdminuss));
    DBG_PRINT_VECTOR(1, "vec_AATc_c", *vec_AATc_c);
    DBG_PRINT_VECTOR(1, "vec_AATc_d", *vec_AATc_d);
    Number AATc_norm = IpCq().CalcNormOfType(NORM_2, *vec_AATc_c, *vec_AATc_d);

    // Compute the step size for the Cauchy step
    Number alpha_cs = Min(curr_omega_, v_ATc_norm*v_ATc_norm/(AATc_norm*AATc_norm));
    Jnlst().Printf(J_MOREDETAILED, J_SOLVE_PD_SYSTEM,
                   "Dogleg step: Cauchy step size alpha_cs = %e\n", alpha_cs);
    DBG_PRINT((1, "alpha_cs = %e v_ATc_norm = %e AATc_norm = %e\n", alpha_cs, v_ATc_norm,AATc_norm));

    // Finally get the Cauchy step
    SmartPtr<Vector> v_cauchy_x =
      curr_jac_cdT_times_curr_cdminuss->MakeNewCopy();
    SmartPtr<Vector> v_cauchy_s =
      curr_slack_scaled_d_minus_s->MakeNewCopy();
    v_cauchy_x->Scal(-alpha_cs);
    v_cauchy_s->Scal(alpha_cs);

    // output
    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_SOLVE_PD_SYSTEM)) {
      Jnlst().Printf(J_MOREVECTOR, J_SOLVE_PD_SYSTEM,
                     "Dogleg step: Cauchy step:\n");
      v_cauchy_x->Print(Jnlst(), J_MOREVECTOR, J_SOLVE_PD_SYSTEM, "v_cauchy_x");
      v_cauchy_s->Print(Jnlst(), J_MOREVECTOR, J_SOLVE_PD_SYSTEM, "v_cauchy_s");
    }

    // Compute the objective function reduction of the normal problem
    // for the Cauchy step
    SmartPtr<const Vector> curr_c = IpCq().curr_c();
    SmartPtr<const Vector> curr_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<Vector> inf_c = curr_c->MakeNew();
    SmartPtr<Vector> inf_d = curr_d_minus_s->MakeNew();
    inf_c->AddTwoVectors(1., *curr_c, -alpha_cs, *vec_AATc_c, 0.);
    inf_d->AddTwoVectors(1., *curr_d_minus_s, -alpha_cs, *vec_AATc_d, 0.);
    Number c_Avc_norm_cauchy = IpCq().CalcNormOfType(NORM_2, *inf_c, *inf_d);
    if (IsValid(normal_tester_)) {
      normal_tester_->Set_c_Avc_norm_cauchy(c_Avc_norm_cauchy);
    }
    Number objred_normal_cs = 0.5*(IpCq().CalcNormOfType(NORM_2, *curr_c, *curr_d_minus_s)-c_Avc_norm_cauchy);
    Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                   "Dogleg: Reduction of normal problem objective function by Cauchy step = %23.16e\n", objred_normal_cs);

    // If the Cauchy step already hits the trust region, we are done
    if (alpha_cs == curr_omega_) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM, "Dogleg step:  Cauchy step already hits trust region.\n");
      normal_x = v_cauchy_x;
      normal_s = v_cauchy_s;
      // unscale the slack-based scaling
      normal_s->ElementWiseMultiply(*InexCq().curr_scaling_slacks());
      IpData().Append_info_string("Nc ");
      return true;
    }
    // ToDo: We don't need this if we do a proper check for Newton step below
    SmartPtr<Vector> v_cauchy_x_bak = v_cauchy_x->MakeNewCopy();
    SmartPtr<Vector> v_cauchy_s_bak = v_cauchy_s->MakeNewCopy();

    ///////////////////// Newton Step

    SmartPtr<Vector> v_newton_x = v_cauchy_x->MakeNew();
    SmartPtr<Vector> v_newton_s = v_cauchy_s->MakeNew();
    bool retval =
      newton_step_->ComputeNewtonNormalStep(*v_newton_x, *v_newton_s);
    if (!retval) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM, "Dogleg step: Newton step could not be calculated, return Cauchy step.\n");
      normal_x = v_cauchy_x_bak;
      normal_s = v_cauchy_s_bak;
      // unscale the slack-based scaling
      normal_s->ElementWiseMultiply(*InexCq().curr_scaling_slacks());
      IpData().Append_info_string("NF ");
      return true;
    }
    // output
    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_SOLVE_PD_SYSTEM)) {
      Jnlst().Printf(J_MOREVECTOR, J_SOLVE_PD_SYSTEM,
                     "Dogleg step: Newton step:\n");
      v_newton_x->Print(Jnlst(), J_MOREVECTOR, J_SOLVE_PD_SYSTEM, "v_newton_x");
      v_newton_s->Print(Jnlst(), J_MOREVECTOR, J_SOLVE_PD_SYSTEM, "v_newton_s");
    }

    /////////////////////  Compute the dogleg step

    // Compute the trust region radius
    const Number tr_radius = curr_omega_ * v_ATc_norm;

    // norm of the Newton step
    Number v_newton_norm = IpCq().CalcNormOfType(NORM_2,
                           *v_newton_x, *v_newton_s);
    Jnlst().Printf(J_MOREDETAILED, J_SOLVE_PD_SYSTEM,
                   "Norm of Newton step = %e, trust region radius = %e\n",
                   v_newton_norm, tr_radius);
    if (v_newton_norm <= tr_radius) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM, "Dogleg step:  Newton step is within trust region.\n");
      normal_x = v_newton_x;
      normal_s = v_newton_s;
      last_tr_inactive_ = true;
      IpData().Append_info_string("Nn ");
    }
    else {
      Number v_cauchy_norm = IpCq().CalcNormOfType(NORM_2,
                             *v_cauchy_x, *v_cauchy_s);
      Number v_cs_dot_n = v_newton_x->Dot(*v_cauchy_x) +
                          v_newton_s->Dot(*v_cauchy_s);
      Number a = v_newton_norm*v_newton_norm - 2*v_cs_dot_n +
                 v_cauchy_norm*v_cauchy_norm;
      Number b = 2*(v_cs_dot_n - v_newton_norm*v_newton_norm);
      Number c = v_newton_norm*v_newton_norm - tr_radius*tr_radius;
      Number lambda = (-b-sqrt(b*b-4.*a*c))/(2.*a);

      DBG_PRINT((1, "v_cauchy_norm = %e v_cs_dot_n = %e v_newton_norm = %e\n",v_cauchy_norm,v_cs_dot_n,v_newton_norm));

      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM, "Dogleg step:  Using convex combination of Cauchy and Newton step with factor lambda = %e\n", lambda);
      v_cauchy_x->AddOneVector(1.-lambda, *v_newton_x, lambda);
      v_cauchy_s->AddOneVector(1.-lambda, *v_newton_s, lambda);
      normal_x = v_cauchy_x;
      normal_s = v_cauchy_s;
      IpData().Append_info_string("Nd ");

      DBG_PRINT((1, "v_normal^2  = %e\n",normal_x->Dot(*normal_x)+normal_s->Dot(*normal_s)));
    }

    DBG_PRINT_VECTOR(1, "normal_x scaled", *normal_x);
    DBG_PRINT_VECTOR(1, "normal_s scaled", *normal_s);

    // Compute the unscaled steps
    normal_s->ElementWiseMultiply(*InexCq().curr_scaling_slacks());
    v_cauchy_s_bak->ElementWiseMultiply(*InexCq().curr_scaling_slacks());

    // We now check if the Dogleg step, shorted by the
    // fraction-to-the-boundary rule, gives at least as much progress
    // as the Cauchy step, also shortened by the
    // fraction-to-the-boundary rule.  If not, we throw away the
    // Newton step component.

    // TODO: Implement efficiently
    const Number tau = IpData().curr_tau();
    Number ftb_cauchy =
      IpCq().primal_frac_to_the_bound(tau, *v_cauchy_x_bak, *v_cauchy_s_bak);
    Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                   "Dogleg: Fraction-to-the-bounary step size for Cauchy step = %23.16e\n", ftb_cauchy);
    inf_c = IpCq().curr_jac_c_times_vec(*v_cauchy_x_bak)->MakeNewCopy();
    inf_c->AddOneVector(1., *curr_c, ftb_cauchy);
    inf_d = curr_d_minus_s->MakeNewCopy();
    inf_d->AddTwoVectors(-ftb_cauchy, *v_cauchy_s_bak,
                         ftb_cauchy, *IpCq().curr_jac_d_times_vec(*v_cauchy_x_bak) , 1.);
    Number objred_ftb_cauchy = 0.5*(IpCq().CalcNormOfType(NORM_2, *curr_c, *curr_d_minus_s)-IpCq().CalcNormOfType(NORM_2, *inf_c, *inf_d));
    Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                   "Dogleg: Reduction of normal problem objective function by ftb cauchy step = %23.16e\n", objred_ftb_cauchy);

    Number ftb_dogleg =
      IpCq().primal_frac_to_the_bound(tau, *normal_x, *normal_s);
    Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                   "Dogleg: Fraction-to-the-bounary step size for Dogleg step = %23.16e\n", ftb_dogleg);
    inf_c = IpCq().curr_jac_c_times_vec(*normal_x)->MakeNewCopy();
    inf_c->AddOneVector(1., *curr_c, ftb_dogleg);
    inf_d = curr_d_minus_s->MakeNewCopy();
    inf_d->AddTwoVectors(-ftb_dogleg, *normal_s,
                         ftb_dogleg, *IpCq().curr_jac_d_times_vec(*normal_x) , 1.);
    Number objred_ftb_dogleg = 0.5*(IpCq().CalcNormOfType(NORM_2, *curr_c, *curr_d_minus_s)-IpCq().CalcNormOfType(NORM_2, *inf_c, *inf_d));
    Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                   "Dogleg: Reduction of normal problem objective function by ftb dogleg step = %23.16e\n", objred_ftb_dogleg);

    Number rhs = 10.*objred_ftb_dogleg;
    Number lhs = objred_ftb_cauchy;
    Number BasVal = curr_c->Nrm2()+curr_d_minus_s->Nrm2();
    bool ok = Compare_le(lhs, rhs, BasVal);
    if (!ok) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM, "Dogleg step: Dogleg step makes less progress than Cauchy step, resetting to Cauchy step.\n");
      normal_x = v_cauchy_x_bak;
      normal_s = v_cauchy_s_bak;
      IpData().Append_info_string("NR ");
    }

    return true;
  }

} // namespace Ipopt
