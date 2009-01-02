// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2008-09-11
//               derived file from IpPenaltyLSAcceptor.cpp (rev 1121)

#include "IpInexactLSAcceptor.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

// for sprintf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  InexactLSAcceptor::InexactLSAcceptor()
  {
    DBG_START_FUN("InexactLSAcceptor::InexactLSAcceptor",
                  dbg_verbosity);
  }

  InexactLSAcceptor::~InexactLSAcceptor()
  {
    DBG_START_FUN("InexactLSAcceptor::~InexactLSAcceptor()",
                  dbg_verbosity);
  }

  void InexactLSAcceptor::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "nu_update_inf_skip_tol",
      "Lower bound on infeasibility to perform penalty parameter update.",
      0.0, true, 1e-9,
      "If the current infeasibility is less than this value, the penalty "
      "parameter update is skipped");
    roptions->AddStringOption2(
      "flexible_penalty_function",
      "yes",
      "no", "do not use the flexible penalty function procedure",
      "yes", "use the flexible penalty function procedure",
      "This determines if the flexible penalty function procedure by "
      "Curtis/Nocedal should be used in the line search.  For now, this only "
      "is implemented for the inexact algorithm.");
    roptions->AddLowerBoundedNumberOption(
      "nu_low_init",
      "Initial value for the lower penalty parameter.",
      0.0, true, 1e-6,
      "This is the initial value for the lower penalty parameter in the "
      "Curtis/Nocedal flexible penalty function line search procedure.  This "
      "must be smaller or equal to the intial value of the upper penalty "
      "parameter, see option \"nu_init\".");
  }

  bool InexactLSAcceptor::InitializeImpl(const OptionsList& options,
                                         const std::string& prefix)
  {
    options.GetNumericValue("nu_init", nu_init_, prefix);
    options.GetNumericValue("nu_inc", nu_inc_, prefix);
    options.GetNumericValue("eta_phi", eta_, prefix);
    options.GetNumericValue("rho", rho_, prefix);
    options.GetNumericValue("tcc_theta", tcc_theta_, prefix);
    options.GetNumericValue("nu_update_inf_skip_tol", nu_update_inf_skip_tol_,
                            prefix);
    options.GetBoolValue("flexible_penalty_function",
                         flexible_penalty_function_, prefix);
    if (flexible_penalty_function_) {
      options.GetNumericValue("nu_low_init", nu_low_init_, prefix);
      ASSERT_EXCEPTION(nu_low_init_<=nu_init_, OPTION_INVALID,
                       "Option \"nu_low_init\" must be smaller or equal to \"nu_init\"");
    }

    // The following options have been declared in FilterLSAcceptor
    Index max_soc;
    options.GetIntegerValue("max_soc", max_soc, prefix);
    ASSERT_EXCEPTION(max_soc==0, OPTION_INVALID,
                     "Option \"max_soc\" must be zero for inexact version.");

    Reset();

    return true;
  }

  void InexactLSAcceptor::InitThisLineSearch(bool in_watchdog)
  {
    DBG_START_METH("InexactLSAcceptor::InitThisLineSearch",
                   dbg_verbosity);

    InexData().set_full_step_accepted(false);

    // Set the values for the reference point
    if (!in_watchdog) {
      reference_theta_ = IpCq().curr_constraint_violation();
      reference_barr_ = IpCq().curr_barrier_obj();

      //////////////////// Update the penalty parameter

      const Number uWu = InexCq().curr_uWu();
      SmartPtr<const Vector> tangential_x = InexData().tangential_x();
      SmartPtr<const Vector> tangential_s = InexData().tangential_s();
      const Number scaled_tangential_norm =
        InexCq().slack_scaled_norm(*tangential_x, *tangential_s);

      SmartPtr<const Vector> delta_x = IpData().delta()->x();
      SmartPtr<const Vector> delta_s = IpData().delta()->s();

      // Get the product of the steps with the Jacobian
      SmartPtr<Vector> cplusAd_c = IpData().curr()->y_c()->MakeNew();
      cplusAd_c->AddTwoVectors(1., *IpCq().curr_jac_c_times_vec(*delta_x),
                               1., *IpCq().curr_c(), 0.);
      SmartPtr<Vector> cplusAd_d = delta_s->MakeNew();
      cplusAd_d->AddTwoVectors(1., *IpCq().curr_jac_d_times_vec(*delta_x),
                               -1., *delta_s,  0.);
      cplusAd_d->Axpy(1., *IpCq().curr_d_minus_s());
      const Number norm_cplusAd =
        IpCq().CalcNormOfType(NORM_2, *cplusAd_c, *cplusAd_d);

      const Number gradBarrTDelta = IpCq().curr_gradBarrTDelta();

      DBG_PRINT((1,"gradBarrTDelta = %e norm_cplusAd = %e reference_theta_ = %e\n", gradBarrTDelta, norm_cplusAd, reference_theta_));

      // update the upper penalty parameter
      Number nu_mid = nu_;
      Number norm_delta_xs = Max(delta_x->Amax(), delta_s->Amax());
      last_nu_ = nu_;
      in_tt2_ = false;
      if (norm_delta_xs == 0.) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "  Zero step, skipping line search\n");
        in_tt2_ = true;
      }
      // TODO: We need to find a proper cut-off value
      else if (reference_theta_ > nu_update_inf_skip_tol_) {
        DBG_PRINT((1,"uWu=%e scaled_tangential_norm=%e\n",uWu ,scaled_tangential_norm ));
        Number numerator = (gradBarrTDelta + Max(0.5*uWu, tcc_theta_*pow(scaled_tangential_norm,2)));
        Number denominator = (1-rho_)*(reference_theta_-norm_cplusAd);
        const Number nu_trial = numerator/denominator;
//DELETEME
        char snu[64];
        sprintf(snu, " nt=%8.2e", nu_trial);
        IpData().Append_info_string(snu);
        Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                       "In penalty parameter update formula:\n  gradBarrTDelta = %e 0.5*uWu = %e tcc_theta_*pow(scaled_tangential_norm,2) = %e numerator = %e\n  reference_theta_ = %e norm_cplusAd + %e denominator = %e nu_trial = %e\n", gradBarrTDelta, 0.5*uWu, tcc_theta_*pow(scaled_tangential_norm,2), numerator, reference_theta_, norm_cplusAd, denominator, nu_trial);

        if (nu_ < nu_trial) {
          nu_ = nu_trial + nu_inc_;
        }
        if (flexible_penalty_function_) {
          last_nu_low_ = nu_low_;
          nu_mid = Max(nu_low_, nu_trial);
          Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                         "nu_low = %8.2e\n", nu_low_);
        }
#if 0
// DELETEME
        if (nu_trial < 0.) {
          nu_ = 1e-6;
        }
#endif
      }
      else {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Warning: Skipping nu update because current constraint violation (%e) less than nu_update_inf_skip_tol.\n", reference_theta_);
        IpData().Append_info_string("nS");
      }
      InexData().set_curr_nu(nu_);
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "  using nu = %23.16e (nu_mid = %23.16e)\n", nu_, nu_mid);

      // Compute the linear model reduction prediction
      DBG_PRINT((1,"gradBarrTDelta=%e reference_theta_=%e norm_cplusAd=%e\n", gradBarrTDelta, reference_theta_, norm_cplusAd));
      reference_pred_ = gradBarrTDelta - nu_mid*(reference_theta_ - norm_cplusAd);

      watchdog_pred_ = -1e300;
    }
    else {
      reference_theta_ = watchdog_theta_;
      reference_barr_ = watchdog_barr_;
      reference_pred_ = watchdog_pred_;
    }
  }

  Number
  InexactLSAcceptor::CalcPred(Number alpha)
  {
    DBG_START_METH("InexactLSAcceptor::CalcPred",
                   dbg_verbosity);

    Number pred = alpha*reference_pred_;

    if (pred > 0.) {
      Jnlst().Printf(J_WARNING, J_LINE_SEARCH, "  pred = %23.16e is positive.  Setting to zero.\n", pred);
      pred = 0.;
    }

    return pred;
  }

  bool
  InexactLSAcceptor::CheckAcceptabilityOfTrialPoint(Number alpha_primal_test)
  {
    DBG_START_METH("InexactLSAcceptor::CheckAcceptabilityOfTrialPoint",
                   dbg_verbosity);

    // If we are in termiation test 2 iteration, we skip the line search
    if (in_tt2_) {
      return true;
    }

    // First compute the barrier function and constraint violation at the
    // current iterate and the trial point

    Number trial_theta = IpCq().trial_constraint_violation();
    Number trial_barr = IpCq().trial_barrier_obj();
    DBG_ASSERT(IsFiniteNumber(trial_barr));

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Checking acceptability for trial step size alpha_primal_test=%13.6e:\n", alpha_primal_test);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of barrier function     = %23.16e  (reference %23.16e):\n", trial_barr, reference_barr_);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of constraint violation = %23.16e  (reference %23.16e):\n", trial_theta, reference_theta_);

    Number pred = CalcPred(alpha_primal_test);
    resto_pred_ = pred;
    DBG_PRINT((1, "nu_ = %e reference_barr_ + nu_*(reference_theta_)=%e trial_barr + nu_*trial_theta=%e\n",nu_,reference_barr_ + nu_*(reference_theta_),trial_barr + nu_*trial_theta));
    Number ared = reference_barr_ + nu_*(reference_theta_) -
                  (trial_barr + nu_*trial_theta);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Checking Armijo Condition with pred = %23.16e and ared = %23.16e\n", pred, ared);

    bool accept =
      Compare_le(eta_*pred, ared, reference_barr_ + nu_*(reference_theta_));
    bool accept_low = false;
    if (flexible_penalty_function_) {
      accept_low = Compare_le(eta_*pred, ared, reference_barr_ + nu_low_*(reference_theta_));
    }

    if (accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "   Success...\n");
    }
    else if (flexible_penalty_function_ && accept_low) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "   Success with nu_low...\n");
      accept = true;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "   Failed...\n");
    }
    if (accept) {
      // HERE WE RESET THE SLACKS.  MAYBE THIS SHOULD BE BEFORE THE
      // FUNCTION EVALUATIONS?
      ResetSlacks();

      if (flexible_penalty_function_) {
        // update the lower penalty parameter if necessary
        if (!accept_low) {
          Number nu_real = -(trial_barr - reference_barr_)/(trial_theta - reference_theta_);
          nu_low_ = Min(nu_, nu_low_ + Max(0.1*(nu_real-nu_low_), nu_inc_));

          Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                         "Updating nu_low to %8.2e with nu_real = %8.2e\n", nu_low_, nu_real);
        }
      }
    }
    return accept;
  }

  void InexactLSAcceptor::ResetSlacks()
  {
    DBG_START_METH("InexactLSAcceptor::ResetSlacks",
                   dbg_verbosity);

    DBG_PRINT_VECTOR(1, "sorig", *IpData().trial()->s());
    DBG_PRINT_VECTOR(1, "dtrial", *IpCq().trial_d());
    SmartPtr<Vector> new_s = IpData().trial()->s()->MakeNew();
    SmartPtr<Vector> tmp_d = IpNLP().d_L()->MakeNew();
    IpNLP().Pd_L()->TransMultVector(1., *IpCq().trial_d(), 0., *tmp_d);
    SmartPtr<Vector> tmp_s = IpNLP().d_L()->MakeNew();
    IpNLP().Pd_L()->TransMultVector(1., *IpData().trial()->s(), 0., *tmp_s);
    tmp_s->ElementWiseMax(*tmp_d);
    IpNLP().Pd_L()->MultVector(1., *tmp_s, 0., *new_s);
    tmp_d = IpNLP().d_U()->MakeNew();
    IpNLP().Pd_U()->TransMultVector(1., *IpCq().trial_d(), 0., *tmp_d);
    tmp_s = IpNLP().d_U()->MakeNew();
    IpNLP().Pd_U()->TransMultVector(1., *IpData().trial()->s(), 0., *tmp_s);
    tmp_s->ElementWiseMin(*tmp_d);
    IpNLP().Pd_U()->MultVector(1., *tmp_s, 1., *new_s);
    SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
    trial->Set_s(*new_s);
    IpData().set_trial(trial);
    DBG_PRINT_VECTOR(1, "new_s", *IpData().trial()->s());
  }

  Number InexactLSAcceptor::CalculateAlphaMin()
  {
    // ToDo: make better
    return 1e-16;
  }

  void InexactLSAcceptor::StartWatchDog()
  {
    DBG_START_FUN("InexactLSAcceptor::StartWatchDog", dbg_verbosity);

    watchdog_theta_ = reference_theta_;
    watchdog_barr_ = reference_barr_;
    watchdog_pred_ = reference_pred_;
  }

  void InexactLSAcceptor::StopWatchDog()
  {
    DBG_START_FUN("InexactLSAcceptor::StopWatchDog", dbg_verbosity);

    reference_theta_ = watchdog_theta_;
    reference_barr_ = watchdog_barr_;
    reference_pred_ = watchdog_pred_;
  }

  void InexactLSAcceptor::Reset()
  {
    DBG_START_FUN("InexactLSAcceptor::Reset", dbg_verbosity);

    nu_ = nu_init_;
    if (flexible_penalty_function_) {
      nu_low_ = nu_low_init_;
    }
    InexData().set_curr_nu(nu_);
  }

  bool
  InexactLSAcceptor::TrySecondOrderCorrection(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_METH("InexactLSAcceptor::TrySecondOrderCorrection",
                   dbg_verbosity);
    return false;
  }

  bool
  InexactLSAcceptor::TryCorrector(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    return false;
  }

  char InexactLSAcceptor::UpdateForNextIteration(Number alpha_primal_test)
  {
    char info_alpha_primal_char = 'k';
    // Augment the filter if required
    if (last_nu_ != nu_) {
      info_alpha_primal_char = 'n';
      char snu[40];
      sprintf(snu, " nu=%8.2e", nu_);
      IpData().Append_info_string(snu);
    }
    if (last_nu_low_ != nu_low_) {
      char snu[40];
      sprintf(snu, " nl=%8.2e", nu_low_);
      IpData().Append_info_string(snu);
      if (info_alpha_primal_char == 'k' ) {
        info_alpha_primal_char = 'l';
      }
      else {
        info_alpha_primal_char = 'b';
      }
    }

    if (alpha_primal_test==1. && watchdog_pred_==-1e300) {
      InexData().set_full_step_accepted(true);

    }
    return info_alpha_primal_char;
  }

  void InexactLSAcceptor::PrepareRestoPhaseStart()
  {
    THROW_EXCEPTION(INTERNAL_ABORT, "Restoration phase called");
  }

  bool
  InexactLSAcceptor::IsAcceptableToCurrentIterate(Number trial_barr,
      Number trial_theta,
      bool called_from_restoration /*=false*/) const
  {
    DBG_START_METH("InexactLSAcceptor::IsAcceptableToCurrentIterate",
                   dbg_verbosity);
    THROW_EXCEPTION(INTERNAL_ABORT, "InexactLSAcceptor::IsAcceptableToCurrentIterate called");
    ASSERT_EXCEPTION(resto_pred_ >= 0., INTERNAL_ABORT,
                     "resto_pred_ not set for check from restoration phase.");

    Number ared = reference_barr_ + nu_*(reference_theta_) -
                  (trial_barr + nu_*trial_theta);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Checking Armijo Condition (for resto) with pred = %23.16e and ared = %23.16e\n",
                   resto_pred_, ared);

    bool accept;
    if (Compare_le(eta_*resto_pred_, ared, reference_barr_ + nu_*(reference_theta_))) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "   Success...\n");
      accept = true;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "   Failed...\n");
      accept = false;
    }
    return accept;
  }

} // namespace Ipopt
