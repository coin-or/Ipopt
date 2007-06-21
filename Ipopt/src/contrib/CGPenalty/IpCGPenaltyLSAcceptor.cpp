// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpCGPenaltyLSAcceptor.cpp 552 2005-10-27 00:58:58Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//           Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.cpp

#include "IpCGPenaltyLSAcceptor.hpp"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpAlgTypes.hpp"

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

  CGPenaltyLSAcceptor::CGPenaltyLSAcceptor(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      pd_solver_(pd_solver)
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::CGPenaltyLSAcceptor",
                  dbg_verbosity);
  }

  CGPenaltyLSAcceptor::~CGPenaltyLSAcceptor()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::~CGPenaltyLSAcceptor()",
                  dbg_verbosity);
  }

  void CGPenaltyLSAcceptor::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedNumberOption(
      "eta_penalty",
      "Relaxation factor in the Armijo condition for the penalty function.",
      0.0, true, 0.5, true, 1e-8);
    roptions->AddLowerBoundedNumberOption(
      "penalty_update_infeasibility_tol",
      "Threshold for infeasibility in penalty parameter update test.",
      0.0, true, 1e-9,
      "If the new constraint violation is smaller than this tolerance, the "
      "penalty parameter is not increased.");
    roptions->AddLowerBoundedNumberOption(
      "eta_min",
      "LIFENG WRITES THIS.",
      0.0, true, 1e-2,
      "");
    roptions->AddLowerBoundedNumberOption(
      "penalty_update_compl_tol",
      "LIFENG WRITES THIS.",
      0.0, true, 0.2,
      "");
    roptions->AddLowerBoundedNumberOption(
      "chi_hat",
      "LIFENG WRITES THIS.",
      0.0, true, 2.,
      "");
    roptions->AddLowerBoundedNumberOption(
      "chi_tilde",
      "LIFENG WRITES THIS.",
      0.0, true, 5.,
      "");
    roptions->AddLowerBoundedNumberOption(
      "chi_cup",
      "LIFENG WRITES THIS.",
      0.0, true, 1.5,
      "");
    roptions->AddLowerBoundedNumberOption(
      "gamma_hat",
      "LIFENG WRITES THIS.",
      0.0, true, 0.04,
      "");
    roptions->AddLowerBoundedNumberOption(
      "gamma_tilde",
      "LIFENG WRITES THIS.",
      0.0, true, 4.,
      "");
    roptions->AddLowerBoundedNumberOption(
      "penalty_max",
      "LIFENG WRITES THIS.",
      0.0, true, 1e20,
      "");
    roptions->AddLowerBoundedNumberOption(
      "epsilon_c",
      "LIFENG WRITES THIS.",
      0.0, true, 1e-2,
      "");
  }

  bool CGPenaltyLSAcceptor::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("eta_penalty", eta_penalty_, prefix);
    options.GetNumericValue("penalty_update_infeasibility_tol",
                            penalty_update_infeasibility_tol_, prefix);
    options.GetNumericValue("eta_min", eta_min_, prefix);
    options.GetNumericValue("penalty_update_compl_tol",
                            penalty_update_compl_tol_, prefix);
    options.GetNumericValue("chi_hat", chi_hat_, prefix);
    options.GetNumericValue("chi_tilde", chi_tilde_, prefix);
    options.GetNumericValue("chi_cup", chi_cup_, prefix);
    options.GetNumericValue("gamma_hat", gamma_hat_, prefix);
    options.GetNumericValue("gamma_tilde", gamma_tilde_, prefix);
    options.GetNumericValue("penalty_max", penalty_max_, prefix);
    options.GetNumericValue("epsilon_c", epsilon_c_, prefix);
    // The following two option is registered by FilterLSAcceptor
    options.GetIntegerValue("max_soc", max_soc_, prefix);
    if (max_soc_>0) {
      ASSERT_EXCEPTION(IsValid(pd_solver_), OPTION_INVALID,
                       "Option \"max_soc\": This option is non-negative, but no linear solver for computing the SOC given to FilterLSAcceptor object.");
    }
    options.GetNumericValue("kappa_soc", kappa_soc_, prefix);

    Reset();

    counter_penalty_updates_ = 0;
    curr_eta_ = -1.;
    IpData().CGPenData().SetPenaltyUninitialized();

    return true;
  }

  void CGPenaltyLSAcceptor::InitThisLineSearch(bool in_watchdog)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::InitThisLineSearch",
                   dbg_verbosity);

    // Set the values for the reference point
    if (!in_watchdog) {
      reference_penalty_function_ = IpCq().CGPenCq().curr_penalty_function();
      if (IpData().CGPenData().HaveCgFastDeltas()) {
        // use the fast step
        reference_direct_deriv_penalty_function_ =
          IpCq().CGPenCq().curr_fast_direct_deriv_penalty_function();
      }
      else {
        reference_direct_deriv_penalty_function_ =
          IpCq().CGPenCq().curr_direct_deriv_penalty_function();
      }
    }
    else {
      reference_penalty_function_ = watchdog_penalty_function_;
      reference_direct_deriv_penalty_function_ =
        watchdog_direct_deriv_penalty_function_;
    }
  }

  bool
  CGPenaltyLSAcceptor::CheckAcceptabilityOfTrialPoint(Number alpha_primal_test)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::CheckAcceptabilityOfTrialPoint",
                   dbg_verbosity);

    bool accept;

    Number trial_penalty_function = IpCq().CGPenCq().trial_penalty_function();
    DBG_ASSERT(IsFiniteNumber(trial_penalty_function));

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Checking acceptability for trial step size alpha_primal_test=%13.6e:\n", alpha_primal_test);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of penalty function     = %23.16e  (reference %23.16e):\n", trial_penalty_function, reference_penalty_function_);
    if (Jnlst().ProduceOutput(J_DETAILED, J_LINE_SEARCH)) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "curr_barr  = %23.16e curr_inf  = %23.16e\n",
                     IpCq().curr_barrier_obj(),
                     IpCq().curr_primal_infeasibility(NORM_2));
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "trial_barr = %23.16e trial_inf = %23.16e\n",
                     IpCq().trial_barrier_obj(),
                     IpCq().trial_primal_infeasibility(NORM_2));
    }

    // Now check the Armijo condition
    accept = Compare_le(trial_penalty_function-reference_penalty_function_,
                        eta_penalty_*alpha_primal_test*reference_direct_deriv_penalty_function_,
                        reference_penalty_function_);

    return accept;
  }

  Number CGPenaltyLSAcceptor::CalculateAlphaMin()
  {
    // ToDo For now we just return zero
    return 0.;
  }

  bool CGPenaltyLSAcceptor::Compare_le(Number lhs, Number rhs, Number BasVal)
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::Compare_le",
                  dbg_verbosity);
    DBG_PRINT((1,"lhs = %27.16e rhs = %27.16e  BasVal = %27.16e\n",lhs,rhs,BasVal));

    Number mach_eps = std::numeric_limits<Number>::epsilon();
    return (lhs - rhs <= 10.*mach_eps*fabs(BasVal));
  }

  void CGPenaltyLSAcceptor::StartWatchDog()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::StartWatchDog", dbg_verbosity);

    watchdog_penalty_function_ = IpCq().CGPenCq().curr_penalty_function();
    watchdog_direct_deriv_penalty_function_ =
      IpCq().CGPenCq().curr_direct_deriv_penalty_function();
    watchdog_delta_cgpen_ = IpData().CGPenData().delta_cgpen();
  }

  void CGPenaltyLSAcceptor::StopWatchDog()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::StopWatchDog", dbg_verbosity);

    reference_penalty_function_ = watchdog_penalty_function_;
    reference_direct_deriv_penalty_function_ =
      watchdog_direct_deriv_penalty_function_;
    IpData().CGPenData().set_delta_cgpen(watchdog_delta_cgpen_);
    watchdog_delta_cgpen_ = NULL;
  }

  void CGPenaltyLSAcceptor::Reset()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::Reset", dbg_verbosity);
    /*
    counter_penalty_updates_ = 0;
    curr_eta_ = -1.;
    IpData().CGPenData().SetPenaltyUninitialized();
    */
  }

  bool
  CGPenaltyLSAcceptor::TrySecondOrderCorrection(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::TrySecondOrderCorrection",
                   dbg_verbosity);

    if (max_soc_==0) {
      return false;
    }

    bool accept = false;
    Index count_soc = 0;

    Number theta_soc_old = 0.;
    Number theta_trial = IpCq().trial_constraint_violation();
    Number alpha_primal_soc = alpha_primal;

    // delta_y_c and delta_y_d are the steps used in the right hand
    // side for the SOC step
    SmartPtr<const Vector> delta_y_c = IpData().delta()->y_c();
    SmartPtr<const Vector> delta_y_d = IpData().delta()->y_d();

    SmartPtr<Vector> c_soc = IpCq().curr_c()->MakeNewCopy();
    SmartPtr<Vector> dms_soc = IpCq().curr_d_minus_s()->MakeNewCopy();

    while (count_soc<max_soc_ && !accept &&
           (count_soc==0 || theta_trial<=kappa_soc_*theta_soc_old) ) {
      theta_soc_old = theta_trial;

      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Trying second order correction number %d\n",
                     count_soc+1);

      // Compute SOC constraint violation
      Number c_over_r = IpCq().CGPenCq().curr_cg_pert_fact();
      c_soc->AddTwoVectors(1.0, *IpCq().trial_c(),
                           -c_over_r, *delta_y_c,
                           alpha_primal_soc);
      dms_soc->AddTwoVectors(1.0, *IpCq().trial_d_minus_s(),
                             -c_over_r, *delta_y_d,
                             alpha_primal_soc);

      // Compute the SOC search direction
      SmartPtr<IteratesVector> delta_soc =
        actual_delta->MakeNewIteratesVector(true);
      SmartPtr<IteratesVector> rhs = actual_delta->MakeNewContainer();
      rhs->Set_x(*IpCq().curr_grad_lag_with_damping_x());
      rhs->Set_s(*IpCq().curr_grad_lag_with_damping_s());
      rhs->Set_y_c(*c_soc);
      rhs->Set_y_d(*dms_soc);
      rhs->Set_z_L(*IpCq().curr_relaxed_compl_x_L());
      rhs->Set_z_U(*IpCq().curr_relaxed_compl_x_U());
      rhs->Set_v_L(*IpCq().curr_relaxed_compl_s_L());
      rhs->Set_v_U(*IpCq().curr_relaxed_compl_s_U());
      pd_solver_->Solve(-1.0, 0.0, *rhs, *delta_soc, true);

      // Update the delta_y_c and delta_y_d vectors in case we do
      // additional SOC steps
      delta_y_c = ConstPtr(delta_soc->y_c());
      delta_y_d = ConstPtr(delta_soc->y_d());

      // Compute step size
      alpha_primal_soc =
        IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                        *delta_soc->x(),
                                        *delta_soc->s());

      // Check if trial point is acceptable
      try {
        // Compute the primal trial point
        IpData().SetTrialPrimalVariablesFromStep(alpha_primal_soc, *delta_soc->x(), *delta_soc->s());

        // in acceptance tests, use original step size!
        accept = CheckAcceptabilityOfTrialPoint(alpha_primal_test);
      }
      catch(IpoptNLP::Eval_Error& e) {
        e.ReportException(Jnlst(), J_DETAILED);
        Jnlst().Printf(J_WARNING, J_MAIN, "Warning: SOC step rejected due to evaluation error\n");
        IpData().Append_info_string("e");
        accept = false;
        // There is no point in continuing SOC procedure
        break;
      }

      if (accept) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Second order correction step accepted with %d corrections.\n", count_soc+1);
        // Accept all SOC quantities
        alpha_primal = alpha_primal_soc;
        actual_delta = delta_soc;
      }
      else {
        count_soc++;
        theta_trial = IpCq().trial_constraint_violation();
      }
    }
    return accept;
  }

  bool
  CGPenaltyLSAcceptor::TryCorrector(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::TryCorrector",
                   dbg_verbosity);

    return false;
  }

  char CGPenaltyLSAcceptor::UpdateForNextIteration(Number alpha_primal_test)
  {
    char info_alpha_primal_char='?';

    if (curr_eta_<0.) {
      // We need to initialize the eta tolerance
      curr_eta_ = Max(eta_min_, Min(gamma_tilde_,
                                    gamma_hat_*IpCq().curr_nlp_error()));
    }

    // Check if the penalty parameter is to be increased
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "Starting tests for penalty parameter update:\n");

    // We use the new infeasibility here...
    Number trial_inf = IpCq().trial_primal_infeasibility(NORM_2);
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "trial infeasibility = %8.2\n", trial_inf);
    bool increase = (trial_inf >= penalty_update_infeasibility_tol_);
    if (!increase) {
      info_alpha_primal_char='i';
    }

    if (increase) {
      Number max_step = Max(IpData().CGPenData().delta_cgpen()->x()->Amax(),
                            IpData().CGPenData().delta_cgpen()->s()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Max norm of step = %8.2\n", max_step);
      increase = (max_step <= curr_eta_);
      if (!increase) {
        info_alpha_primal_char='d';
      }
    }

    // Lifeng: Should we use the new complementarity here?  If so, I
    // have to restructure BacktrackingLineSearch
    Number mu = IpData().curr_mu();
    if (increase) {
      Number min_compl = mu;
      Number max_compl = mu;
      if (IpNLP().x_L()->Dim()>0) {
        SmartPtr<const Vector> compl_x_L = IpCq().curr_compl_x_L();
        min_compl = Min(min_compl, compl_x_L->Min());
        max_compl = Max(max_compl, compl_x_L->Max());
      }
      if (IpNLP().x_U()->Dim()>0) {
        SmartPtr<const Vector> compl_x_U = IpCq().curr_compl_x_U();
        min_compl = Min(min_compl, compl_x_U->Min());
        max_compl = Max(max_compl, compl_x_U->Max());
      }
      if (IpNLP().d_L()->Dim()>0) {
        SmartPtr<const Vector> compl_s_L = IpCq().curr_compl_s_L();
        min_compl = Min(min_compl, compl_s_L->Min());
        max_compl = Max(max_compl, compl_s_L->Max());
      }
      if (IpNLP().d_U()->Dim()>0) {
        SmartPtr<const Vector> compl_s_U = IpCq().curr_compl_s_U();
        min_compl = Min(min_compl, compl_s_U->Min());
        max_compl = Max(max_compl, compl_s_U->Max());
      }
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Minimal compl = %8.2\n", min_compl);
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Maximal compl = %8.2\n", max_compl);
      increase = (min_compl >= mu*penalty_update_compl_tol_ &&
                  max_compl <= mu/penalty_update_compl_tol_);
      if (!increase) {
        info_alpha_primal_char='c';
      }
    }

    // Lifeng: Here I'm using the information from the current step
    // and the current infeasibility
    if (increase) {
      SmartPtr<Vector> vec = IpData().curr()->y_c()->MakeNewCopy();
      vec->AddTwoVectors(1., *IpData().CGPenData().delta_cgpen()->y_c(),
                         -1./IpCq().CGPenCq().curr_cg_pert_fact(), *IpCq().curr_c(),
                         1.);
      Number omega_test = vec->Amax();
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "omega_test for c = %8.2\n", omega_test);
      increase = (omega_test < curr_eta_);
      if (increase) {
        SmartPtr<Vector> vec = IpData().curr()->y_d()->MakeNewCopy();
        vec->AddTwoVectors(1., *IpData().delta()->y_d(),
                           -1./IpCq().CGPenCq().curr_cg_pert_fact(), *IpCq().curr_d_minus_s(),
                           1.);
        omega_test = vec->Amax();
        Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                       "omega_test for d = %8.2\n", omega_test);
        increase = (omega_test < curr_eta_);
      }
      if (!increase) {
        info_alpha_primal_char='m';
      }
    }

    if (increase) {
      // Ok, now we should increase the penalty parameter
      counter_penalty_updates_++;

      // Update the eta tolerance
      curr_eta_ = Max(eta_min_, curr_eta_/2.);
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Updating eta to = %8.2\n", curr_eta_);
      Number penalty = IpData().CGPenData().curr_penalty();
      Number y_full_step_max;
      SmartPtr<Vector> vec = IpData().curr()->y_c()->MakeNew();
      vec->AddTwoVectors(1., *IpData().curr()->y_c(),
                         1., *IpData().CGPenData().delta_cgpen()->y_c(), 0.);
      y_full_step_max = vec->Amax();
      vec = IpData().curr()->y_d()->MakeNew();
      vec->AddTwoVectors(1., *IpData().curr()->y_d(),
                         1., *IpData().CGPenData().delta_cgpen()->y_d(), 0.);
      y_full_step_max = Max(y_full_step_max, vec->Amax());
      if (IpCq().curr_primal_infeasibility(NORM_2) >= epsilon_c_) {
        penalty = Max(chi_hat_*penalty, y_full_step_max + 1.);
        info_alpha_primal_char = 'l';
      }
      else {
        penalty = Max(chi_tilde_*penalty, chi_cup_*y_full_step_max);
        info_alpha_primal_char = 's';
      }
      if (penalty > penalty_max_) {
        THROW_EXCEPTION(IpoptException, "Penalty parameter becomes too large.");
      }
      IpData().CGPenData().Set_penalty(penalty);

      char spen[40];
      sprintf(spen, " penalty=%8.2e", penalty);
      IpData().Append_info_string(spen);
    }

    return info_alpha_primal_char;
  }

  void CGPenaltyLSAcceptor::PrepareRestoPhaseStart()
  {
    DBG_ASSERT(false && "PrepareRestoPhaseStart not yet implemented");
  }


} // namespace Ipopt
