// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//           Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.cpp

#include "IpCGPenaltyLSAcceptor.hpp"
#include "IpCGPenaltyData.hpp"
#include "IpCGPenaltyCq.hpp"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpAlgTypes.hpp"
#include "IpIpoptAlg.hpp"

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
      PiecewisePenalty_(1),
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
    roptions->AddStringOption2(
      "never_use_piecewise_penalty_ls",
      "Toggle to switch off the piecewise penalty method",
      "no",
      "no", "always use the piecewise penalty method",
      "yes", "never use the piecewise penalty method",
      "");
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
      0.0, true, 1e1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "pen_theta_max_fact",
      "Determines upper bound for constraint violation in the filter.",
      0.0, true, 1e4,
      "The algorithmic parameter theta_max is determined as theta_max_fact "
      "times the maximum of 1 and the constraint violation at initial point.  "
      "Any point with a constraint violation larger than theta_max is "
      "unacceptable to the filter (see Eqn. (21) in implementation paper).");
    roptions->AddLowerBoundedNumberOption(
      "penalty_update_compl_tol",
      "LIFENG WRITES THIS.",
      0.0, true, 1e1,
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
      "epsilon_c",
      "LIFENG WRITES THIS.",
      0.0, true, 1e-2,
      "");
    roptions->AddLowerBoundedNumberOption(
      "piecewisepenalty_gamma_obj",
      "LIFENG WRITES THIS.",
      0.0, true, 1e-13,
      "");
    roptions->AddLowerBoundedNumberOption(
      "piecewisepenalty_gamma_infeasi",
      "LIFENG WRITES THIS.",
      0.0, true, 1e-13,
      "");
    roptions->AddLowerBoundedNumberOption(
      "min_alpha_primal",
      "LIFENG WRITES THIS.",
      0.0, true, 1e-13,
      "");
    roptions->AddLowerBoundedNumberOption(
      "theta_min",
      "LIFENG WRITES THIS.",
      0.0, true, 1e-6,
      "");
    roptions->AddLowerBoundedNumberOption(
      "mult_diverg_feasibility_tol",
      "tolerance for deciding if the multipliers are diverging",
      0, true, 1e-7,
      "");
    roptions->AddLowerBoundedNumberOption(
      "mult_diverg_y_tol",
      "tolerance for deciding if the multipliers are diverging",
      0, true, 1e8,
      "");

  }

  bool CGPenaltyLSAcceptor::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetBoolValue("never_use_piecewise_penalty_ls",
                         never_use_piecewise_penalty_ls_, prefix);
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
    options.GetNumericValue("epsilon_c", epsilon_c_, prefix);
    options.GetNumericValue("piecewisepenalty_gamma_obj",
                            piecewisepenalty_gamma_obj_, prefix);
    options.GetNumericValue("piecewisepenalty_gamma_infeasi",
                            piecewisepenalty_gamma_infeasi_, prefix);
    options.GetNumericValue("pen_theta_max_fact", pen_theta_max_fact_, prefix);
    options.GetNumericValue("min_alpha_primal", min_alpha_primal_, prefix);
    options.GetNumericValue("theta_min", theta_min_, prefix);
    options.GetNumericValue("mult_diverg_feasibility_tol", mult_diverg_feasibility_tol_, prefix);
    options.GetNumericValue("mult_diverg_y_tol", mult_diverg_y_tol_, prefix);
    // The following option has been registered by FilterLSAcceptor
    options.GetIntegerValue("max_soc", max_soc_, prefix);
    // The following option has been registered by CGSearhDirCalc
    options.GetNumericValue("penalty_max", penalty_max_, prefix);

    if (max_soc_>0) {
      ASSERT_EXCEPTION(IsValid(pd_solver_), OPTION_INVALID,
                       "Option \"max_soc\": This option is non-negative, but no linear solver for computing the SOC given to FilterLSAcceptor object.");
    }
    options.GetNumericValue("kappa_soc", kappa_soc_, prefix);

    pen_theta_max_ = -1.;
    pen_curr_mu_ = IpData().curr_mu();
    counter_first_type_penalty_updates_ = 0;
    counter_second_type_penalty_updates_ = 0;
    curr_eta_ = -1.;
    CGPenData().SetPenaltyUninitialized();
    ls_counter_ = 0;
    best_KKT_error_ = -1.;
    accepted_by_Armijo_ = true;
    //never_do_restor_ = true;
    jump_for_tiny_step_ = 0;

    return true;
  }

  void CGPenaltyLSAcceptor::InitThisLineSearch(bool in_watchdog)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::InitThisLineSearch",
                   dbg_verbosity);

    accepted_by_Armijo_ = true;
    ls_counter_ = 0;

    // If the algorithm restarts from a previous iteration, reset line search parameters.
    if (CGPenData().restor_iter() == IpData().iter_count()) {
      Reset();
    }
    // Every time mu is decreased, reset line search parameters.
    if (pen_curr_mu_ > IpData().curr_mu()) {
      Reset();
    }
    if (reset_piecewise_penalty_) {
      Number curr_barr = IpCq().curr_barrier_obj();
      Number curr_infeasi =  IpCq().curr_constraint_violation();
      PiecewisePenalty_.InitPiecewisePenaltyList(0.,curr_barr, curr_infeasi);
      reset_piecewise_penalty_ = false;
    }
    // Set the values for the reference point
    if (!in_watchdog) {
      reference_penalty_function_ = CGPenCq().curr_penalty_function();
      reference_theta_ = IpCq().curr_constraint_violation();
      if (CGPenData().HaveCgFastDeltas()) {
        // use the fast step
        reference_direct_deriv_penalty_function_ =
          CGPenCq().curr_fast_direct_deriv_penalty_function();
      }
      else {
        reference_direct_deriv_penalty_function_ =
          CGPenCq().curr_direct_deriv_penalty_function();
      }
    }
    else {
      //reference_theta_ = watchdog_theta_;
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

    Number curr_barr = IpCq().curr_barrier_obj();
    Number curr_infeasi =  IpCq().curr_constraint_violation();
    //Number trial_barr = IpCq().trial_barrier_obj();
    Number trial_infeasi =  IpCq().trial_constraint_violation();
    bool accept = false;
    ls_counter_++;
    if (ls_counter_ == 1) {
      CGPenData().SetPrimalStepSize(alpha_primal_test);
    }
    if (jump_for_tiny_step_ == 1) {
      jump_for_tiny_step_ = 0;
      Reset();
      IpData().Append_info_string("jump");
      return true;
    }
    /*
    if (jump_for_tiny_step_ == 1) {
      
      PiecewisePenalty_.InitPiecewisePenaltyList(0.,trial_barr, trial_infeasi);    
      jump_for_tiny_step_ = 0;
    }
    */
    // Initialize the piecewise penalty list in case that it's empty.
    if ( PiecewisePenalty_.IsPiecewisePenaltyListEmpty() ) {
      PiecewisePenalty_.InitPiecewisePenaltyList(0.,curr_barr, curr_infeasi);
    }
    // Initialize the max infeasibility that's allowed for every iteration,
    if (pen_theta_max_ < 0.) {
      pen_theta_max_ = pen_theta_max_fact_*Max(1.0, reference_theta_);
    }
    // Check if the constraint violation is becoming too large. If the violation
    // is bigger than pen_theta_max_, the trial point is rejected.
    if (pen_theta_max_>0 &&
        trial_infeasi > pen_theta_max_) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "trial_infeasi = %e is larger than theta_max = %e\n",
                     trial_infeasi, pen_theta_max_);
      return false;
    }
    // Check Armijo conditions.
    if (!accept) {
      accept = ArmijoHolds(alpha_primal_test);
    }
    // Check PLPF criteria.
    if (!accept && !never_use_piecewise_penalty_ls_) {
      accept = IsAcceptableToPiecewisePenalty(alpha_primal_test);
      if (accept) {
        accepted_by_Armijo_ = false;
      }
    }


    if (alpha_primal_test < min_alpha_primal_ && !accept) {
      accept = true;
    }
    if (accept) {
      if (ls_counter_ > 15 && alpha_primal_test < 1e-5 && jump_for_tiny_step_ == 0) {
        jump_for_tiny_step_ = 1;
      }
      ls_counter_ = 0;
    }


    return accept;

  }


  Number CGPenaltyLSAcceptor::CalculateAlphaMin()
  {
    // ToDo For now we just return zero
    Number min_step_size = 0.;
    return min_step_size;
  }

  bool CGPenaltyLSAcceptor::IsAcceptableToPiecewisePenalty(Number alpha_primal_test)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::IsAcceptableToPiecewisePenalty",
                   dbg_verbosity);

    // If the current infeasibility is small, we require the barrier to be decreased.
    bool accept = false;
    Number infeasibility = IpCq().curr_primal_infeasibility(NORM_MAX);
    SmartPtr<const Vector> dx = IpData().delta()->x();
    SmartPtr<const Vector> ds = IpData().delta()->s();
    Number curr_barr = IpCq().curr_barrier_obj();
    Number trial_barr = IpCq().trial_barrier_obj();
    Number nrm_dx_ds = pow(dx->Nrm2(),2.) + pow(ds->Nrm2(),2.);
    if (infeasibility < theta_min_) {
      Number biggest_barr = PiecewisePenalty_.BiggestBarr();
      accept = Compare_le(trial_barr-biggest_barr, -alpha_primal_test*
                          piecewisepenalty_gamma_obj_*nrm_dx_ds, curr_barr);
      if (!accept) {
        return false;
      }
    }
    Number Fzconst, Fzlin;
    Fzconst = IpCq().trial_barrier_obj() + alpha_primal_test *
              piecewisepenalty_gamma_obj_ * nrm_dx_ds;
    Fzlin = IpCq().trial_constraint_violation() + alpha_primal_test *
            piecewisepenalty_gamma_infeasi_ * IpCq().curr_constraint_violation();
    accept  = PiecewisePenalty_.Acceptable(Fzconst, Fzlin);
    return accept;
  }

  bool CGPenaltyLSAcceptor::ArmijoHolds(Number alpha_primal_test)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::ArmijoHolds",
                   dbg_verbosity);
    bool accept = false;
    Number trial_penalty_function = CGPenCq().trial_penalty_function();
    DBG_ASSERT(IsFiniteNumber(trial_penalty_function));
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Checking acceptability for trial step size alpha_primal_test=%13.6e:\n", alpha_primal_test);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   " New values of penalty function     = %23.16e  (reference %23.16e):\n", trial_penalty_function, reference_penalty_function_);
    if (Jnlst().ProduceOutput(J_DETAILED, J_LINE_SEARCH)) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "curr_barr  = %23.16e curr_inf  = %23.16e\n",
                     IpCq().curr_barrier_obj(),
                     IpCq().curr_constraint_violation());
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "trial_barr = %23.16e trial_inf = %23.16e\n",
                     IpCq().trial_barrier_obj(),
                     IpCq().trial_constraint_violation());
    }
    // Now check the Armijo condition
    accept = Compare_le(trial_penalty_function-reference_penalty_function_,
                        eta_penalty_*alpha_primal_test*reference_direct_deriv_penalty_function_,
                        reference_penalty_function_);
    return accept;
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
    watchdog_penalty_function_ = CGPenCq().curr_penalty_function();
    watchdog_direct_deriv_penalty_function_ =
      CGPenCq().curr_direct_deriv_penalty_function();
    watchdog_delta_cgpen_ = CGPenData().delta_cgpen();
  }

  void CGPenaltyLSAcceptor::StopWatchDog()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::StopWatchDog", dbg_verbosity);
    reference_penalty_function_ = watchdog_penalty_function_;
    reference_direct_deriv_penalty_function_ =
      watchdog_direct_deriv_penalty_function_;
    CGPenData().set_delta_cgpen(watchdog_delta_cgpen_);
    watchdog_delta_cgpen_ = NULL;
  }

  void CGPenaltyLSAcceptor::Reset()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::Reset", dbg_verbosity);
    reset_piecewise_penalty_ = true;
    pen_theta_max_ = -1.;
    curr_eta_ = -1.;
    pen_curr_mu_ = IpData().curr_mu();
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
    Number theta_soc_old2 = 0.;
    Number theta_trial = IpCq().trial_constraint_violation();
    Number theta_trial2 = IpCq().curr_primal_infeasibility(NORM_2);
    Number alpha_primal_soc = alpha_primal;
    // delta_y_c and delta_y_d are the steps used in the right hand
    // side for the SOC step
    SmartPtr<const Vector> delta_y_c = IpData().delta()->y_c();
    SmartPtr<const Vector> delta_y_d = IpData().delta()->y_d();
    SmartPtr<Vector> c_soc = IpCq().curr_c()->MakeNewCopy();
    SmartPtr<Vector> dms_soc = IpCq().curr_d_minus_s()->MakeNewCopy();
    while (count_soc<max_soc_ && !accept &&
           (count_soc==0 || (theta_trial<=kappa_soc_*theta_soc_old ||
                             theta_trial2<=kappa_soc_*theta_soc_old2)) ) {
      theta_soc_old = theta_trial;
      theta_soc_old2 = theta_trial2;
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Trying second order correction number %d\n",
                     count_soc+1);
      // Compute SOC constraint violation
      /*
      Number c_over_r = 0.;
      if (IpData().BiggerJacPert()){
      c_over_r = IpCq().curr_cg_pert_fact();
      }*/
      c_soc->AddTwoVectors(1.0, *IpCq().trial_c(),
                           -CGPenData().CurrPenaltyPert(), *delta_y_c,
                           alpha_primal_soc);
      dms_soc->AddTwoVectors(1.0, *IpCq().trial_d_minus_s(),
                             -CGPenData().CurrPenaltyPert(), *delta_y_d,
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
      catch (IpoptNLP::Eval_Error& e) {
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
        theta_trial2 = IpCq().trial_primal_infeasibility(NORM_2);
      }
    }
    if (accept) {
      ls_counter_ = 0;
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
    char info_alpha_primal_char='n';
    if (pen_curr_mu_ > IpData().curr_mu()) {
      pen_curr_mu_ = IpData().curr_mu();
      best_KKT_error_ = -1.;
    }
    // See if the current iterate has the least KKT errors so far
    // If so, store the current iterate
    if (CurrentIsBest()) {
      StoreBestPoint();
    }
    // update piecewise penalty parameters
    PiecewisePenalty_.Print( Jnlst() );
    if (!accepted_by_Armijo_) {
      PiecewisePenalty_.UpdateEntry(IpCq().trial_barrier_obj(),
                                    IpCq().trial_constraint_violation());
    }
    PiecewisePenalty_.Print( Jnlst() );

    // update regular penalty parameter
    if (CGPenData().CurrPenaltyPert() != 0) {
      info_alpha_primal_char = UpdatePenaltyParameter();
    }
    return info_alpha_primal_char;
  }

  bool CGPenaltyLSAcceptor::RestoredIterate()
  {
    bool restored_iterate = false;
    // See if we want to restor a previous iterate.
    // This is the case when the multipliers are blowing up.
    //if (!CGPenData().NeverTryPureNewton()){
    if (CGPenData().restor_counter()<3.) {
      if (MultipliersDiverged()) {
        if (RestoreBestPoint()) {
          // so far the restoration method is wrong. Although it restores
          // both the primal and dual iterates here, the dual iterates
          // may be overwritten later when performing dual steps
          // in the line search method
          Index restor_iter = IpData().iter_count() + 1;
          Number restor_counter = CGPenData().restor_counter();
          CGPenData().SetRestorCounter(restor_counter+1);
          CGPenData().SetNeverTryPureNewton(true);
          CGPenData().SetRestorIter(restor_iter);
          restored_iterate = true;
        }
      }
    }
    //}

    return restored_iterate;
  }

  bool CGPenaltyLSAcceptor::NeverRestorationPhase()
  {
    return true;
  }
  bool CGPenaltyLSAcceptor::CurrentIsBest()
  {
    DBG_START_METH("CGPenaltyLSAcceptor::CurrentIsBest",
                   dbg_verbosity);

    Number dual_inf = IpCq().unscaled_curr_dual_infeasibility(NORM_MAX);
    Number constr_viol = IpCq().unscaled_curr_nlp_constraint_violation(NORM_MAX);
    Number compl_inf = IpCq().unscaled_curr_complementarity(0., NORM_MAX);
    Number KKT_error = Max(dual_inf,Max(constr_viol,compl_inf));
    DBG_PRINT((1, "dual_inf = %e\n", dual_inf));
    DBG_PRINT((1, "constr_viol = %e\n", constr_viol));
    DBG_PRINT((1, "compl_inf = %e\n", compl_inf));
    DBG_PRINT((1, "best_KKT_error_= %e\n", best_KKT_error_));
    bool best = false;
    if (KKT_error < best_KKT_error_ ||
        best_KKT_error_ < 0.) {
      best_KKT_error_ = KKT_error;
      best = true;
    }
    return best;
  }

  void CGPenaltyLSAcceptor::StoreBestPoint()
  {
    DBG_START_METH("CGPenaltyLSAcceptor::StoreBestPoint",
                   dbg_verbosity);

    best_iterate_ = IpData().curr();
  }

  bool CGPenaltyLSAcceptor::RestoreBestPoint()
  {
    DBG_START_METH("CGPenaltyLSAcceptor::RestoreBestPoint",
                   dbg_verbosity);

    if (!IsValid(best_iterate_)) {
      return false;
    }
    SmartPtr<IteratesVector> prev_iterate = best_iterate_->MakeNewContainer();
    IpData().set_trial(prev_iterate);
    return true;
  }

  bool CGPenaltyLSAcceptor::MultipliersDiverged()
  {
    DBG_START_METH("CGPenaltyLSAcceptor::MultipliersDiverged",
                   dbg_verbosity);

    bool diverged = false;
    Number curr_inf = IpCq().curr_primal_infeasibility(NORM_2);
    Number trial_inf = IpCq().trial_primal_infeasibility(NORM_2);
    if (curr_inf > mult_diverg_feasibility_tol_
        && trial_inf > mult_diverg_feasibility_tol_
        && IpCq().curr_dual_infeasibility(NORM_MAX) > 1e4) {
      Number y_ref_big_step = mult_diverg_y_tol_;
      Number y_ref_tiny_step = 1e4;
      Number alpha_ref = 1e-4;
      /*
      if (IpCq().trial_grad_f()->Amax() < IpData().tol()){
        y_ref_big_step *= 1e4;
      y_ref_tiny_step *= 1e4;
      }
      */
      Number y_Amax = CGPenCq().curr_scaled_y_Amax();
      if ((y_Amax > y_ref_big_step &&
           (IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim() +
            IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim() +
            IpData().curr()->y_d()->Dim() == 0 ||
            CGPenData().PrimalStepSize() < 1e-2)) ||
          (CGPenData().PrimalStepSize() < alpha_ref &&
           y_Amax > y_ref_tiny_step)) {
        diverged = true;
      }
    }
    return diverged;
  }

  char CGPenaltyLSAcceptor::UpdatePenaltyParameter()
  {
    DBG_START_METH("CGPenaltyLSAcceptor::UpdatePenaltyParameter",
                   dbg_verbosity);
    char info_alpha_primal_char = 'n';
    // We use the new infeasibility here...
    Number trial_inf = IpCq().trial_primal_infeasibility(NORM_2);
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "trial infeasibility = %8.2\n", trial_inf);
    if (curr_eta_<0.) {
      // We need to initialize the eta tolerance
      curr_eta_ = Max(eta_min_, Min(gamma_tilde_,
                                    gamma_hat_*IpCq().curr_nlp_error()));
    }
    // Check if the penalty parameter is to be increased
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "Starting tests for penalty parameter update:\n");
    bool increase = (trial_inf >= penalty_update_infeasibility_tol_);
    if (!increase) {
      info_alpha_primal_char='i';
    }
    if (increase) {
      Number max_step = Max(CGPenData().delta_cgpen()->x()->Amax(),
                            CGPenData().delta_cgpen()->s()->Amax());
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
      vec->AddTwoVectors(1., *CGPenData().delta_cgpen()->y_c(),
                         -1./CGPenCq().curr_cg_pert_fact(), *IpCq().curr_c(),
                         1.);
      Number omega_test = vec->Amax();
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "omega_test for c = %8.2\n", omega_test);
      increase = (omega_test < curr_eta_);
      if (increase) {
        SmartPtr<Vector> vec = IpData().curr()->y_d()->MakeNewCopy();
        vec->AddTwoVectors(1., *IpData().delta()->y_d(),
                           -1./CGPenCq().curr_cg_pert_fact(), *IpCq().curr_d_minus_s(),
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
      counter_first_type_penalty_updates_++;
      // Update the eta tolerance
      curr_eta_ = Max(eta_min_, curr_eta_/2.);
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Updating eta to = %8.2\n", curr_eta_);
      Number penalty = CGPenData().curr_kkt_penalty();
      Number y_full_step_max;
      SmartPtr<Vector> vec = IpData().curr()->y_c()->MakeNew();
      vec->AddTwoVectors(1., *IpData().curr()->y_c(),
                         1., *CGPenData().delta_cgpen()->y_c(), 0.);
      y_full_step_max = vec->Amax();
      vec = IpData().curr()->y_d()->MakeNew();
      vec->AddTwoVectors(1., *IpData().curr()->y_d(),
                         1., *CGPenData().delta_cgpen()->y_d(), 0.);
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
      CGPenData().Set_kkt_penalty(penalty);
      if (CGPenData().NeverTryPureNewton()) {
        CGPenData().Set_penalty(penalty);
      }
    }

    // Second heuristic update for penalty parameters
    if (IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim() > 0 &&
        !never_use_piecewise_penalty_ls_) {
      Number scaled_y_Amax = CGPenCq().curr_scaled_y_Amax();
      if (scaled_y_Amax <= 1e4 ||
          counter_second_type_penalty_updates_ < 5) {
        Number result;
        SmartPtr<const Vector> ty_c = IpData().curr()->y_c();
        SmartPtr<const Vector> ty_d = IpData().curr()->y_d();
        SmartPtr<const Vector> dy_c = IpData().delta()->y_c();
        SmartPtr<const Vector> dy_d = IpData().delta()->y_d();
        Number curr_inf = IpCq().curr_primal_infeasibility(NORM_2);
        result = dy_c->Dot(*IpCq().curr_c()) + dy_d->Dot(*IpCq().curr_d_minus_s());
        if (!CGPenData().HaveCgFastDeltas()) {
          result += ty_c->Dot(*IpCq().curr_c()) + ty_d->Dot(*IpCq().curr_d_minus_s());
        }
        Number k_pen = CGPenData().curr_kkt_penalty();
        if (result > 0.5*k_pen*curr_inf || result < -0.5*k_pen*curr_inf) {
          Number nrm2_y = CGPenCq().curr_added_y_nrm2();
          result = 5.*nrm2_y;
          CGPenData().Set_kkt_penalty(result);
          if (CGPenData().NeverTryPureNewton()) {
            CGPenData().Set_penalty(result);
          }
          if (scaled_y_Amax > 1e4) {
            counter_second_type_penalty_updates_++;
          }
        }
      }
    }
    return info_alpha_primal_char;
  }

  bool CGPenaltyLSAcceptor::DoFallback()
  {
    // LIFENG: REPLACE ME!
    if (RestoreBestPoint()) {
      Index restor_iter = IpData().iter_count() + 1;
      CGPenData().SetRestorIter(restor_iter);
      CGPenData().SetNeverTryPureNewton(true);
      IpData().Append_info_string("help");
      return true;
    }
    else {
      return false;
    }
  }

  void CGPenaltyLSAcceptor::PrepareRestoPhaseStart()
  {
    DBG_ASSERT(false && "PrepareRestoPhaseStart not yet implemented");
  }


} // namespace Ipopt
