// Copyright (C) 2004, 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
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

#ifdef IP_DEBUG
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
  }

  bool CGPenaltyLSAcceptor::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("eta_penalty", eta_penalty_, prefix);

    // The following could become regular options
    penalty_update_infeasibility_tol_ = 1e-9;
    eta_min_ = 1e-2;
    penalty_update_compl_tol_ = 0.2;
    chi_hat_ = 2.;
    chi_tilde_ = 5.;
    chi_cup_ = 1.5;
    gamma_hat_ = 0.04;
    gamma_tilde_ = 4.;
    penalty_max_ = 1e20;
    epsilon_c_ = 1e-2;

    Reset();

    counter_penalty_updates_ = 0;
    curr_eta_ = -1.;
    IpData().SetPenaltyUninitialized();

    return true;
  }

  void CGPenaltyLSAcceptor::InitThisLineSearch(bool in_watchdog)
  {
    DBG_START_METH("CGPenaltyLSAcceptor::InitThisLineSearch",
                   dbg_verbosity);

    // Set the values for the reference point
    if (!in_watchdog) {
      reference_penalty_function_ = IpCq().curr_penalty_function();
      // ToDo we need to make this cleaner!
      if (IpData().HaveCgPenDeltas()) {
	// use the fast step
	reference_direct_deriv_penalty_function_ =
	  IpCq().curr_fast_direct_deriv_penalty_function();
	// Overwrite the original step (hmm.....)
	SmartPtr<const IteratesVector> delta_cgpen = IpData().delta_cgpen();
	IpData().set_delta(delta_cgpen);
      }
      else {
	reference_direct_deriv_penalty_function_ =
	  IpCq().curr_direct_deriv_penalty_function();
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

    Number trial_penalty_function = IpCq().trial_penalty_function();
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

    watchdog_penalty_function_ = IpCq().curr_penalty_function();
    watchdog_direct_deriv_penalty_function_ =
      IpCq().curr_direct_deriv_penalty_function();
  }

  void CGPenaltyLSAcceptor::StopWatchDog()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::StopWatchDog", dbg_verbosity);

    reference_penalty_function_ = watchdog_penalty_function_;
    reference_direct_deriv_penalty_function_ =
      watchdog_direct_deriv_penalty_function_;
  }

  void CGPenaltyLSAcceptor::Reset()
  {
    DBG_START_FUN("CGPenaltyLSAcceptor::Reset", dbg_verbosity);
    /*
    counter_penalty_updates_ = 0;
    curr_eta_ = -1.;
    IpData().SetPenaltyUninitialized();
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

    return false;
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
    char info_alpha_primal_char;

    if (curr_eta_<0.) {
      // We need to initialize the eta tolerance
      // Lifeng, at the moment I'm using the scaled 1-norm for E_0
      curr_eta_ = Max(eta_min_, Min(gamma_tilde_,
                                    gamma_hat_*IpCq().curr_primal_dual_system_error(0.)));
    }

    // Check if the penalty parameter is to be increased

    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
		   "Starting tests for penalty parameter update:\n");
    // Lifeng: Should we use the new infeasibility here?
    Number curr_inf = IpCq().curr_primal_infeasibility(NORM_2);
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
		   "current infeasibility = %8.2\n", curr_inf);
    bool increase = (curr_inf >= penalty_update_infeasibility_tol_);
    if (!increase) {
      info_alpha_primal_char='i';
    }

    if (increase) {
      Number max_step = Max(IpData().delta()->x()->Amax(),
			    IpData().delta()->s()->Amax());
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
      vec->AddTwoVectors(1., *IpData().delta()->y_c(),
                         -1./IpCq().curr_cg_pert_fact(), *IpCq().curr_c(),
                         1.);
      Number omega_test = vec->Amax();
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
		      "omega_test for c = %8.2\n", omega_test);
      increase = (omega_test < curr_eta_);
      if (increase) {
        SmartPtr<Vector> vec = IpData().curr()->y_d()->MakeNewCopy();
        vec->AddTwoVectors(1., *IpData().delta()->y_d(),
                           -1./IpCq().curr_cg_pert_fact(), *IpCq().curr_d_minus_s(),
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
      Number penalty = IpData().curr_penalty();
      Number y_full_step_max;
      SmartPtr<Vector> vec = IpData().curr()->y_c()->MakeNew();
      vec->AddTwoVectors(1., *IpData().curr()->y_c(),
                         1., *IpData().delta()->y_c(), 0.);
      y_full_step_max = vec->Amax();
      vec = IpData().curr()->y_d()->MakeNew();
      vec->AddTwoVectors(1., *IpData().curr()->y_d(),
                         1., *IpData().delta()->y_d(), 0.);
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
      IpData().Set_penalty(penalty);

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
