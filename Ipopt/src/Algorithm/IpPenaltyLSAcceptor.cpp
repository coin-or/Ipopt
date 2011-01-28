// Copyright (C) 2008, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2008-04-04
//               derived file from IpFilterLSAcceptor.cpp

#include "IpPenaltyLSAcceptor.hpp"

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

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  PenaltyLSAcceptor::PenaltyLSAcceptor(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      pd_solver_(pd_solver)
  {
    DBG_START_FUN("PenaltyLSAcceptor::PenaltyLSAcceptor",
                  dbg_verbosity);
  }

  PenaltyLSAcceptor::~PenaltyLSAcceptor()
  {
    DBG_START_FUN("PenaltyLSAcceptor::~PenaltyLSAcceptor()",
                  dbg_verbosity);
  }

  void PenaltyLSAcceptor::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "nu_init",
      "Initial value of the penalty parameter.",
      0.0, true, 1e-6,
      "");
    roptions->AddLowerBoundedNumberOption(
      "nu_inc",
      "Increment of the penalty parameter.",
      0.0, true, 1e-4,
      "");
    roptions->AddBoundedNumberOption(
      "rho",
      "Value in penalty parameter update formula.",
      0.0, true, 1.0, true, 1e-1,
      "");
    
  }

  bool PenaltyLSAcceptor::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    options.GetNumericValue("nu_init", nu_init_, prefix);
    options.GetNumericValue("nu_inc", nu_inc_, prefix);
    options.GetNumericValue("eta_phi", eta_, prefix);
    options.GetNumericValue("rho", rho_, prefix);

    // The following options have been declared in FilterLSAcceptor
    options.GetIntegerValue("max_soc", max_soc_, prefix);
    if (max_soc_>0) {
      ASSERT_EXCEPTION(IsValid(pd_solver_), OPTION_INVALID,
                       "Option \"max_soc\": This option is non-negative, but no linear solver for computing the SOC given to PenaltyLSAcceptor object.");
    }
    options.GetNumericValue("kappa_soc", kappa_soc_, prefix);

    Reset();

    return true;
  }

  void PenaltyLSAcceptor::InitThisLineSearch(bool in_watchdog)
  {
    DBG_START_METH("PenaltyLSAcceptor::InitThisLineSearch",
                   dbg_verbosity);

    // Set the values for the reference point
    if (!in_watchdog) {
      reference_theta_ = IpCq().curr_constraint_violation();
      reference_barr_ = IpCq().curr_barrier_obj();
      reference_gradBarrTDelta_ = IpCq().curr_gradBarrTDelta();

      // Compute d_x^T (W + \Sigma_x + delta_x I) d_x +
      //         d_s^T (\Sigma_s + delta_s I) d_s
      Number pd_pert_x;
      Number pd_pert_s;
      Number pd_pert_c;
      Number pd_pert_d;
      IpData().getPDPert(pd_pert_x, pd_pert_s, pd_pert_c, pd_pert_d);

      SmartPtr<const Vector> delta_x = IpData().delta()->x();
      SmartPtr<Vector> tmp = delta_x->MakeNew();
      IpData().W()->MultVector(1., *delta_x, 0., *tmp);
      reference_dWd_ = tmp->Dot(*delta_x);
      tmp->Copy(*delta_x);
      tmp->ElementWiseMultiply(*IpCq().curr_sigma_x());
      reference_dWd_ += tmp->Dot(*delta_x);
      if (pd_pert_x != 0.) {
	Number dnrm2 = delta_x->Nrm2();
	reference_dWd_ += pd_pert_x * dnrm2*dnrm2;
      }
      SmartPtr<const Vector> delta_s = IpData().delta()->s();
      tmp = delta_s->MakeNewCopy();
      tmp->ElementWiseMultiply(*IpCq().curr_sigma_s());
      reference_dWd_ += tmp->Dot(*delta_s);
      if (pd_pert_s != 0.) {
	Number dnrm2 = delta_s->Nrm2();
	reference_dWd_ += pd_pert_s * dnrm2*dnrm2;
      }
      // Set back to zero, if negative
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "  dWd = %23.16e\n",
		     reference_dWd_);
      if (reference_dWd_ <= 0.) {
	reference_dWd_ = 0.;
      }
      //reference_dWd_ = 0.;

      // Get the product of the steps with the Jacobian
      reference_JacC_delta_ = IpCq().curr_jac_c_times_vec(*delta_x);
      tmp = delta_s->MakeNew();
      tmp->AddTwoVectors(1., *IpCq().curr_jac_d_times_vec(*delta_x),
			 -1., *delta_s,  0.);
      reference_JacD_delta_ = ConstPtr(tmp);

      reference_pred_ = -1.;
      resto_pred_ = -1;

      // update the penalty parameter
      last_nu_ = nu_;
      if (reference_theta_ > 0.) {
	Number nu_plus = (reference_gradBarrTDelta_+reference_dWd_/2.)/((1.-rho_)*reference_theta_);
	if (nu_ < nu_plus) {
	  nu_ = nu_plus + nu_inc_;
	}
      }
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "  using nu = %23.16e\n", nu_);
    }
    else {
      reference_theta_ = watchdog_theta_;
      reference_barr_ = watchdog_barr_;
      reference_pred_ = watchdog_pred_;
    }
  }

  Number
  PenaltyLSAcceptor::CalcPred(Number alpha)
  {
    DBG_START_METH("PenaltyLSAcceptor::CalcPred",
                   dbg_verbosity);

    SmartPtr<const Vector> curr_c = IpCq().curr_c();
    SmartPtr<const Vector> curr_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<Vector> tmp_c = curr_c->MakeNew();
    SmartPtr<Vector> tmp_d = curr_d_minus_s->MakeNew();
    tmp_c->AddTwoVectors(1., *curr_c, alpha, *reference_JacC_delta_, 0.);
    tmp_d->AddTwoVectors(1., *curr_d_minus_s, alpha, *reference_JacD_delta_, 0.);

    Number theta2 = IpCq().CalcNormOfType(IpCq().constr_viol_normtype(),
					  *tmp_c, *tmp_d);
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH, "  theta2 = %23.16e\n", theta2);

    Number pred = -alpha*reference_gradBarrTDelta_ - alpha*alpha/2.*reference_dWd_ +
      nu_*(reference_theta_ - theta2);

    if (pred < 0.) {
      Jnlst().Printf(J_WARNING, J_LINE_SEARCH, "  pred = %23.16e is negative.  Setting to zero.\n", pred);
      pred = 0.;
    }

    return pred;
  }

  bool
  PenaltyLSAcceptor::CheckAcceptabilityOfTrialPoint(Number alpha_primal_test)
  {
    DBG_START_METH("PenaltyLSAcceptor::CheckAcceptabilityOfTrialPoint",
                   dbg_verbosity);

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

    Number pred;
    if (reference_pred_ < 0.) {
      pred = CalcPred(alpha_primal_test);
    }
    else {
      pred = reference_pred_;
    }
    resto_pred_ = pred;
    Number ared = reference_barr_ + nu_*(reference_theta_) -
      (trial_barr + nu_*trial_theta);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
		   "  Checking Armijo Condition with pred = %23.16e and ared = %23.16e\n", pred, ared);

    bool accept;
    if (Compare_le(eta_*pred, ared, reference_barr_ + nu_*(reference_theta_))) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "   Success...\n");
      accept = true;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "   Failed...\n");
      accept = false;
    }
    return accept;
  }

  Number PenaltyLSAcceptor::CalculateAlphaMin()
  {
    // ToDo: make better
    return 1e-16;
  }

  void PenaltyLSAcceptor::StartWatchDog()
  {
    DBG_START_FUN("PenaltyLSAcceptor::StartWatchDog", dbg_verbosity);
    THROW_EXCEPTION(OPTION_INVALID, "Watchdog not implemented for penalty function line search.  Set watchdog_shortened_iter_trigger to 0.");
  }

  void PenaltyLSAcceptor::StopWatchDog()
  {
    DBG_START_FUN("PenaltyLSAcceptor::StopWatchDog", dbg_verbosity);
    THROW_EXCEPTION(OPTION_INVALID, "Watchdog not implemented for penalty function line search.  Set watchdog_shortened_iter_trigger to 0.");
  }

  void PenaltyLSAcceptor::Reset()
  {
    DBG_START_FUN("PenaltyLSAcceptor::Reset", dbg_verbosity);

    nu_ = nu_init_;
  }

  bool
  PenaltyLSAcceptor::TrySecondOrderCorrection(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_METH("PenaltyLSAcceptor::TrySecondOrderCorrection",
                   dbg_verbosity);

    if (max_soc_==0) {
      return false;
    }

    bool accept = false;
    Index count_soc = 0;

    Number theta_soc_old = 0.;
    Number theta_trial = IpCq().trial_constraint_violation();
    Number alpha_primal_soc = alpha_primal;

    SmartPtr<Vector> c_soc = IpCq().curr_c()->MakeNew();
    SmartPtr<Vector> dms_soc = IpCq().curr_d_minus_s()->MakeNew();
    c_soc->Copy(*IpCq().curr_c());
    dms_soc->Copy(*IpCq().curr_d_minus_s());
    while (count_soc<max_soc_ && !accept &&
           (count_soc==0 || theta_trial<=kappa_soc_*theta_soc_old) ) {
      theta_soc_old = theta_trial;

      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Trying second order correction number %d\n",
                     count_soc+1);

      // Compute SOC constraint violation
      c_soc->AddOneVector(1.0, *IpCq().trial_c(), alpha_primal_soc);
      dms_soc->AddOneVector(1.0, *IpCq().trial_d_minus_s(), alpha_primal_soc);

      // Compute the SOC search direction
      SmartPtr<IteratesVector> delta_soc = actual_delta->MakeNewIteratesVector(true);
      SmartPtr<IteratesVector> rhs = actual_delta->MakeNewContainer();
      rhs->Set_x(*IpCq().curr_grad_lag_with_damping_x());
      rhs->Set_s(*IpCq().curr_grad_lag_with_damping_s());
      rhs->Set_y_c(*c_soc);
      rhs->Set_y_d(*dms_soc);
      rhs->Set_z_L(*IpCq().curr_relaxed_compl_x_L());
      rhs->Set_z_U(*IpCq().curr_relaxed_compl_x_U());
      rhs->Set_v_L(*IpCq().curr_relaxed_compl_s_L());
      rhs->Set_v_U(*IpCq().curr_relaxed_compl_s_U());

      bool retval = pd_solver_->Solve(-1.0, 0.0, *rhs, *delta_soc, true);
      if (!retval) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "The linear system could not be solved for the corrector step.\n");
        return false;
      }

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
      }
    }
    return accept;
  }

  bool
  PenaltyLSAcceptor::TryCorrector(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    return false;
  }

  char PenaltyLSAcceptor::UpdateForNextIteration(Number alpha_primal_test)
  {
    // delete some stuff
    reference_JacC_delta_ = NULL;
    reference_JacD_delta_ = NULL;

    char info_alpha_primal_char = ' ';
    // Augment the filter if required
    if (last_nu_ != nu_) {
      info_alpha_primal_char = 'n';
      char snu[40];
      sprintf(snu, " nu=%8.2e", nu_);
      IpData().Append_info_string(snu);
    }
    else {
      info_alpha_primal_char = 'k';
    }
    return info_alpha_primal_char;
  }

  void PenaltyLSAcceptor::PrepareRestoPhaseStart()
  {
  }

  bool
  PenaltyLSAcceptor::IsAcceptableToCurrentIterate(Number trial_barr,
						  Number trial_theta,
						  bool called_from_restoration /*=false*/) const
  {
    DBG_START_METH("PenaltyLSAcceptor::IsAcceptableToCurrentIterate",
                   dbg_verbosity);
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
