// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2005-10-13
//               derived from IpIpoptAlg.cpp
//
//           Lifeng Chen/Zaiwen Wen      Columbia Univ

#include "IpCGPenaltyData.hpp"
#include "IpCGPenaltyCq.hpp"
#include "IpCGSearchDirCalc.hpp"

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

  CGSearchDirCalculator::CGSearchDirCalculator(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      pd_solver_(pd_solver)
  {
    DBG_START_FUN("CGSearchDirCalculator::CGSearchDirCalculator",
                  dbg_verbosity);
    DBG_ASSERT(IsValid(pd_solver_));
  }

  CGSearchDirCalculator::~CGSearchDirCalculator()
  {
    DBG_START_FUN("CGSearchDirCalculator::~CGSearchDirCalculator()",
                  dbg_verbosity);
  }

  void CGSearchDirCalculator::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "penalty_init_max",
      "Maximal value for the intial penalty parameter (for Chen-Goldfarb line search).",
      0, true, 1e5,
      "");
    roptions->AddLowerBoundedNumberOption(
      "penalty_init_min",
      "Minimal value for the intial penalty parameter for line search(for Chen-Goldfarb line search).",
      0., true, 1.,
      "");
    roptions->AddLowerBoundedNumberOption(
      "penalty_max",
      "Maximal value for the penalty parameter (for Chen-Goldfarb line search).",
      0, true, 1e30,
      "");
    roptions->AddLowerBoundedNumberOption(
      "pen_des_fact",
      "a parameter used in penalty parameter computation (for Chen-Goldfarb line search).",
      0.0, true, 2e-1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "kappa_x_dis",
      "a parameter used to check if the fast direction can be used as"
      "the line search direction (for Chen-Goldfarb line search).",
      0.0, true, 1e2,
      "");
    roptions->AddLowerBoundedNumberOption(
      "kappa_y_dis",
      "a parameter used to check if the fast direction can be used as"
      "the line search direction (for Chen-Goldfarb line search).",
      0.0, true, 1e4,
      "");
    roptions->AddLowerBoundedNumberOption(
      "vartheta",
      "a parameter used to check if the fast direction can be used as"
      "the line search direction (for Chen-Goldfarb line search).",
      0.0, true, 0.5,
      "");
    roptions->AddLowerBoundedNumberOption(
      "delta_y_max",
      "a parameter used to check if the fast direction can be used as"
      "the line search direction (for Chen-Goldfarb line search).",
      0.0, true, 1e12,
      "");
    roptions->AddLowerBoundedNumberOption(
      "fast_des_fact",
      "a parameter used to check if the fast direction can be used as"
      "the line search direction (for Chen-Goldfarb line search).",
      0.0, true, 1e-1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "pen_init_fac",
      "a parameter used to choose initial penalty parameters"
      "when the regularized Newton method is used.",
      0.0, true, 5e1,
      "");
    roptions->AddStringOption2(
      "never_use_fact_cgpen_direction",
      "Toggle to switch off the fast Chen-Goldfarb direction",
      "no",
      "no", "always compute the fast direction",
      "yes", "never compute the fast direction",
      "");
  }

  bool CGSearchDirCalculator::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("penalty_init_max", penalty_init_max_, prefix);
    options.GetNumericValue("penalty_init_min", penalty_init_min_, prefix);
    options.GetNumericValue("penalty_max", penalty_max_, prefix);
    options.GetNumericValue("kappa_x_dis", kappa_x_dis_, prefix);
    options.GetNumericValue("kappa_y_dis", kappa_y_dis_, prefix);
    options.GetNumericValue("vartheta", vartheta_, prefix);
    options.GetNumericValue("delta_y_max", delta_y_max_, prefix);
    options.GetNumericValue("fast_des_fact", fast_des_fact_, prefix);
    options.GetNumericValue("pen_des_fact", pen_des_fact_, prefix);
    options.GetNumericValue("pen_init_fac", pen_init_fac_, prefix);
    options.GetBoolValue("never_use_fact_cgpen_direction",
                         never_use_fact_cgpen_direction_, prefix);
    options.GetNumericValue("penalty_init_min", penalty_init_min_, prefix);
    nonmonotone_pen_update_counter_ = 0;

    return pd_solver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                  options, prefix);
  }

  bool CGSearchDirCalculator::ComputeSearchDirection()
  {
    DBG_START_METH("CGSearchDirCalculator::ComputeSearchDirection",
                   dbg_verbosity);

    bool improve_solution = false;

    // So far, the adaptive mu strategies do not yet work with the
    // penalty function line search
    DBG_ASSERT(!IpData().FreeMuMode());

    SmartPtr<IteratesVector> rhs = IpData().curr()->MakeNewContainer();
    rhs->Set_x(*IpCq().curr_grad_lag_with_damping_x());
    rhs->Set_s(*IpCq().curr_grad_lag_with_damping_s());
    rhs->Set_z_L(*IpCq().curr_relaxed_compl_x_L());
    rhs->Set_z_U(*IpCq().curr_relaxed_compl_x_U());
    rhs->Set_v_L(*IpCq().curr_relaxed_compl_s_L());
    rhs->Set_v_U(*IpCq().curr_relaxed_compl_s_U());

    /** Initialize the penalty parameter */
    if (!CGPenData().PenaltyInitialized() ||
        !CGPenData().KKTPenaltyInitialized()) {
      Number penalty_init = penalty_init_min_;
      Number kkt_penalty_init = penalty_init_min_;
      if (!CGPenData().NeverTryPureNewton()) {
        Number y_max = Max(IpData().curr()->y_c()->Amax(),
                           IpData().curr()->y_d()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                       "Initializing penalty parameter for KKT matrix...\n");
        Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                       "Max(||y_c||_inf,||y_d||_inf = %8.2e\n",
                       y_max);
        penalty_init =  Max(penalty_init_min_, Min(y_max, penalty_init_max_));
      }
      else {
        // penalty_init = CGPenCq().compute_curr_cg_penalty_scale();
        // penalty_init = Max(penalty_init_min_, Min(penalty_init, penalty_init_max_));
        // For the moment,let's just not do scale
        penalty_init = 1e2*IpCq().curr_primal_infeasibility(NORM_2);
        penalty_init = Min(1e5, Max(1e1,penalty_init));
        kkt_penalty_init = penalty_init;
      }
      CGPenData().Set_penalty(penalty_init);
      CGPenData().Set_kkt_penalty(kkt_penalty_init);
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Initial value of the penalty parameter for line search = %8.2e\n",
                     penalty_init);
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Initial value of the kkt penalty parameter for scaling the linear system = %8.2e\n",
                     kkt_penalty_init);
    }
    else {
      /*
      if (CGPenData().NeverTryPureNewton()) {
        // If pure Newton method is not used, adjust the penalty para for
        // purposes of well scaling the KKT system.
        if (CGPenData().restor_iter() == IpData().iter_count()) {
          Number restor_penalty_init = CGPenCq().compute_curr_cg_penalty_scale();
          restor_penalty_init = Max(penalty_init_min_,
                                    Min(restor_penalty_init, penalty_init_max_));
          CGPenData().Set_penalty(restor_penalty_init);
          CGPenData().Set_kkt_penalty(restor_penalty_init);
        }
        else {
          Number penalty_new = CGPenCq().compute_curr_cg_penalty_scale();
          penalty_new = Max(CGPenData().curr_penalty(),
                            Min(penalty_new, penalty_max_));
          CGPenData().Set_penalty(penalty_new);
          CGPenData().Set_kkt_penalty(penalty_new);

        }
      }*/

      if (CGPenData().restor_iter() == IpData().iter_count()) {
        Number i = CGPenData().restor_counter();
        Number fac = pen_init_fac_*pow(1e-1, i);
        //Number restor_penalty_init = fac*IpCq().curr_primal_infeasibility(NORM_2);
        Number restor_penalty_init = fac;
        restor_penalty_init = Min(1e6, Max(1e1,restor_penalty_init));
        CGPenData().Set_penalty(restor_penalty_init);
        CGPenData().Set_kkt_penalty(restor_penalty_init);
      }
    }


    // Initialize iteration data
    CGPenData().SetCurrPenaltyPert(0.);
    CGPenData().SetPrimalStepSize(1.);

    // Compute the fast direction
    rhs->Set_y_c(*IpCq().curr_c());
    rhs->Set_y_d(*IpCq().curr_d_minus_s());
    // Get space for the search direction
    SmartPtr<IteratesVector> delta_fast =
      IpData().curr()->MakeNewIteratesVector(true);
    bool allow_inexact = false;
    bool retval = pd_solver_->Solve(-1.0, 0.0, *rhs, *delta_fast,
                                    allow_inexact, improve_solution);
    if (!retval) {
      return false;
    }
    // Store the fast search direction in the IpData->CGPenData object
    CGPenData().set_delta_cgfast(delta_fast);
    CGPenData().SetHaveCgFastDeltas(true);

    bool keep_fast_delta = true;

    // Get space for the cg_pen search direction
    SmartPtr<IteratesVector> delta_cgpen =
      IpData().curr()->MakeNewIteratesVector(true);

    if (CGPenData().CurrPenaltyPert() == 0.) {
      // if there is no perturbation on the Jacob, delta_fast = delta_cgpen.
      delta_cgpen->AddOneVector(1., *CGPenData().delta_cgfast(), 0.);
      CGPenData().set_delta_cgpen(delta_cgpen);
      CGPenData().SetHaveCgPenDeltas(true);
    }
    else {
      SmartPtr<Vector> rhs_c = IpData().curr()->y_c()->MakeNew();
      Number curr_pen_pert = CGPenData().CurrPenaltyPert();
      rhs_c->AddTwoVectors(1., *IpCq().curr_c(),
                           -curr_pen_pert, *IpData().curr()->y_c(), 0.);
      rhs->Set_y_c(*rhs_c);
      SmartPtr<Vector> rhs_d = IpData().curr()->y_d()->MakeNew();
      rhs_d->AddTwoVectors(1., *IpCq().curr_d_minus_s(),
                           -curr_pen_pert, *IpData().curr()->y_d(), 0.);
      rhs->Set_y_d(*rhs_d);
      DBG_PRINT_VECTOR(2, "rhs_cgpen", *rhs);
      allow_inexact = false;
      retval = pd_solver_->Solve(-1.0, 0.0, *rhs, *delta_cgpen, allow_inexact,
                                 improve_solution);
      if (!retval) {
        return false;
      }
      // Store the original search direction in the IpData object
      CGPenData().set_delta_cgpen(delta_cgpen);
      CGPenData().SetHaveCgPenDeltas(true);
      // Now we check whether the fast direction is good compatible with
      // the merit function
      // do the || tilde d_x - bar d_x ||_2 <= k_x_dis ||tilde d_x||_2^{vartheta} test
      SmartPtr<const Vector> delta_fast_x = CGPenData().delta_cgfast()->x();
      SmartPtr<const Vector> delta_fast_s = CGPenData().delta_cgfast()->s();
      SmartPtr<const Vector> delta_x = CGPenData().delta_cgpen()->x();
      SmartPtr<const Vector> delta_s = CGPenData().delta_cgpen()->s();
      Number tilde_dx_nrm = sqrt(pow(delta_fast_x->Nrm2(), 2.)
                                 + pow(delta_fast_s->Nrm2(), 2.));
      Number diff_dx_nrm = sqrt(pow(delta_fast_x->Nrm2(), 2.)
                                + pow(delta_fast_s->Nrm2(), 2.)
                                - 2.*delta_x->Dot(*delta_fast_x)
                                - 2.*delta_s->Dot(*delta_fast_s)
                                + pow(delta_x->Nrm2(), 2.)
                                + pow(delta_s->Nrm2(), 2.));
      Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                     "Testing if fast direction can be used.\n"
                     "  diff_dx_nrm = %8.2e tilde_dx_norm = %8.2e\n",
                     diff_dx_nrm, tilde_dx_nrm);
      tilde_dx_nrm = Max(tilde_dx_nrm, pow(tilde_dx_nrm, vartheta_));
      if (diff_dx_nrm > kappa_x_dis_*tilde_dx_nrm) {
        keep_fast_delta = false;
      }
      if (keep_fast_delta) {
        // do the || tilde d_y  ||_2 <= Max(d_y_max, k_y_dis ||  y + bar d_y  ||_2)
        SmartPtr<const Vector> y_c = IpData().curr()->y_c();
        SmartPtr<const Vector> y_d = IpData().curr()->y_d();
        SmartPtr<const Vector> delta_fast_y_c = CGPenData().delta_cgfast()->y_c();
        SmartPtr<const Vector> delta_fast_y_d = CGPenData().delta_cgfast()->y_d();
        SmartPtr<const Vector> delta_y_c = CGPenData().delta_cgpen()->y_c();
        SmartPtr<const Vector> delta_y_d = CGPenData().delta_cgpen()->y_d();
        Number tilde_dy_nrm = sqrt(pow(delta_fast_y_c->Nrm2(), 2.)
                                   + pow(delta_fast_y_d->Nrm2(), 2.));
        Number bar_y_nrm = sqrt(pow(y_c->Nrm2(), 2.)
                                + pow(y_d->Nrm2(), 2.)
                                + 2.*y_c->Dot(*delta_y_c)
                                + 2.*y_d->Dot(*delta_y_d)
                                + pow(delta_y_c->Nrm2(), 2.)
                                + pow(delta_y_d->Nrm2(), 2.));
        Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                       "Testing if fast direction can be used.\n"
                       "  tilde_dy_nrm = %8.2e bar_y_nrm = %8.2e\n",
                       tilde_dy_nrm, bar_y_nrm);
        if (tilde_dy_nrm > Max(delta_y_max_, kappa_y_dis_*bar_y_nrm)) {
          keep_fast_delta = false;
        }
      }
      if (keep_fast_delta) {
        // For now, I just check if the directional derivative for the
        // penalty functions are not too much off
        Number dT_times_BarH_times_d = CGPenCq().dT_times_barH_times_d();
        Number fast_direct_deriv = CGPenCq().curr_fast_direct_deriv_penalty_function();
        Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                       "dT_times_BarH_times_d = %23.15e  fast_direct_deriv = %23.15e\n",
                       dT_times_BarH_times_d, fast_direct_deriv);
        if (fast_direct_deriv > fast_des_fact_*dT_times_BarH_times_d) {
          keep_fast_delta = false;
          IpData().Append_info_string("g");
        }
      }
    }

    // Set search direction
    SmartPtr<IteratesVector> delta =
      IpData().curr()->MakeNewIteratesVector(true);
    if (!keep_fast_delta) {
      CGPenData().SetHaveCgFastDeltas(false);
      delta->AddOneVector(1., *CGPenData().delta_cgpen(), 0.);
    }
    else {
      CGPenData().SetHaveCgFastDeltas(true);
      delta->AddOneVector(1., *CGPenData().delta_cgfast(), 0.);
    }
    IpData().set_delta(delta);
    // If the pure Newton method is used, update the penalty parameter for line search
    if (!CGPenData().NeverTryPureNewton()) {
      Number penalty = CGPenCq().compute_curr_cg_penalty(pen_des_fact_);
      Number curr_penalty = CGPenData().curr_penalty();
      Number curr_kkt_penalty = CGPenData().curr_kkt_penalty();
      if (penalty>curr_penalty) {
        penalty = Max(penalty, curr_penalty+1.);
      }
      else {
        if (curr_penalty <= curr_kkt_penalty ||
            CGPenData().CurrPenaltyPert() == 0.) {
          penalty = curr_penalty;
        }
        else {
          penalty = curr_kkt_penalty;
          nonmonotone_pen_update_counter_++;
        }
      }
      CGPenData().Set_penalty(penalty);
      if (penalty > curr_kkt_penalty &&
          nonmonotone_pen_update_counter_ > 50) {
        CGPenData().Set_kkt_penalty(penalty);
      }
    }

    return true;
  }

} // namespace Ipopt
