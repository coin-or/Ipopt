// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpFilterLineSearch.hpp"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
# include <ctype.h>
#else
# include <cmath>
# include <cctype>
#endif

namespace Ipopt
{

  static const Index dbg_verbosity = 0;

  FilterLineSearch::FilterLineSearch(const SmartPtr<RestorationPhase>& resto_phase,
                                     const SmartPtr<PDSystemSolver>& pd_solver)
      :
      LineSearch(),
      resto_phase_(resto_phase),
      pd_solver_(pd_solver),
      theta_min_(-1.0),
      theta_max_(-1.0),
      filter_(2)
  {}

  FilterLineSearch::~FilterLineSearch()
  {}

  bool FilterLineSearch::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    Number value = 0.0;
    Int ivalue = 0;

    // Check for the algorithm options
    if (options.GetNumericValue("theta_max_fact", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"theta_max_fact\": This value must be larger than 0.");
      theta_max_fact_ = value;
    }
    else {
      theta_max_fact_ = 1e4;
    }

    if (options.GetNumericValue("theta_min_fact", value, prefix)) {
      ASSERT_EXCEPTION(value > 0 && value < theta_max_fact_, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"theta_min_fact\": This value must be larger than 0 and less than theta_max_fact.");
      theta_min_fact_ = value;
    }
    else {
      theta_min_fact_ = 1e-4;
    }

    if (options.GetNumericValue("eta_phi", value, prefix)) {
      ASSERT_EXCEPTION(value > 0 && value < 0.5, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"eta_phi\": This value must be between 0 and 0.5.");
      eta_phi_ = value;
    }
    else {
      eta_phi_ = 1e-4;
    }

    if (options.GetNumericValue("delta", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"delta\": This value must be larger than 0.");
      delta_ = value;
    }
    else {
      delta_ = 1.0;
    }

    if (options.GetNumericValue("s_phi", value, prefix)) {
      ASSERT_EXCEPTION(value > 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"s_phi\": This value must be larger than 1.");
      s_phi_ = value;
    }
    else {
      s_phi_ = 2.3;
    }

    if (options.GetNumericValue("s_theta", value, prefix)) {
      ASSERT_EXCEPTION(value > 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"s_theta\": This value must be larger than 1.0.");
      s_theta_ = value;
    }
    else {
      s_theta_ = 1.1;
    }


    if (options.GetNumericValue("gamma_phi", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"gamma_phi\": This value must be between 0 and 1.");
      gamma_phi_ = value;
    }
    else {
      gamma_phi_ = 1e-5;
    }

    if (options.GetNumericValue("gamma_theta", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"gamma_theta\": This value must be between 0 and 1.");
      gamma_theta_ = value;
    }
    else {
      gamma_theta_ = 1e-5;
    }

    if (options.GetNumericValue("alpha_min_frac", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"alpha_min_frac\": This value must be > 0 and <= 1.");
      alpha_min_frac_ = value;
    }
    else {
      alpha_min_frac_ = 0.05;
    }

    if (options.GetNumericValue("alpha_red_factor", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"alpha_red_factor\": This value must be larger than 0.");
      alpha_red_factor_ = value;
    }
    else {
      alpha_red_factor_ = 0.5;
    }

    if (options.GetIntegerValue("max_soc", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue >= 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"max_soc\": This value must be non-negative.");
      max_soc_ = ivalue;
    }
    else {
      max_soc_ = 4;
    }
    if (max_soc_>0) {
      ASSERT_EXCEPTION(IsValid(pd_solver_), OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"max_soc\": This option is non-negative, but no linear solver for computing the SOC given to FilterLineSearch object.");
    }

    if (options.GetNumericValue("kappa_soc", value, prefix)) {
      ASSERT_EXCEPTION(value > 0., OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"kappa_soc\": This value must be larger than 0.");
      kappa_soc_ = value;
    }
    else {
      kappa_soc_ = 0.99;
    }

    if (options.GetNumericValue("obj_max_inc", value, prefix)) {
      ASSERT_EXCEPTION(value > 1., OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"obj_max_inc_\": This value must be larger than 1.");
      obj_max_inc_ = value;
    }
    else {
      obj_max_inc_ = 5.;
    }

    if (options.GetIntegerValue("magic_steps", ivalue, prefix)) {
      magic_steps_ = (ivalue != 0);
    }
    else {
      magic_steps_ = false;
    }

    if (options.GetIntegerValue("corrector_type", ivalue, prefix)) {
      corrector_type_ = ivalue;
    }
    else {
      corrector_type_ = 0;
    }

    if (options.GetIntegerValue("skip_corr_if_neg_curv", ivalue, prefix)) {
      skip_corr_if_neg_curv_ = (ivalue != 0);
    }
    else {
      skip_corr_if_neg_curv_ = true;
    }

    if (options.GetIntegerValue("ls_always_accept", ivalue, prefix)) {
      ls_always_accept_ = (ivalue != 0);
    }
    else {
      ls_always_accept_ = false;
    }

    if (options.GetNumericValue("corrector_compl_avrg_red_fact", value, prefix)) {
      ASSERT_EXCEPTION(value > 0., OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"corrector_compl_avrg_red_fact_\": This value must be positive1.");
      corrector_compl_avrg_red_fact_ = value;
    }
    else {
      corrector_compl_avrg_red_fact_ = 1.;
    }

    bool retvalue = true;
    if (IsValid(resto_phase_)) {
      retvalue = resto_phase_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                          options, prefix);
    }

    // ToDo decide if also the PDSystemSolver should be initialized here...

    return retvalue;
  }

  void FilterLineSearch::FindAcceptableTrialPoint()
  {
    DBG_START_METH("FilterLineSearch::FindAcceptableTrialPoint",
                   dbg_verbosity);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "--> Starting filter line search in iteration %d <--\n",
                   IpData().iter_count());

    // Get the search directions (this will store the actual search
    // direction, possibly including higher order corrections)
    SmartPtr<const Vector> actual_delta_x = IpData().delta_x();
    SmartPtr<const Vector> actual_delta_s = IpData().delta_s();
    SmartPtr<const Vector> actual_delta_y_c = IpData().delta_y_c();
    SmartPtr<const Vector> actual_delta_y_d = IpData().delta_y_d();
    SmartPtr<const Vector> actual_delta_z_L = IpData().delta_z_L();
    SmartPtr<const Vector> actual_delta_z_U = IpData().delta_z_U();
    SmartPtr<const Vector> actual_delta_v_L = IpData().delta_v_L();
    SmartPtr<const Vector> actual_delta_v_U = IpData().delta_v_U();

    bool goto_resto = false;
    if (actual_delta_x->Asum() + actual_delta_s->Asum() +
        actual_delta_y_c->Asum() + actual_delta_y_d->Asum() +
        actual_delta_z_L->Asum() + actual_delta_z_U->Asum() +
        actual_delta_v_L->Asum() + actual_delta_v_U->Asum() == 0. ) {
      // In this case, a search direction could not be computed, and
      // we should immediately go to the restoration phase.  ToDo: Cue
      // off of a something else than the norm of the search direction
      goto_resto = true;
    }

    bool accept = false;
    bool corr_taken = false;
    bool soc_taken = false;
    Index n_steps = 0;
    Number alpha_primal = 0.;

    if (!goto_resto) {
      // Compute smallest step size allowed
      Number alpha_min = CalculateAlphaMin();
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "minimal step size ALPHA_MIN = %E\n", alpha_min);

      // Start line search from primal fraction-to-the-boundary value
      Number alpha_primal_max =
        IpCq().curr_primal_frac_to_the_bound(IpData().curr_tau());
      alpha_primal = alpha_primal_max;

      // Step size used in ftype and armijo tests
      Number alpha_primal_test = alpha_primal;

      filter_.Print(Jnlst());

      if (corrector_type_!=0 &&
          (!skip_corr_if_neg_curv_ || IpData().info_regu_x()==0.) ) {
        // Before we do the actual backtracking line search for the
        // regular primal-dual search direction, let's see if a step
        // including a higher-order correctior is already acceptable
        accept = TryCorrector(alpha_primal_test,
                              alpha_primal,
                              actual_delta_x,
                              actual_delta_s,
                              actual_delta_y_c,
                              actual_delta_y_d,
                              actual_delta_z_L,
                              actual_delta_z_U,
                              actual_delta_v_L,
                              actual_delta_v_U);
      }
      if (accept) {
        corr_taken = true;
      }

      if (!accept) {
        // Loop over decreaseing step sizes until acceptable point is
        // found or until step size becomes too small

        while (alpha_primal>alpha_min ||
               n_steps == 0) { // always allow the "full" step if it is
          // acceptable (even if alpha_primal<=alpha_min)
          Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                         "Starting checks for alpha (primal) = %8.2e\n",
                         alpha_primal);

          try {
            // Compute the primal trial point
            IpData().SetTrialPrimalVariablesFromStep(alpha_primal,
                *actual_delta_x,
                *actual_delta_s);

            if (magic_steps_) {
              PerformMagicStep();
            }

            // If it is acceptable, stop the search
            accept = CheckAcceptabilityOfTrialPoint(alpha_primal_test);
          }
          catch(IpoptNLP::Eval_Error& e) {
            e.ReportException(Jnlst());
            Jnlst().Printf(J_WARNING, J_MAIN,
                           "Warning: Cutting back alpha due to evaluation error\n");
            accept = false;
          }

          if (accept) {
            break;
          }

          Number theta_curr = IpCq().curr_constraint_violation();
          Number theta_trial = IpCq().trial_constraint_violation();
          if (alpha_primal==alpha_primal_max &&       // i.e. first trial point
              theta_curr<=theta_trial && max_soc_>0) {
            // Try second order correction
            accept = TrySecondOrderCorrection(alpha_primal_test,
                                              alpha_primal,
                                              actual_delta_x,
                                              actual_delta_s,
                                              actual_delta_y_c,
                                              actual_delta_y_d,
                                              actual_delta_z_L,
                                              actual_delta_z_U,
                                              actual_delta_v_L,
                                              actual_delta_v_U);
          }
          if (accept) {
            soc_taken = true;
            break;
          }

          // Point is not yet acceptable, try a shorter one
          alpha_primal *= alpha_red_factor_;
          n_steps++;
        }
      }

      char info_alpha_primal_char;
      // Augment the filter if required
      if (!IsFtype(alpha_primal_test) || !ArmijoHolds(alpha_primal_test)) {
        AugmentFilter();
        info_alpha_primal_char = 'h';
      }
      else {
        info_alpha_primal_char = 'f';
      }
      if (soc_taken) {
        info_alpha_primal_char = toupper(info_alpha_primal_char);
      }
      IpData().Set_info_alpha_primal_char(info_alpha_primal_char);
      IpData().Set_info_ls_count(n_steps+1);
      if (corr_taken) {
        IpData().Append_info_string("C");
      }
    }

    // If line search has been aborted because the step size becomes too small,
    // go to the restoration phase
    if (!accept) {
      if (IsValid(resto_phase_)) {
        if (IpCq().curr_constraint_violation()==0.) {
          THROW_EXCEPTION(IpoptException, "Restoration phase called, but constraint violation is zero.");
        }

        // Augment the filter with the current point
        AugmentFilter();

        // Set the info fields for the first output line in the
        // restoration phase which reflects why the restoration phase
        // was called
        IpData().Set_info_alpha_primal(alpha_primal);
        IpData().Set_info_alpha_dual(0.);
        IpData().Set_info_alpha_primal_char('R');
        IpData().Set_info_ls_count(n_steps+1);

        accept = resto_phase_->PerformRestoration();
        if (!accept) {
          DBG_ASSERT(false && "Failed restoration phase!!!");
        }
      }
      else {
        //ToDo
        DBG_ASSERT(false && "No Restoration Phase given to this Filter Line Search Object!");
      }
    }
    else {
      // we didn't do the restoration phase and are now updating the
      // trial point
      IpData().SetTrialEqMultipilersFromStep(alpha_primal,
                                             *actual_delta_y_c,
                                             *actual_delta_y_d);
      Number alpha_dual_max =
        IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                      *actual_delta_z_L, *actual_delta_z_U,
                                      *actual_delta_v_L, *actual_delta_v_U);

      IpData().SetTrialBoundMutlipliersFromStep(alpha_dual_max,
          *actual_delta_z_L, *actual_delta_z_U,
          *actual_delta_v_L, *actual_delta_v_U);

      // Set some information for iteration summary output
      IpData().Set_info_alpha_primal(alpha_primal);
      IpData().Set_info_alpha_dual(alpha_dual_max);
    }
  }

  bool FilterLineSearch::IsFtype(Number alpha_primal_test)
  {
    DBG_START_METH("FilterLineSearch::IsFtype",
                   dbg_verbosity);
    Number curr_theta = IpCq().curr_constraint_violation();

    return (IpCq().curr_gradBarrTDelta() < 0.0 &&
            alpha_primal_test*pow(-IpCq().curr_gradBarrTDelta(),s_phi_) >
            delta_*pow(curr_theta,s_theta_));
  }

  void FilterLineSearch::AugmentFilter()
  {
    DBG_START_METH("FilterLineSearch::AugmentFilter",
                   dbg_verbosity);
    Number curr_barr = IpCq().curr_barrier_obj();
    Number curr_theta = IpCq().curr_constraint_violation();

    Number phi_add = curr_barr - gamma_phi_*curr_theta;
    Number theta_add = (1.-gamma_theta_)*curr_theta;

    filter_.AddEntry(phi_add, theta_add, IpData().iter_count());
  }

  bool
  FilterLineSearch::CheckAcceptabilityOfTrialPoint(Number alpha_primal_test)
  {
    DBG_START_METH("FilterLineSearch::CheckAcceptabilityOfTrialPoint",
                   dbg_verbosity);

    if (ls_always_accept_) {
      //       // We call the evaluation of a trial point once here, because
      //       // otherwise the SafeSlack mechanism in
      //       // IpoptCalculatedQuantities complains (it currently only ever
      //       // corrects trial point slacks)
      //       Number trial_barr = IpCq().trial_barrier_obj();
      return true;
    }

    bool accept;

    // First compute the barrier function and constraint violation at the
    // current iterate and the trial point
    Number curr_barr = IpCq().curr_barrier_obj();
    Number curr_theta = IpCq().curr_constraint_violation();

    Number trial_theta = IpCq().trial_constraint_violation();
    // Check if constraint violation is becoming too large
    if (theta_max_ < 0.0) {
      theta_max_ = theta_max_fact_*Max(1.0, curr_theta);
    }
    if (theta_min_ < 0.0) {
      theta_min_ = theta_min_fact_*Max(1.0, curr_theta);
    }

    if (theta_max_>0 && trial_theta>theta_max_) {
      return false;
    }

    Number trial_barr = IpCq().trial_barrier_obj();

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Checking acceptability for trial step size alpha_primal_test=%13.6e:\n", alpha_primal_test);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of barrier function     = %23.16e  (current %23.16e):\n", trial_barr, curr_barr);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of constraint violation = %23.16e  (current %23.16e):\n", trial_theta, curr_theta);

    // Check if point is acceptable w.r.t current iterate
    if (IsFtype(alpha_primal_test) && curr_theta <= theta_min_) {
      // Armijo condition for the barrier function has to be satisfied
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Checking Armijo Condition...\n");
      accept = ArmijoHolds(alpha_primal_test);
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Checking sufficient reduction...\n");
      accept = IsAcceptableToCurrentIterate(trial_barr, trial_theta);
    }

    if (!accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Failed...\n");
      return accept;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Succeeded...\n");
    }

    // Now check if that pair is acceptable to the filter
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Checking filter acceptability...\n");
    accept = IsAcceptableToCurrentFilter(trial_barr, trial_theta);
    if (!accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Failed...\n");
      return accept;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Succeeded...\n");
    }

    return accept;
  }

  bool FilterLineSearch::ArmijoHolds(Number alpha_primal_test)
  {
    return Compare_le(IpCq().trial_barrier_obj()-IpCq().curr_barrier_obj(),
                      eta_phi_*alpha_primal_test*IpCq().curr_gradBarrTDelta(),
                      IpCq().curr_barrier_obj());
  }

  Number FilterLineSearch::CalculateAlphaMin()
  {
    Number gBD = IpCq().curr_gradBarrTDelta();
    Number curr_theta = IpCq().curr_constraint_violation();
    Number alpha_min = gamma_theta_;

    if (gBD < 0) {
      alpha_min = Min( gamma_theta_,
                       gamma_phi_*curr_theta/(-gBD));
      if (curr_theta <= theta_min_) {
        alpha_min = Min( alpha_min,
                         delta_*pow(curr_theta,s_theta_)/pow(-gBD,s_phi_)
                       );
      }
    }

    return alpha_min_frac_*alpha_min;
  }

  bool FilterLineSearch::IsAcceptableToCurrentIterate(Number trial_barr, Number trial_theta) const
  {
    DBG_START_METH("FilterLineSearch::IsAcceptableToCurrentIterate",
                   dbg_verbosity);
    Number curr_barr = IpCq().curr_barrier_obj();

    // Check if the barrier objective function is increasing to
    // rapidly (according to option obj_max_inc)
    if (trial_barr > curr_barr) {
      Number basval = 1.;
      if (fabs(curr_barr)>10.) {
        basval = log10(fabs(curr_barr));
      }
      if (log10(trial_barr-curr_barr)>obj_max_inc_*basval) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Rejecting trial point because barrier objective function increasing too rapidly (from %27.15e to %27.15e)\n",curr_barr,trial_barr);
        return false;
      }
    }

    Number curr_theta = IpCq().curr_constraint_violation();
    DBG_PRINT((1,"trial_barr  = %e curr_barr  = %e\n", trial_barr, curr_barr));
    DBG_PRINT((1,"trial_theta = %e curr_theta = %e\n", trial_theta, curr_theta));
    return (Compare_le(trial_theta, (1.-gamma_theta_)*curr_theta, curr_theta)
            || Compare_le(trial_barr-curr_barr, -gamma_phi_*curr_theta, curr_barr));
  }

  bool FilterLineSearch::IsAcceptableToCurrentFilter(Number trial_barr, Number trial_theta) const
  {
    return filter_.Acceptable(trial_barr, trial_theta);
  }

  bool FilterLineSearch::Compare_le(Number lhs, Number rhs, Number BasVal)
  {
    DBG_START_FUN("FilterLineSearch::Compare_le",
                  dbg_verbosity);
    DBG_PRINT((1,"lhs = %27.16e rhs = %27.16e  BasVal = %27.16e\n",lhs,rhs,BasVal));
    // ToDo: Comparison based on machine precision
    return (lhs - rhs <= 1e-15*fabs(BasVal));
  }

  void FilterLineSearch::Reset()
  {
    DBG_START_FUN("FilterLineSearch::Reset", dbg_verbosity);
    filter_.Clear();
  }

  bool
  FilterLineSearch::TrySecondOrderCorrection(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<const Vector>& actual_delta_x,
    SmartPtr<const Vector>& actual_delta_s,
    SmartPtr<const Vector>& actual_delta_y_c,
    SmartPtr<const Vector>& actual_delta_y_d,
    SmartPtr<const Vector>& actual_delta_z_L,
    SmartPtr<const Vector>& actual_delta_z_U,
    SmartPtr<const Vector>& actual_delta_v_L,
    SmartPtr<const Vector>& actual_delta_v_U)
  {
    DBG_START_METH("FilterLineSearch::TrySecondOrderCorrection",
                   dbg_verbosity);

    bool accept = false;
    Index count_soc = 0;

    Number theta_soc_old = 0.;
    Number theta_curr = IpCq().curr_constraint_violation();
    Number theta_trial = 0.;
    Number alpha_primal_soc = alpha_primal;

    SmartPtr<Vector> c_soc = IpCq().curr_c()->MakeNew();
    SmartPtr<Vector> dms_soc = IpCq().curr_d_minus_s()->MakeNew();
    c_soc->Copy(*IpCq().curr_c());
    dms_soc->Copy(*IpCq().curr_d_minus_s());
    while (count_soc<max_soc_ &&
           theta_trial<=kappa_soc_*theta_soc_old &&
           !accept) {
      if (count_soc==0) {
        theta_soc_old = theta_curr;
      }
      else {
        theta_soc_old = theta_trial;
      }

      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Trying second order correction number %d\n",
                     count_soc+1);

      // Compute SOC constraint violation
      c_soc->Scal(alpha_primal_soc);
      dms_soc->Scal(alpha_primal_soc);
      c_soc->Axpy(1.0, *IpCq().trial_c());
      dms_soc->Axpy(1.0, *IpCq().trial_d_minus_s());

      // Compute the SOC search direction
      SmartPtr<Vector> delta_soc_x = actual_delta_x->MakeNew();
      SmartPtr<Vector> delta_soc_s = actual_delta_s->MakeNew();
      SmartPtr<Vector> delta_soc_y_c = actual_delta_y_c->MakeNew();
      SmartPtr<Vector> delta_soc_y_d = actual_delta_y_d->MakeNew();
      SmartPtr<Vector> delta_soc_z_L = actual_delta_z_L->MakeNew();
      SmartPtr<Vector> delta_soc_z_U = actual_delta_z_U->MakeNew();
      SmartPtr<Vector> delta_soc_v_L = actual_delta_v_L->MakeNew();
      SmartPtr<Vector> delta_soc_v_U = actual_delta_v_U->MakeNew();

      SmartPtr<const Vector> rhs_grad_lag_x
      = IpCq().curr_grad_lag_x();
      SmartPtr<const Vector> rhs_grad_lag_s
      = IpCq().curr_grad_lag_s();
      SmartPtr<const Vector> rhs_rel_compl_x_L
      = IpCq().curr_relaxed_compl_x_L();
      SmartPtr<const Vector> rhs_rel_compl_x_U
      = IpCq().curr_relaxed_compl_x_U();
      SmartPtr<const Vector> rhs_rel_compl_s_L
      = IpCq().curr_relaxed_compl_s_L();
      SmartPtr<const Vector> rhs_rel_compl_s_U
      = IpCq().curr_relaxed_compl_s_U();
      pd_solver_->Solve(-1.0, 0.0,
                        *rhs_grad_lag_x,
                        *rhs_grad_lag_s,
                        *c_soc,
                        *dms_soc,
                        *rhs_rel_compl_x_L,
                        *rhs_rel_compl_x_U,
                        *rhs_rel_compl_s_L,
                        *rhs_rel_compl_s_U,
                        *delta_soc_x,
                        *delta_soc_s,
                        *delta_soc_y_c,
                        *delta_soc_y_d,
                        *delta_soc_z_L,
                        *delta_soc_z_U,
                        *delta_soc_v_L,
                        *delta_soc_v_U,
                        true);

      // Compute step size
      alpha_primal_soc =
        IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                        *delta_soc_x,
                                        *delta_soc_s);

      // Check if trial point is acceptable
      try {
        // Compute the primal trial point
        IpData().SetTrialPrimalVariablesFromStep(alpha_primal_soc,
            *delta_soc_x,
            *delta_soc_s);

        // in acceptance tests, use original step size!
        accept = CheckAcceptabilityOfTrialPoint(alpha_primal_test);
      }
      catch(IpoptNLP::Eval_Error& e) {
        e.ReportException(Jnlst());
        Jnlst().Printf(J_WARNING, J_MAIN, "Warning: SOC step rejected due to evaluation error\n");
        accept = false;
      }

      if (accept) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Second order correction step accepted with %d corrections.\n", count_soc+1);
        // Accept all SOC quantities
        alpha_primal = alpha_primal_soc;
        actual_delta_x = ConstPtr(delta_soc_x);
        actual_delta_s = ConstPtr(delta_soc_s);
        actual_delta_y_c = ConstPtr(delta_soc_y_c);
        actual_delta_y_d = ConstPtr(delta_soc_y_d);
        actual_delta_z_L = ConstPtr(delta_soc_z_L);
        actual_delta_z_U = ConstPtr(delta_soc_z_U);
        actual_delta_v_L = ConstPtr(delta_soc_v_L);
        actual_delta_v_U = ConstPtr(delta_soc_v_U);
      }
      else {
        count_soc++;
        theta_trial = IpCq().trial_constraint_violation();
      }
    }
    return accept;
  }

  bool
  FilterLineSearch::TryCorrector(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<const Vector>& actual_delta_x,
    SmartPtr<const Vector>& actual_delta_s,
    SmartPtr<const Vector>& actual_delta_y_c,
    SmartPtr<const Vector>& actual_delta_y_d,
    SmartPtr<const Vector>& actual_delta_z_L,
    SmartPtr<const Vector>& actual_delta_z_U,
    SmartPtr<const Vector>& actual_delta_v_L,
    SmartPtr<const Vector>& actual_delta_v_U)
  {
    DBG_START_METH("FilterLineSearch::TryCorrector",
                   dbg_verbosity);

    bool accept = false;

    // Compute the corrector step based on corrector_type parameter
    SmartPtr<Vector> delta_corr_x = actual_delta_x->MakeNew();
    SmartPtr<Vector> delta_corr_s = actual_delta_s->MakeNew();
    SmartPtr<Vector> delta_corr_y_c = actual_delta_y_c->MakeNew();
    SmartPtr<Vector> delta_corr_y_d = actual_delta_y_d->MakeNew();
    SmartPtr<Vector> delta_corr_z_L = actual_delta_z_L->MakeNew();
    SmartPtr<Vector> delta_corr_z_U = actual_delta_z_U->MakeNew();
    SmartPtr<Vector> delta_corr_v_L = actual_delta_v_L->MakeNew();
    SmartPtr<Vector> delta_corr_v_U = actual_delta_v_U->MakeNew();

    switch (corrector_type_) {
      case 1 : {
        // 1: Standard MPC corrector

        // 	// I think it doesn't make sense to try the MPC corrector step
        // 	// in the fixed mu mode, since we are not heading to mu=0
        // 	if (!IpData().FreeMuMode()) {
        // 	  return false;
        // 	}

        // ToDo: For now, recompute the affine scaling step.  Later we
        // have to find a way so that it doesn't have to be recomputed
        // of it has been obtained as a probing step
        SmartPtr<Vector> delta_aff_x = actual_delta_x->MakeNew();
        SmartPtr<Vector> delta_aff_s = actual_delta_s->MakeNew();
        SmartPtr<Vector> delta_aff_y_c = actual_delta_y_c->MakeNew();
        SmartPtr<Vector> delta_aff_y_d = actual_delta_y_d->MakeNew();
        SmartPtr<Vector> delta_aff_z_L = actual_delta_z_L->MakeNew();
        SmartPtr<Vector> delta_aff_z_U = actual_delta_z_U->MakeNew();
        SmartPtr<Vector> delta_aff_v_L = actual_delta_v_L->MakeNew();
        SmartPtr<Vector> delta_aff_v_U = actual_delta_v_U->MakeNew();

        // Now solve the primal-dual system to get the step
        pd_solver_->Solve(-1.0, 0.0,
                          *IpCq().curr_grad_lag_x(),
                          *IpCq().curr_grad_lag_s(),
                          *IpCq().curr_c(),
                          *IpCq().curr_d_minus_s(),
                          *IpCq().curr_compl_x_L(),
                          *IpCq().curr_compl_x_U(),
                          *IpCq().curr_compl_s_L(),
                          *IpCq().curr_compl_s_U(),
                          *delta_aff_x,
                          *delta_aff_s,
                          *delta_aff_y_c,
                          *delta_aff_y_d,
                          *delta_aff_z_L,
                          *delta_aff_z_U,
                          *delta_aff_v_L,
                          *delta_aff_v_U,
                          true);

        DBG_PRINT_VECTOR(2, "delta_aff_x", *delta_aff_x);
        DBG_PRINT_VECTOR(2, "delta_aff_s", *delta_aff_s);
        DBG_PRINT_VECTOR(2, "delta_aff_y_c", *delta_aff_y_c);
        DBG_PRINT_VECTOR(2, "delta_aff_y_d", *delta_aff_y_d);
        DBG_PRINT_VECTOR(2, "delta_aff_z_L", *delta_aff_z_L);
        DBG_PRINT_VECTOR(2, "delta_aff_z_U", *delta_aff_z_U);
        DBG_PRINT_VECTOR(2, "delta_aff_v_L", *delta_aff_v_L);
        DBG_PRINT_VECTOR(2, "delta_aff_v_U", *delta_aff_v_U);

        delta_corr_x->Copy(*actual_delta_x);
        delta_corr_s->Copy(*actual_delta_s);
        delta_corr_y_c->Copy(*actual_delta_y_c);
        delta_corr_y_d->Copy(*actual_delta_y_d);
        delta_corr_z_L->Copy(*actual_delta_z_L);
        delta_corr_z_U->Copy(*actual_delta_z_U);
        delta_corr_v_L->Copy(*actual_delta_v_L);
        delta_corr_v_U->Copy(*actual_delta_v_U);

        SmartPtr<Vector> rhs_x = actual_delta_x->MakeNew();
        SmartPtr<Vector> rhs_s = actual_delta_s->MakeNew();
        SmartPtr<Vector> rhs_c = actual_delta_y_c->MakeNew();
        SmartPtr<Vector> rhs_d = actual_delta_y_d->MakeNew();
        SmartPtr<Vector> rhs_compl_x_L = actual_delta_z_L->MakeNew();
        SmartPtr<Vector> rhs_compl_x_U = actual_delta_z_U->MakeNew();
        SmartPtr<Vector> rhs_compl_s_L = actual_delta_v_L->MakeNew();
        SmartPtr<Vector> rhs_compl_s_U = actual_delta_v_U->MakeNew();

        rhs_x->Set(0.);
        rhs_s->Set(0.);
        rhs_c->Set(0.);
        rhs_d->Set(0.);
        IpNLP().Px_L()->TransMultVector(-1., *delta_aff_x, 0., *rhs_compl_x_L);
        rhs_compl_x_L->ElementWiseMultiply(*delta_aff_z_L);
        IpNLP().Px_U()->TransMultVector(1., *delta_aff_x, 0., *rhs_compl_x_U);
        rhs_compl_x_U->ElementWiseMultiply(*delta_aff_z_U);
        IpNLP().Pd_L()->TransMultVector(-1., *delta_aff_s, 0., *rhs_compl_s_L);
        rhs_compl_s_L->ElementWiseMultiply(*delta_aff_v_L);
        IpNLP().Pd_U()->TransMultVector(1., *delta_aff_s, 0., *rhs_compl_s_U);
        rhs_compl_s_U->ElementWiseMultiply(*delta_aff_v_U);

        pd_solver_->Solve(1.0, 1.0,
                          *rhs_x,
                          *rhs_s,
                          *rhs_c,
                          *rhs_d,
                          *rhs_compl_x_L,
                          *rhs_compl_x_U,
                          *rhs_compl_s_L,
                          *rhs_compl_s_U,
                          *delta_corr_x,
                          *delta_corr_s,
                          *delta_corr_y_c,
                          *delta_corr_y_d,
                          *delta_corr_z_L,
                          *delta_corr_z_U,
                          *delta_corr_v_L,
                          *delta_corr_v_U,
                          true);

        DBG_PRINT_VECTOR(2, "delta_corr_x", *delta_corr_x);
        DBG_PRINT_VECTOR(2, "delta_corr_s", *delta_corr_s);
        DBG_PRINT_VECTOR(2, "delta_corr_y_c", *delta_corr_y_c);
        DBG_PRINT_VECTOR(2, "delta_corr_y_d", *delta_corr_y_d);
        DBG_PRINT_VECTOR(2, "delta_corr_z_L", *delta_corr_z_L);
        DBG_PRINT_VECTOR(2, "delta_corr_z_U", *delta_corr_z_U);
        DBG_PRINT_VECTOR(2, "delta_corr_v_L", *delta_corr_v_L);
        DBG_PRINT_VECTOR(2, "delta_corr_v_U", *delta_corr_v_U);
      }
      break;
      case 2 : {
        // 2: Second order correction for primal-dual step to
        // primal-dual mu

        delta_corr_x->Copy(*actual_delta_x);
        delta_corr_s->Copy(*actual_delta_s);
        delta_corr_y_c->Copy(*actual_delta_y_c);
        delta_corr_y_d->Copy(*actual_delta_y_d);
        delta_corr_z_L->Copy(*actual_delta_z_L);
        delta_corr_z_U->Copy(*actual_delta_z_U);
        delta_corr_v_L->Copy(*actual_delta_v_L);
        delta_corr_v_U->Copy(*actual_delta_v_U);

        SmartPtr<Vector> rhs_x = actual_delta_x->MakeNew();
        SmartPtr<Vector> rhs_s = actual_delta_s->MakeNew();
        SmartPtr<Vector> rhs_c = actual_delta_y_c->MakeNew();
        SmartPtr<Vector> rhs_d = actual_delta_y_d->MakeNew();
        SmartPtr<Vector> rhs_compl_x_L = actual_delta_z_L->MakeNew();
        SmartPtr<Vector> rhs_compl_x_U = actual_delta_z_U->MakeNew();
        SmartPtr<Vector> rhs_compl_s_L = actual_delta_v_L->MakeNew();
        SmartPtr<Vector> rhs_compl_s_U = actual_delta_v_U->MakeNew();

        rhs_x->Set(0.);
        rhs_s->Set(0.);
        rhs_c->Set(0.);
        rhs_d->Set(0.);

        Number mu = IpData().curr_mu();
        SmartPtr<Vector> tmp;

        rhs_compl_x_L->Copy(*IpCq().curr_slack_x_L());
        IpNLP().Px_L()->TransMultVector(-1., *actual_delta_x,
                                        -1., *rhs_compl_x_L);
        tmp = actual_delta_z_L->MakeNew();
        tmp->Copy(*IpData().curr_z_L());
        tmp->Axpy(1., *actual_delta_z_L);
        rhs_compl_x_L->ElementWiseMultiply(*tmp);
        rhs_compl_x_L->AddScalar(mu);

        rhs_compl_x_U->Copy(*IpCq().curr_slack_x_U());
        IpNLP().Px_U()->TransMultVector(1., *actual_delta_x,
                                        -1., *rhs_compl_x_U);
        tmp = actual_delta_z_U->MakeNew();
        tmp->Copy(*IpData().curr_z_U());
        tmp->Axpy(1., *actual_delta_z_U);
        rhs_compl_x_U->ElementWiseMultiply(*tmp);
        rhs_compl_x_U->AddScalar(mu);

        rhs_compl_s_L->Copy(*IpCq().curr_slack_s_L());
        IpNLP().Pd_L()->TransMultVector(-1., *actual_delta_s,
                                        -1., *rhs_compl_s_L);
        tmp = actual_delta_v_L->MakeNew();
        tmp->Copy(*IpData().curr_v_L());
        tmp->Axpy(1., *actual_delta_v_L);
        rhs_compl_s_L->ElementWiseMultiply(*tmp);
        rhs_compl_s_L->AddScalar(mu);

        rhs_compl_s_U->Copy(*IpCq().curr_slack_s_U());
        IpNLP().Pd_U()->TransMultVector(1., *actual_delta_s,
                                        -1., *rhs_compl_s_U);
        tmp = actual_delta_v_U->MakeNew();
        tmp->Copy(*IpData().curr_v_U());
        tmp->Axpy(1., *actual_delta_v_U);
        rhs_compl_s_U->ElementWiseMultiply(*tmp);
        rhs_compl_s_U->AddScalar(mu);

        DBG_PRINT_VECTOR(2, "rhs_compl_x_L", *rhs_compl_x_L);
        DBG_PRINT_VECTOR(2, "rhs_compl_x_U", *rhs_compl_x_U);
        DBG_PRINT_VECTOR(2, "rhs_compl_s_L", *rhs_compl_s_L);
        DBG_PRINT_VECTOR(2, "rhs_compl_s_U", *rhs_compl_s_U);

        pd_solver_->Solve(1.0, 1.0,
                          *rhs_x,
                          *rhs_s,
                          *rhs_c,
                          *rhs_d,
                          *rhs_compl_x_L,
                          *rhs_compl_x_U,
                          *rhs_compl_s_L,
                          *rhs_compl_s_U,
                          *delta_corr_x,
                          *delta_corr_s,
                          *delta_corr_y_c,
                          *delta_corr_y_d,
                          *delta_corr_z_L,
                          *delta_corr_z_U,
                          *delta_corr_v_L,
                          *delta_corr_v_U,
                          true);

        DBG_PRINT_VECTOR(2, "delta_corr_x", *delta_corr_x);
        DBG_PRINT_VECTOR(2, "delta_corr_s", *delta_corr_s);
        DBG_PRINT_VECTOR(2, "delta_corr_y_c", *delta_corr_y_c);
        DBG_PRINT_VECTOR(2, "delta_corr_y_d", *delta_corr_y_d);
        DBG_PRINT_VECTOR(2, "delta_corr_z_L", *delta_corr_z_L);
        DBG_PRINT_VECTOR(2, "delta_corr_z_U", *delta_corr_z_U);
        DBG_PRINT_VECTOR(2, "delta_corr_v_L", *delta_corr_v_L);
        DBG_PRINT_VECTOR(2, "delta_corr_v_U", *delta_corr_v_U);
      }
      break;
      default:
      DBG_ASSERT("Unknown corrector_type value.");
    }

    // Compute step size
    Number alpha_primal_corr =
      IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                      *delta_corr_x,
                                      *delta_corr_s);
    // Compute the primal trial point
    IpData().SetTrialPrimalVariablesFromStep(alpha_primal_corr,
        *delta_corr_x,
        *delta_corr_s);

    // Check if we want to not even try the filter criterion
    Number alpha_dual_max =
      IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                    *delta_corr_z_L, *delta_corr_z_U,
                                    *delta_corr_v_L, *delta_corr_v_U);

    IpData().SetTrialBoundMutlipliersFromStep(alpha_dual_max,
        *delta_corr_z_L, *delta_corr_z_U,
        *delta_corr_v_L, *delta_corr_v_U);

    Number trial_avrg_compl = IpCq().trial_avrg_compl();
    Number curr_avrg_compl = IpCq().curr_avrg_compl();
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "avrg_compl(curr) = %e, avrg_compl(trial) = %e\n",
                   curr_avrg_compl, trial_avrg_compl);
    if (corrector_type_==1 &&
        trial_avrg_compl>=corrector_compl_avrg_red_fact_*curr_avrg_compl) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Rejecting corrector step, because trial complementarity is too large.\n" );
      return false;
    }

    // Check if trial point is acceptable
    try {
      // in acceptance tests, use original step size!
      accept = CheckAcceptabilityOfTrialPoint(alpha_primal_test);
    }
    catch(IpoptNLP::Eval_Error& e) {
      e.ReportException(Jnlst());
      Jnlst().Printf(J_WARNING, J_MAIN, "Warning: Corrector step rejected due to evaluation error\n");
      accept = false;
    }

    if (accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Corrector step accepted with alpha_primal = %e\n",
                     alpha_primal_corr);
      // Accept all SOC quantities
      alpha_primal = alpha_primal_corr;
      actual_delta_x = ConstPtr(delta_corr_x);
      actual_delta_s = ConstPtr(delta_corr_s);
      actual_delta_y_c = ConstPtr(delta_corr_y_c);
      actual_delta_y_d = ConstPtr(delta_corr_y_d);
      actual_delta_z_L = ConstPtr(delta_corr_z_L);
      actual_delta_z_U = ConstPtr(delta_corr_z_U);
      actual_delta_v_L = ConstPtr(delta_corr_v_L);
      actual_delta_v_U = ConstPtr(delta_corr_v_U);

      if (Jnlst().ProduceOutput(J_MOREVECTOR, J_MAIN)) {
        Jnlst().Printf(J_MOREVECTOR, J_MAIN,
                       "*** Accepted corrector for Iteration: %d\n",
                       IpData().iter_count());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_x", *delta_corr_x);
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_s", *delta_corr_s);
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_y_c", *delta_corr_y_c);
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_y_d", *delta_corr_y_d);
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_z_L", *delta_corr_z_L);
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_z_U", *delta_corr_z_U);
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_v_L", *delta_corr_v_L);
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr_v_U", *delta_corr_v_U);
      }
    }

    return accept;
  }

  void
  FilterLineSearch::PerformMagicStep()
  {
    DBG_START_METH("FilterLineSearch::PerformMagicStep",
                   dbg_verbosity);

    DBG_PRINT((1,"Incoming barr = %e and constrviol %e\n",
               IpCq().trial_barrier_obj(),
               IpCq().trial_constraint_violation()));
    DBG_PRINT_VECTOR(2, "s in", *IpData().trial_s());
    DBG_PRINT_VECTOR(2, "d minus s in", *IpCq().trial_d_minus_s());
    DBG_PRINT_VECTOR(2, "slack_s_L in", *IpCq().trial_slack_s_L());
    DBG_PRINT_VECTOR(2, "slack_s_U in", *IpCq().trial_slack_s_U());

    SmartPtr<const Vector> d_L = IpNLP().d_L();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<Vector> delta_s_magic_L = d_L->MakeNew();
    delta_s_magic_L->Set(0.);
    SmartPtr<Vector> tmp = d_L->MakeNew();
    Pd_L->TransMultVector(1., *IpCq().trial_d_minus_s(), 0., *tmp);
    delta_s_magic_L->ElementWiseMax(*tmp);

    SmartPtr<const Vector> d_U = IpNLP().d_U();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<Vector> delta_s_magic_U = d_U->MakeNew();
    delta_s_magic_U->Set(0.);
    tmp = d_U->MakeNew();
    Pd_U->TransMultVector(1., *IpCq().trial_d_minus_s(), 0., *tmp);
    delta_s_magic_U->ElementWiseMin(*tmp);

    SmartPtr<Vector> delta_s_magic = IpData().trial_s()->MakeNew();
    Pd_L->MultVector(1., *delta_s_magic_L, 0., *delta_s_magic);
    Pd_U->MultVector(1., *delta_s_magic_U, 1., *delta_s_magic);
    delta_s_magic_L = NULL; // free memory
    delta_s_magic_U = NULL; // free memory

    // Now find those entries with both lower and upper bounds, there
    // the step is too large
    // ToDo this should only be done if there are inequality
    // constraints with two bounds
    // also this can be done in a smaller space (d_L or d_U whichever
    // is smaller)
    tmp = delta_s_magic->MakeNew();
    tmp->Copy(*IpData().trial_s());
    Pd_L->MultVector(1., *d_L, -2., *tmp);
    Pd_U->MultVector(1., *d_U, 1., *tmp);
    SmartPtr<Vector> tmp2 = tmp->MakeNew();
    tmp2->Copy(*tmp);
    tmp2->ElementWiseAbs();
    tmp->Axpy(-2., *delta_s_magic);
    tmp->ElementWiseAbs();
    // now, tmp2 = |d_L + d_u - 2*s| and tmp = |d_L + d_u - 2*(s+Delta s)|
    // we want to throw out those for which tmp2 > tmp
    tmp->Axpy(-1., *tmp2);
    tmp->ElementWiseSgn();
    tmp2->Set(0.);
    tmp2->ElementWiseMax(*tmp);
    tmp = d_L->MakeNew();
    Pd_L->TransMultVector(1., *tmp2, 0., *tmp);
    Pd_L->MultVector(1., *tmp, 0., *tmp2);
    tmp = d_U->MakeNew();
    Pd_U->TransMultVector(1., *tmp2, 0., *tmp);
    Pd_U->MultVector(1., *tmp, 0., *tmp2);
    DBG_PRINT_VECTOR(2, "tmp indicator", *tmp2)
    // tmp2 now is one for those entries with both bounds, for which
    // no step should be taken

    tmp = delta_s_magic->MakeNew();
    tmp->Copy(*delta_s_magic);
    tmp->ElementWiseMultiply(*tmp2);
    delta_s_magic->Axpy(-1., *tmp);

    Number delta_s_magic_max = delta_s_magic->Amax();
    Number mach_eps = std::numeric_limits<Number>::epsilon();
    if (delta_s_magic_max>0.) {
      if (delta_s_magic_max > 10*mach_eps*IpData().trial_s()->Amax()) {
        IpData().Append_info_string("M");
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Magic step with max-norm %.6e taken.\n", delta_s_magic->Amax());
        Jnlst().PrintVector(J_MOREVECTOR, J_LINE_SEARCH,
                            "delta_s_magic", *delta_s_magic);
      }

      // now finally compute the new overall slacks
      delta_s_magic->Axpy(1., *IpData().trial_s());
      IpData().SetTrialSVariables(*delta_s_magic);
    }

    DBG_PRINT((1,"Outgoing barr = %e and constrviol %e\n", IpCq().trial_barrier_obj(), IpCq().trial_constraint_violation()));
    DBG_PRINT_VECTOR(2, "s out", *IpData().trial_s());
    DBG_PRINT_VECTOR(2, "d minus s out", *IpCq().trial_d_minus_s());
    DBG_PRINT_VECTOR(2, "slack_s_L out", *IpCq().trial_slack_s_L());
    DBG_PRINT_VECTOR(2, "slack_s_U out", *IpCq().trial_slack_s_U());
  }

} // namespace Ipopt
