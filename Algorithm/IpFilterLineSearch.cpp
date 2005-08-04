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
#include "IpAlgTypes.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
# include <ctype.h>
#else
# include <cmath>
# include <cctype>
#endif

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  DefineIpoptType(FilterLineSearch);

  FilterLineSearch::FilterLineSearch(const SmartPtr<RestorationPhase>& resto_phase,
                                     const SmartPtr<PDSystemSolver>& pd_solver,
                                     const SmartPtr<ConvergenceCheck>& conv_check)
      :
      LineSearch(),
      theta_max_(-1.0),
      theta_min_(-1.0),
      filter_(2),
      resto_phase_(resto_phase),
      pd_solver_(pd_solver),
      conv_check_(conv_check)
  {
    DBG_START_FUN("FilterLineSearch::FilterLineSearch",
                  dbg_verbosity);
  }

  FilterLineSearch::~FilterLineSearch()
  {
    DBG_START_FUN("FilterLineSearch::~FilterLineSearch()",
                  dbg_verbosity);
  }

  void FilterLineSearch::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "theta_max_fact",
      "Determines upper bound for constraint violation in the filter.",
      0.0, true, 1e4,
      "The algorithmic parameter theta_max is determined as theta_max_fact "
      "times the maximum of 1 and the constraint violation at initial point.  "
      "Any point with a constraint violation larger than theta_max is "
      "unacceptable to the filter (see Eqn. (21) in implementation paper).");
    roptions->AddLowerBoundedNumberOption(
      "theta_min_fact",
      "Determines constraint violation threshold in switching rule.",
      0.0, true, 1e-4,
      "The algorithmic parameter theta_min is determined as theta_max_fact "
      "times the maximum of 1 and the constraint violation at initial point.  "
      "The switching rules treats an iteration as h-type iteration whenever "
      "the current constraint violation is larger than theta_min (see "
      "paragraph before Eqn. (19) in implementation paper).");
    roptions->AddBoundedNumberOption(
      "eta_phi",
      "Relaxation factor in the Armijo condition.",
      0.0, true, 0.5, true, 1e-8,
      "(See Eqn. (20) in implementation paper)");
    roptions->AddLowerBoundedNumberOption(
      "delta", "Multiplier for constraint violation in switching rule.",
      0.0, true, 1.0,
      "(See Eqn. (19) in implementation paper)");
    roptions->AddLowerBoundedNumberOption(
      "s_phi",
      "Exponent for linear barrier function model in switching rule.",
      1.0, true, 2.3,
      "(See Eqn. (19) in implementation paper)");
    roptions->AddLowerBoundedNumberOption(
      "s_theta",
      "Exponent for current constraint violation in switching rule.",
      1.0, true, 1.1,
      "(See Eqn. (19) in implementation paper)");
    roptions->AddBoundedNumberOption(
      "gamma_phi",
      "Relaxation factor in filter margin for barrier function.",
      0.0, true, 1.0, true, 1e-8,
      "(See Eqn. (18a) in implementation paper)");
    roptions->AddBoundedNumberOption(
      "gamma_theta",
      "Relaxation factor in filter margin for constraint violation.",
      0.0, true, 1.0, true, 1e-5,
      "(See Eqn. (18b) in implementation paper)");
    roptions->AddBoundedNumberOption(
      "alpha_min_frac",
      "Safety factor for minimal step size (switch to restoration phase).",
      0.0, true, 1.0, true, 0.05,
      "(This is gamma_alpha in Eqn. (20) in implementation paper)");
    roptions->AddBoundedNumberOption(
      "alpha_red_factor",
      "Fractional reduction of trial step size in the backtracking line search.",
      0.0, true, 1.0, true, 0.5,
      "Determines the fraction by how much the trial step size is reduced in "
      "every step of the backtracking line search.");
    roptions->AddLowerBoundedIntegerOption(
      "max_soc",
      "Maximal number of second order correction trial steps.",
      0, 4,
      "Determines the maximal number of second order correction trial steps "
      "that should be performed.  Choosing 0 disables the second order "
      "corrections. (This is p^{max} of Step A-5.9 of "
      "Algorithm A in implementation paper.)");
    roptions->AddLowerBoundedNumberOption(
      "kappa_soc",
      "Factor in sufficient reduction rule for second order correction.",
      0.0, true, 0.99,
      "Determines by how much a second order correction step must reduce the "
      "constraint violation so that further correction steps are attempted.  "
      "(See Step A-5.9 of Algorithm A in implementation paper.)");
    roptions->AddLowerBoundedNumberOption(
      "obj_max_inc",
      "Determines upper bound on acceptable increase of barrier objective function.",
      1.0, true, 5.0,
      "A trial point leading to more orders of magnitude increase in the "
      "barrier objective function are rejected.");
    roptions->AddStringOption2(
      "magic_steps",
      "Enables magic steps.",
      "no",
      "no", "don't take magic steps",
      "yes", "take magic steps",
      "DOESN'T REALLY WORK YET!");
    roptions->AddStringOption3(
      "corrector_type",
      "Type of corrector steps.",
      "none",
      "none", "no corrector",
      "affine", "corrector step towards mu=0",
      "primal-dual", "corrector step towards current mu",
      "Determines what kind of corrector steps should be tried.");

    roptions->AddStringOption2(
      "skip_corr_if_neg_curv",
      "Skip the corrector step in negative curvature iteration.",
      "yes",
      "no", "don't skip",
      "yes", "skip",
      "The corrector step is not tried if during the computation of "
      "the search direction in the current iteration negative curvature has "
      "been encountered.");

    roptions->AddStringOption2(
      "skip_corr_in_monotone_mode",
      "Skip the corrector step during monotone barrier parameter mode.",
      "yes",
      "no", "don't skip",
      "yes", "skip",
      "The corrector step is not tried if the algorithm is currently in the "
      "monotone mode (see also option \"barrier_strategy\").");

    roptions->AddStringOption2(
      "accept_every_trial_step",
      "always accept the frist trial step",
      "no",
      "no", "don't arbitrarily accept the full step",
      "yes", "always accept the full step",
      "Setting this option to \"yes\" essentially disables the line search "
      "and makes the algorithm take aggressive steps.");

    roptions->AddStringOption7(
      "alpha_for_y",
      "Step size for constraint multipliers.",
      "primal",
      "primal", "use primal step size",
      "bound_mult", "use step size for the bound multipliers",
      "min", "use the min of primal and bound multipliers",
      "max", "use the max of primal and bound multipliers",
      "full", "take a full step of size one",
      "min_dual_infeas", "choose step size minimizing new dual infeasibility",
      "safe_min_dual_infeas", "like \"min_dual_infeas\", but safeguarded by \"min\" and \"max\"",
      "Determines which step size (alpha_y) should be used to update the "
      "constraint multipliers.");

    roptions->AddLowerBoundedNumberOption(
      "corrector_compl_avrg_red_fact",
      "Complementarity tolerance factor for accepting corrector step",
      0.0, true, 1.0,
      "Determines the factor by which complementarity is allowed to increase "
      "for a corrector step to be accepted.");

    roptions->AddStringOption2(
      "expect_infeasible_problem",
      "Enable heuristics to quickly detect an infeasible problem.",
      "no",
      "no", "the problem probably be feasible",
      "yes", "the problem has a good chance to be infeasible",
      "This options is meant to activate heuristics that may speed up the "
      "infeasibility determination if you expect the problem to be "
      "infeasible.  In the filter line search procedure, the restoration "
      "phase is called more qucikly than usually, and more reduction in "
      "the constraint violation is enforced. If the problem is square, this "
      "is enabled automatically.");
    roptions->AddLowerBoundedNumberOption(
      "expect_infeasible_problem_ctol",
      "Threshold for disabling \"expect_infeasible_problem\" option",
      0.0, false, 1e-3,
      "If the constraint violation becomes small than this threshold, "
      "the \"expect_infeasible_problem\" heuristics in the filter line "
      "search will are disabled. If the problem is square, this is set to 0.");
    roptions->AddLowerBoundedNumberOption(
      "soft_resto_pderror_reduction_factor",
      "Required reduction in primal-dual error in soft restoration phase.",
      0.0, false, (1.0 - 1e-4),
      "For the soft restoration phase (which attempts to reduce the "
      "primal-dual error with regular steps), this indicates by which "
      "factor the primal-dual error has to be reduced in order to continue "
      "with the soft restoration phase. If the regular primal-dual step, "
      "damped onl to satisfty the fraction-to-the-boundary rule, is not "
      "decreasing the error by this factor, then the regular restoration "
      "phase is called.  Choosing \"0\" here disables the soft "
      "restoration phase.");
    roptions->AddStringOption2(
      "start_with_resto",
      "Tells algorithm to switch to restoration phase in first iteration.",
      "no",
      "no", "don't force start in restoration phase",
      "yes", "force start in restoration phase",
      "Setting this option to yes forces the algorithm to switch to the "
      "restoration phase in the first iteration.  If the initial point "
      "is feasible, the algorithm will abort with a failure.");
    roptions->AddLowerBoundedNumberOption(
      "tiny_step_tol",
      "Tolerance for detecting numerically insignificant steps.",
      0.0, false, 10.0*std::numeric_limits<double>::epsilon(),
      "If the search direction in the primal variables (x and s) is, in "
      "relative terms for each component, less than this values, the "
      "algorithm accepts the full step without line search.  The default "
      "value is 10 times machine precision.");
    roptions->AddLowerBoundedIntegerOption(
      "watchdog_shortened_iter_trigger",
      "Number of shortened iterations that trigger the watchdog.",
      0, 10,
      "If the number of iterations in which the backtracking line search "
      "did not accept the first trial point exceedes this number, the "
      "watchdog procedure is activated.  Choosing \"0\" here disables the "
      "watchdog procedure.");
    roptions->AddLowerBoundedIntegerOption(
      "watchdog_trial_iter_max",
      "Maximal number of watchdog iterations.",
      1, 3,
      "Determines the number of trial iterations before the watchdog "
      "procedure is aborted and the algorithm returns to the stored point.");
  }

  bool FilterLineSearch::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    options.GetNumericValue("theta_max_fact", theta_max_fact_, prefix);
    options.GetNumericValue("theta_min_fact", theta_min_fact_, prefix);
    ASSERT_EXCEPTION(theta_min_fact_ < theta_max_fact_, OptionsList::OPTION_OUT_OF_RANGE,
                     "Option \"theta_min_fact\": This value must be larger than 0 and less than theta_max_fact.");
    options.GetNumericValue("eta_phi", eta_phi_, prefix);
    options.GetNumericValue("delta", delta_, prefix);
    options.GetNumericValue("s_phi", s_phi_, prefix);
    options.GetNumericValue("s_theta", s_theta_, prefix);
    options.GetNumericValue("gamma_phi", gamma_phi_, prefix);
    options.GetNumericValue("gamma_theta", gamma_theta_, prefix);
    options.GetNumericValue("alpha_min_frac", alpha_min_frac_, prefix);
    options.GetNumericValue("alpha_red_factor", alpha_red_factor_, prefix);
    options.GetIntegerValue("max_soc", max_soc_, prefix);
    if (max_soc_>0) {
      ASSERT_EXCEPTION(IsValid(pd_solver_), OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"max_soc\": This option is non-negative, but no linear solver for computing the SOC given to FilterLineSearch object.");
    }
    options.GetNumericValue("kappa_soc", kappa_soc_, prefix);
    options.GetNumericValue("obj_max_inc", obj_max_inc_, prefix);
    options.GetBoolValue("magic_steps", magic_steps_, prefix);
    Index enum_int;
    options.GetEnumValue("corrector_type", enum_int, prefix);
    corrector_type_ = CorrectorTypeEnum(enum_int);
    options.GetBoolValue("skip_corr_if_neg_curv", skip_corr_if_neg_curv_, prefix);
    options.GetBoolValue("skip_corr_in_monotone_mode", skip_corr_in_monotone_mode_, prefix);
    options.GetBoolValue("accept_every_trial_step", accept_every_trial_step_, prefix);
    options.GetEnumValue("alpha_for_y", enum_int, prefix);
    alpha_for_y_ = AlphaForYEnum(enum_int);
    options.GetNumericValue("corrector_compl_avrg_red_fact", corrector_compl_avrg_red_fact_, prefix);
    options.GetNumericValue("expect_infeasible_problem_ctol", expect_infeasible_problem_ctol_, prefix);
    options.GetBoolValue("expect_infeasible_problem", expect_infeasible_problem_, prefix);

    options.GetBoolValue("start_with_resto", start_with_resto_, prefix);

    bool retvalue = true;
    if (IsValid(resto_phase_)) {
      retvalue = resto_phase_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                          options, prefix);
    }

    options.GetNumericValue("soft_resto_pderror_reduction_factor",
                            soft_resto_pderror_reduction_factor_, prefix);
    options.GetNumericValue("tiny_step_tol", tiny_step_tol_, prefix);
    options.GetIntegerValue("watchdog_trial_iter_max", watchdog_trial_iter_max_, prefix);
    options.GetIntegerValue("watchdog_shortened_iter_trigger", watchdog_shortened_iter_trigger_, prefix);

    // ToDo decide if also the PDSystemSolver should be initialized here...

    rigorous_ = true;
    skipped_line_search_ = false;
    tiny_step_last_iteration_ = false;

    Reset();

    count_successive_shortened_steps_ = 0;

    acceptable_iterate_ = NULL;

    return retvalue;
  }

  void FilterLineSearch::FindAcceptableTrialPoint()
  {
    DBG_START_METH("FilterLineSearch::FindAcceptableTrialPoint",
                   dbg_verbosity);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "--> Starting filter line search in iteration %d <--\n",
                   IpData().iter_count());

    // If the problem is square, we want to enable the
    // expect_infeasible_problem option automatically so that the
    // restoration phase is entered soon
    if (IpCq().IsSquareProblem()) {
      expect_infeasible_problem_ = true;
      expect_infeasible_problem_ctol_ = 0.;
    }

    // Store current iterate if the optimality error is on acceptable
    // level to restored if things fail later
    if (CurrentIsAcceptable()) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Storing current iterate as backup acceptable point.\n");
      StoreAcceptablePoint();
    }

    // First assume that line search will find an acceptable trial point
    skipped_line_search_ = false;

    // Set the values for the reference point
    if (!in_watchdog_) {
      reference_theta_ = IpCq().curr_constraint_violation();
      reference_barr_ = IpCq().curr_barrier_obj();
      reference_gradBarrTDelta_ = IpCq().curr_gradBarrTDelta();
    }
    else {
      reference_theta_ = watchdog_theta_;
      reference_barr_ = watchdog_barr_;
      reference_gradBarrTDelta_ = watchdog_gradBarrTDelta_;
    }

    // Get the search directions (this will store the actual search
    // direction, possibly including higher order corrections)
    SmartPtr<IteratesVector> actual_delta = IpData().delta()->MakeNewContainer();

    bool goto_resto = false;
    if (actual_delta->Asum() == 0. ) {
      // In this case, a search direction could not be computed, and
      // we should immediately go to the restoration phase.  ToDo: Cue
      // off of a something else than the norm of the search direction
      goto_resto = true;
    }

    if (start_with_resto_) {
      // If the use requested to start with the restoration phase,
      // skip the lin e search and do exactly that.  Reset the flag so
      // that this happens only once.
      goto_resto = true;
      start_with_resto_= false;
    }

    bool accept = false;
    bool corr_taken = false;
    bool soc_taken = false;
    Index n_steps = 0;
    Number alpha_primal = 0.;

    // Check if search direction becomes too small
    // ToDo: move this into place independent of this particular line search?
    bool tiny_step = (!goto_resto && DetectTinyStep());

    if (in_watchdog_ && (goto_resto || tiny_step)) {
      // If the step could not be computed or is too small and the
      // watchdog is active, stop the watch dog and resume everything
      // from reference point
      StopWatchDog(actual_delta);
      goto_resto = false;
      tiny_step = false;
    }

    // Check if we want to wake up the watchdog
    if (watchdog_shortened_iter_trigger_ > 0 &&
        !in_watchdog_ && !goto_resto && !tiny_step &&
        !in_soft_resto_phase_ && !expect_infeasible_problem_ &&
        watchdog_shortened_iter_ >= watchdog_shortened_iter_trigger_) {
      StartWatchDog();
    }

    if (tiny_step) {
      alpha_primal =
        IpCq().curr_primal_frac_to_the_bound(IpData().curr_tau());
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Tiny step detected. Use step size alpha = %e unchecked\n",
                     alpha_primal);
      IpData().SetTrialPrimalVariablesFromStep(alpha_primal, *IpData().delta()->x(), *IpData().delta()->s());
      IpData().Set_info_ls_count(0);

      if (tiny_step_last_iteration_) {
        IpData().Set_info_alpha_primal_char('T');
        IpData().Set_tiny_step_flag(true);
      }
      else {
        IpData().Set_info_alpha_primal_char('t');
      }

      tiny_step_last_iteration_ = true;
      accept = true;
    }
    else {
      tiny_step_last_iteration_ = false;
    }

    if (!goto_resto && !tiny_step) {

      if (in_soft_resto_phase_) {
        bool satisfies_original_filter = false;
        // ToDo use tiny_step in TrySoftRestoStep?
        accept = TrySoftRestoStep(actual_delta,
                                  satisfies_original_filter);
        if (accept) {
          IpData().Set_info_alpha_primal_char('s');
          if (satisfies_original_filter) {
            in_soft_resto_phase_ = false;
            IpData().Set_info_alpha_primal_char('S');
          }
        }
      }
      else {
        bool done = false;
        bool skip_first_trial_point = false;
        bool evaluation_error;
        while (!done) {
          accept = DoBacktrackingLineSearch(skip_first_trial_point,
                                            alpha_primal,
                                            corr_taken,
                                            soc_taken,
                                            n_steps,
                                            evaluation_error,
                                            actual_delta);
          DBG_PRINT((1, "evaluation_error = %d\n", evaluation_error));
          if (in_watchdog_) {
            if (accept) {
              in_watchdog_ = false;
              IpData().Append_info_string("W");
              Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                             "Watch dog procedure successful!\n");
              done = true;
            }
            else {
              watchdog_trial_iter_++;
              if (evaluation_error ||
                  watchdog_trial_iter_ > watchdog_trial_iter_max_) {
                StopWatchDog(actual_delta);
                skip_first_trial_point = true;
              }
              else {
                done = true;
                accept = true;
              }
            }
          }
          else {
            done = true;
          }
        }
      } /* else: if (in_soft_resto_phase_) { */
    } /* if (!goto_resto && !tiny_step) { */

    // If line search has been aborted because the step size becomes too small,
    // go to the restoration phase
    if (!accept) {
      // If we are not asked to do a rigorous line search, do no call
      // the restoration phase.
      if (!rigorous_) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Skipping call of restoration phase...\n");
        skipped_line_search_=true;
      }
      else {
        // Check if we should start the soft restoration phase
        if (!in_soft_resto_phase_ && soft_resto_pderror_reduction_factor_>0.
            && !goto_resto && !expect_infeasible_problem_) {
          Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                         "--> Starting soft restoration phase <--\n");
          // Augment the filter with the current point
          AugmentFilter();

          // Try the current search direction for the soft restoration phase
          bool satisfies_original_filter;
          accept = TrySoftRestoStep(actual_delta,
                                    satisfies_original_filter);
          // If it has been accepted: If the original filter is also
          // satisfied, we can just take that step and continue with
          // the regular algorithm, otherwise we stay in the soft
          // restoration phase
          if (accept) {
            if (satisfies_original_filter) {
              IpData().Set_info_alpha_primal_char('S');
            }
            else {
              in_soft_resto_phase_ = true;
              IpData().Set_info_alpha_primal_char('s');
            }
          }
        }

        if (!accept) {
          if (!in_soft_resto_phase_) {
            // Augment the filter with the current point if we are
            // already in the soft restoration phase, this has been
            // done earlier
            AugmentFilter();
          }
          if (CurrentIsAcceptable()) {
            THROW_EXCEPTION(ACCEPTABLE_POINT_REACHED,
                            "Restoration phase called at acceptable point.");
          }

          if (!IsValid(resto_phase_)) {
            THROW_EXCEPTION(IpoptException, "No Restoration Phase given to this Filter Line Search Object!");
          }
          // ToDo make the 1e-2 below a parameter?
          if (IpCq().curr_constraint_violation()<=
              1e-2*IpData().tol()) {
            bool found_acceptable = RestoreAcceptablePoint();
            if (found_acceptable) {
              THROW_EXCEPTION(ACCEPTABLE_POINT_REACHED,
                              "Restoration phase called at almost feasible point, but acceptable point could be restored.\n");
            }
            else {
              // ToDo does that happen too often?
              THROW_EXCEPTION(RESTORATION_FAILED,
                              "Restoration phase called, but point is almost feasible.");
            }
          }

          // Set the info fields for the first output line in the
          // restoration phase which reflects why the restoration phase
          // was called
          IpData().Set_info_alpha_primal(alpha_primal);
          IpData().Set_info_alpha_dual(0.);
          IpData().Set_info_alpha_primal_char('R');
          IpData().Set_info_ls_count(n_steps+1);

          accept = resto_phase_->PerformRestoration();
          if (!accept) {
            bool found_acceptable = RestoreAcceptablePoint();
            if (found_acceptable) {
              THROW_EXCEPTION(ACCEPTABLE_POINT_REACHED,
                              "Restoration phase failed, but acceptable point could be restore.\n");
            }
            else {
              THROW_EXCEPTION(RESTORATION_FAILED,
                              "Failed restoration phase!!!");
            }
          }
          count_successive_shortened_steps_ = 0;
          if (expect_infeasible_problem_) {
            expect_infeasible_problem_ = false;
          }
          in_soft_resto_phase_ = false;
          watchdog_shortened_iter_ = 0;
        }
      }
    }
    else if (!in_soft_resto_phase_ || tiny_step) {
      // we didn't do the restoration phase and are now updating the
      // dual variables of the trial point
      Number alpha_dual_max =
        IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                      *actual_delta->z_L(), *actual_delta->z_U(),
                                      *actual_delta->v_L(), *actual_delta->v_U());

      PerformDualStep(alpha_primal, alpha_dual_max, actual_delta);

      if (n_steps==0) {
        // accepted this if a full step was
        // taken
        count_successive_shortened_steps_ = 0;
        watchdog_shortened_iter_ = 0;
      }
      else {
        count_successive_shortened_steps_++;
        watchdog_shortened_iter_++;
      }

      if (expect_infeasible_problem_ &&
          IpCq().curr_constraint_violation() <= expect_infeasible_problem_ctol_) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Constraint violation is with %e less than expect_infeasible_problem_ctol.\nDisable expect_infeasible_problem_heuristic.\n", IpCq().curr_constraint_violation());
        expect_infeasible_problem_ = false;
      }
    }
  }

  bool FilterLineSearch::DoBacktrackingLineSearch(bool skip_first_trial_point,
      Number& alpha_primal,
      bool& corr_taken,
      bool& soc_taken,
      Index& n_steps,
      bool& evaluation_error,
      SmartPtr<IteratesVector>& actual_delta)
  {
    evaluation_error = false;
    bool accept = false;

    DBG_START_METH("FilterLineSearch::DoBacktrackingLineSearch",
                   dbg_verbosity);

    // Compute primal fraction-to-the-boundary value
    Number alpha_primal_max =
      IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                      *actual_delta->x(),
                                      *actual_delta->s());

    // Compute smallest step size allowed
    Number alpha_min = alpha_primal_max;
    if (!in_watchdog_) {
      alpha_min = CalculateAlphaMin();
    }
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "minimal step size ALPHA_MIN = %E\n", alpha_min);

    // Start line search from maximal step size
    alpha_primal = alpha_primal_max;

    // Step size used in ftype and armijo tests
    Number alpha_primal_test = alpha_primal;
    if (in_watchdog_) {
      alpha_primal_test = watchdog_alpha_primal_test_;
    }

    if (skip_first_trial_point) {
      alpha_primal *= alpha_red_factor_;
    }

    filter_.Print(Jnlst());

    if (corrector_type_!=NO_CORRECTOR && !skip_first_trial_point &&
        (!skip_corr_if_neg_curv_ || IpData().info_regu_x()==0.) &&
        (!skip_corr_in_monotone_mode_ || IpData().FreeMuMode()) ) {
      // Before we do the actual backtracking line search for the
      // regular primal-dual search direction, let's see if a step
      // including a higher-order correctior is already acceptable
      accept = TryCorrector(alpha_primal_test,
                            alpha_primal,
                            actual_delta);
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
          IpData().SetTrialPrimalVariablesFromStep(alpha_primal, *actual_delta->x(), *actual_delta->s());

          if (magic_steps_) {
            PerformMagicStep();
          }

          // If it is acceptable, stop the search
          alpha_primal_test = alpha_primal;
          accept = CheckAcceptabilityOfTrialPoint(alpha_primal_test);
        }
        catch(IpoptNLP::Eval_Error& e) {
          e.ReportException(Jnlst());
          Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                         "Warning: Cutting back alpha due to evaluation error\n");
          accept = false;
          evaluation_error = true;
        }

        if (accept) {
          break;
        }

        if (in_watchdog_) {
          break;
        }

        // Decide if we want to go to the restoration phase in a
        // short cut to check if the problem is infeasible
        if (expect_infeasible_problem_) {
          if (count_successive_shortened_steps_>=5) {
            break;
          }
        }

        // try second order correction step if the function could
        // be evaluated
        // DoTo: check if we want to do SOC when watchdog is active
        if (!evaluation_error) {
          Number theta_curr = IpCq().curr_constraint_violation();
          Number theta_trial = IpCq().trial_constraint_violation();
          if (alpha_primal==alpha_primal_max &&       // i.e. first trial point
              theta_curr<=theta_trial && max_soc_>0) {
            // Try second order correction
            accept = TrySecondOrderCorrection(alpha_primal_test,
                                              alpha_primal,
                                              actual_delta);
          }
          if (accept) {
            soc_taken = true;
            break;
          }
        }

        // Point is not yet acceptable, try a shorter one
        alpha_primal *= alpha_red_factor_;
        n_steps++;
      }
    } /* if (!accept) */

    char info_alpha_primal_char;
    if (!accept && in_watchdog_) {
      info_alpha_primal_char = 'w';
    }
    else {
      // Augment the filter if required
      if (!IsFtype(alpha_primal_test) ||
          !ArmijoHolds(alpha_primal_test)) {
        AugmentFilter();
        info_alpha_primal_char = 'h';
      }
      else {
        info_alpha_primal_char = 'f';
      }
    }
    if (soc_taken) {
      info_alpha_primal_char = toupper(info_alpha_primal_char);
    }
    IpData().Set_info_alpha_primal_char(info_alpha_primal_char);
    IpData().Set_info_ls_count(n_steps+1);
    if (corr_taken) {
      IpData().Append_info_string("C");
    }

    return accept;
  }

  bool FilterLineSearch::IsFtype(Number alpha_primal_test)
  {
    DBG_START_METH("FilterLineSearch::IsFtype",
                   dbg_verbosity);
    DBG_ASSERT(reference_theta_>0. || reference_gradBarrTDelta_ < 0.0);
    return (reference_gradBarrTDelta_ < 0.0 &&
            alpha_primal_test*pow(-reference_gradBarrTDelta_,s_phi_) >
            delta_*pow(reference_theta_,s_theta_));
  }

  void FilterLineSearch::AugmentFilter()
  {
    DBG_START_METH("FilterLineSearch::AugmentFilter",
                   dbg_verbosity);

    Number phi_add = reference_barr_ - gamma_phi_*reference_theta_;
    Number theta_add = (1.-gamma_theta_)*reference_theta_;

    filter_.AddEntry(phi_add, theta_add, IpData().iter_count());
  }

  bool
  FilterLineSearch::CheckAcceptabilityOfTrialPoint(Number alpha_primal_test)
  {
    DBG_START_METH("FilterLineSearch::CheckAcceptabilityOfTrialPoint",
                   dbg_verbosity);

    if (accept_every_trial_step_) {

      // We call the evaluation at the trial point here, so that an
      // exception will the thrown if there are problem during the
      // evaluation of the functions (in that case, we want to further
      // reduce the step size
      /* Number trial_barr = */ IpCq().trial_barrier_obj();
      /* Number trial_theta = */
      IpCq().trial_constraint_violation();
      return true;
    }

    bool accept;

    // First compute the barrier function and constraint violation at the
    // current iterate and the trial point

    Number trial_theta = IpCq().trial_constraint_violation();
    // Check if constraint violation is becoming too large
    if (theta_max_ < 0.0) {
      theta_max_ = theta_max_fact_*Max(1.0, reference_theta_);
    }
    if (theta_min_ < 0.0) {
      theta_min_ = theta_min_fact_*Max(1.0, reference_theta_);
    }

    if (theta_max_>0 && trial_theta>theta_max_) {
      return false;
    }

    Number trial_barr = IpCq().trial_barrier_obj();
    DBG_ASSERT(FiniteNumber(trial_barr));

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Checking acceptability for trial step size alpha_primal_test=%13.6e:\n", alpha_primal_test);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of barrier function     = %23.16e  (reference %23.16e):\n", trial_barr, reference_barr_);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of constraint violation = %23.16e  (reference %23.16e):\n", trial_theta, reference_theta_);

    // Check if point is acceptable w.r.t current iterate
    if (IsFtype(alpha_primal_test) && reference_theta_ <= theta_min_) {
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
    /*
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "ArmijoHolds test with trial_barr = %25.16e reference_barr = %25.16e\n        alpha_primal_test = %25.16e reference_gradBarrTDelta = %25.16e\n", IpCq().trial_barrier_obj(), reference_barr_,alpha_primal_test,reference_gradBarrTDelta_);
    */
    return Compare_le(IpCq().trial_barrier_obj()-reference_barr_,
                      eta_phi_*alpha_primal_test*reference_gradBarrTDelta_,
                      reference_barr_);
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

  bool FilterLineSearch::IsAcceptableToCurrentIterate(Number trial_barr,
      Number trial_theta,
      bool called_from_restoration /*=false*/) const
  {
    DBG_START_METH("FilterLineSearch::IsAcceptableToCurrentIterate",
                   dbg_verbosity);

    // Check if the barrier objective function is increasing to
    // rapidly (according to option obj_max_inc)
    if (!called_from_restoration && trial_barr > reference_barr_) {
      Number basval = 1.;
      if (fabs(reference_barr_)>10.) {
        basval = log10(fabs(reference_barr_));
      }
      if (log10(trial_barr-reference_barr_)>obj_max_inc_+basval) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Rejecting trial point because barrier objective function increasing too rapidly (from %27.15e to %27.15e)\n",reference_barr_,trial_barr);
        return false;
      }
    }

    DBG_PRINT((1,"trial_barr  = %e reference_barr  = %e\n", trial_barr, reference_barr_));
    DBG_PRINT((1,"trial_theta = %e reference_theta = %e\n", trial_theta, reference_theta_));
    return (Compare_le(trial_theta, (1.-gamma_theta_)*reference_theta_, reference_theta_)
            || Compare_le(trial_barr-reference_barr_, -gamma_phi_*reference_theta_, reference_barr_));
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

  bool FilterLineSearch::TrySoftRestoStep(SmartPtr<IteratesVector>& actual_delta,
                                          bool &satisfies_original_filter)
  {
    DBG_START_FUN("FilterLineSearch::TrySoftRestoStep", dbg_verbosity);

    satisfies_original_filter = false;

    // ToDo: Need to decide if we want to try a corrector step first

    // Compute the maximal step sizes (we use identical step sizes for
    // primal and dual variables
    Number alpha_primal_max =
      IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                      *actual_delta->x(),
                                      *actual_delta->s());
    Number alpha_dual_max =
      IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                    *actual_delta->z_L(),
                                    *actual_delta->z_U(),
                                    *actual_delta->v_L(),
                                    *actual_delta->v_U());
    Number alpha_max =  Min(alpha_primal_max, alpha_dual_max);

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Trying soft restoration phase step with step length %13.6e\n",
                   alpha_max);

    // Set the trial point
    IpData().SetTrialPrimalVariablesFromStep(alpha_max, *actual_delta->x(), *actual_delta->s());
    PerformDualStep(alpha_max, alpha_max,
                    actual_delta);

    // Check if that point is acceptable with respect to the current
    // original filter

    Number trial_barr;
    Number trial_theta;
    try {
      trial_barr = IpCq().trial_barrier_obj();
      trial_theta = IpCq().trial_constraint_violation();
    }
    catch(IpoptNLP::Eval_Error& e) {
      e.ReportException(Jnlst());
      Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                     "Warning: Evaluation error during soft restoration phase step.\n");
      return false;
    }
    if (theta_max_<=0 || trial_theta<=theta_max_) {
      if (IsAcceptableToCurrentIterate(trial_barr, trial_theta)) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "  Trial step acceptable with respect to original filter.\n");
        satisfies_original_filter = true;
        return true;
      }
    }

    // Evaluate the optimality error at the new point
    Number mu = .0;
    if (!IpData().FreeMuMode()) {
      mu = IpData().curr_mu();
    }
    Number trial_pderror;
    Number curr_pderror;
    try {
      trial_pderror = IpCq().trial_primal_dual_system_error(mu);
      curr_pderror = IpCq().curr_primal_dual_system_error(mu);
    }
    catch(IpoptNLP::Eval_Error& e) {
      e.ReportException(Jnlst());
      Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                     "Warning: Evaluation error during soft restoration phase step.\n");
      return false;
    }

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Primal-dual error at current point:  %23.16e\n", curr_pderror);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Primal-dual error at trial point  :  %23.16e\n", trial_pderror);
    // Check if there is sufficient reduction in the optimality error
    if (trial_pderror <= soft_resto_pderror_reduction_factor_*curr_pderror) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "  Trial step accepted.\n");
      return true;
    }

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Trial step rejected.\n");
    return false;
  }

  void FilterLineSearch::StartWatchDog()
  {
    DBG_START_FUN("FilterLineSearch::StartWatchDog", dbg_verbosity);

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Starting Watch Dog\n");

    in_watchdog_ = true;
    watchdog_iterate_ = IpData().curr();
    watchdog_delta_ = IpData().delta();
    watchdog_trial_iter_ = 0;
    watchdog_alpha_primal_test_ =
      IpCq().curr_primal_frac_to_the_bound(IpData().curr_tau());
    watchdog_theta_ = IpCq().curr_constraint_violation();
    watchdog_barr_ = IpCq().curr_barrier_obj();
    watchdog_gradBarrTDelta_ = IpCq().curr_gradBarrTDelta();
  }

  void FilterLineSearch::StopWatchDog(SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_FUN("FilterLineSearch::StopWatchDog", dbg_verbosity);

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Stopping Watch Dog\n");

    IpData().Append_info_string("w");

    in_watchdog_ = false;

    // Reset all fields in IpData to reference point
    SmartPtr<IteratesVector> old_trial = watchdog_iterate_->MakeNewContainer();
    IpData().set_trial(old_trial);
    IpData().AcceptTrialPoint();
    actual_delta = watchdog_delta_->MakeNewContainer();
    IpData().SetHaveAffineDeltas(false);

    // reset the stored watchdog iterates
    watchdog_iterate_ = NULL;
    watchdog_delta_ = NULL;

    watchdog_shortened_iter_ = 0;

    reference_theta_ = watchdog_theta_;
    reference_barr_ = watchdog_barr_;
    reference_gradBarrTDelta_ = watchdog_gradBarrTDelta_;
  }

  void FilterLineSearch::Reset()
  {
    DBG_START_FUN("FilterLineSearch::Reset", dbg_verbosity);
    in_soft_resto_phase_ = false;

    // Inactivate the watchdog and release all stored data
    in_watchdog_ = false;
    watchdog_iterate_ = NULL;
    watchdog_delta_ = NULL;
    watchdog_shortened_iter_ = 0;

    filter_.Clear();
  }

  void FilterLineSearch::PerformDualStep(Number alpha_primal,
                                         Number alpha_dual,
                                         SmartPtr<IteratesVector>& delta)
  {
    DBG_START_FUN("FilterLineSearch::PerformDualStep", dbg_verbosity);

    // set the bound multipliers from the step
    IpData().SetTrialBoundMultipliersFromStep(alpha_dual, *delta->z_L(), *delta->z_U(), *delta->v_L(), *delta->v_U());

    Number alpha_y;
    switch (alpha_for_y_) {
      case PRIMAL_ALPHA_FOR_Y:
      alpha_y = alpha_primal;
      break;
      case DUAL_ALPHA_FOR_Y:
      alpha_y = alpha_dual;
      break;
      case MIN_ALPHA_FOR_Y:
      alpha_y = Min(alpha_dual, alpha_primal);
      break;
      case MAX_ALPHA_FOR_Y:
      alpha_y = Max(alpha_dual, alpha_primal);
      break;
      case FULL_STEP_FOR_Y:
      alpha_y = 1;
      break;
      case MIN_DUAL_INFEAS_ALPHA_FOR_Y:
      case SAFE_MIN_DUAL_INFEAS_ALPHA_FOR_Y:
      // Here we compute the step size for y so that the dual
      // infeasibility is minimized along delta_y

      // compute the dual infeasibility at new point with old y
      SmartPtr<IteratesVector> temp_trial
      = IpData().trial()->MakeNewContainer();
      temp_trial->Set_y_c(*IpData().curr()->y_c());
      temp_trial->Set_y_d(*IpData().curr()->y_d());
      IpData().set_trial(temp_trial);
      SmartPtr<const Vector> dual_inf_x = IpCq().trial_grad_lag_x();
      SmartPtr<const Vector> dual_inf_s = IpCq().trial_grad_lag_s();

      SmartPtr<Vector> new_jac_times_delta_y =
        IpData().curr()->x()->MakeNew();
      new_jac_times_delta_y->AddTwoVectors(1., *IpCq().trial_jac_cT_times_vec(*delta->y_c()),
                                           1., *IpCq().trial_jac_dT_times_vec(*delta->y_d()),
                                           0.);

      Number a = pow(new_jac_times_delta_y->Nrm2(), 2.) +
                 pow(delta->y_d()->Nrm2(), 2.);
      Number b = dual_inf_x->Dot(*new_jac_times_delta_y) -
                 dual_inf_s->Dot(*delta->y_d());

      Number alpha = - b/a;

      if (alpha_for_y_==SAFE_MIN_DUAL_INFEAS_ALPHA_FOR_Y) {
        alpha_y = Min(Max(alpha_primal, alpha_dual), Max(alpha, Min(alpha_primal, alpha_dual)));
      }
      else {
        alpha_y = Min(1., Max(0., alpha));
      }
      break;
    }

    // Set the eq multipliers from the step now that alpha_y
    // has been calculated.
    DBG_PRINT((1, "alpha_y = %e\n", alpha_y));
    DBG_PRINT_VECTOR(2, "delta_y_c", *delta->y_c());
    DBG_PRINT_VECTOR(2, "delta_y_d", *delta->y_d());
    IpData().SetTrialEqMultipliersFromStep(alpha_y, *delta->y_c(), *delta->y_d());

    // Set some information for iteration summary output
    IpData().Set_info_alpha_primal(alpha_primal);
    IpData().Set_info_alpha_dual(alpha_dual);
  }

  bool
  FilterLineSearch::TrySecondOrderCorrection(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_METH("FilterLineSearch::TrySecondOrderCorrection",
                   dbg_verbosity);

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
      pd_solver_->Solve(-1.0, 0.0, *rhs, *delta_soc, true);

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
        e.ReportException(Jnlst());
        Jnlst().Printf(J_WARNING, J_MAIN, "Warning: SOC step rejected due to evaluation error\n");
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
  FilterLineSearch::TryCorrector(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_METH("FilterLineSearch::TryCorrector",
                   dbg_verbosity);

    Index n_bounds = IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim()
                     + IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim();
    if (n_bounds==0) {
      // Nothing to be done
      return false;
    }

    bool accept = false;

    // Compute the corrector step based on corrector_type parameter
    // create a new iterates vector and allocate space for all the entries
    SmartPtr<IteratesVector> delta_corr = actual_delta->MakeNewIteratesVector(true);

    switch (corrector_type_) {
      case AFFINE_CORRECTOR : {
        // 1: Standard MPC corrector

        if (!IpData().HaveAffineDeltas()) {
          Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                         "Solving the Primal Dual System for the affine step\n");
          // First get the right hand side
          SmartPtr<IteratesVector> rhs_aff = delta_corr->MakeNewContainer();

          rhs_aff->Set_x(*IpCq().curr_grad_lag_x());
          rhs_aff->Set_s(*IpCq().curr_grad_lag_s());
          rhs_aff->Set_y_c(*IpCq().curr_c());
          rhs_aff->Set_y_d(*IpCq().curr_d_minus_s());
          rhs_aff->Set_z_L(*IpCq().curr_compl_x_L());
          rhs_aff->Set_z_U(*IpCq().curr_compl_x_U());
          rhs_aff->Set_v_L(*IpCq().curr_compl_s_L());
          rhs_aff->Set_v_U(*IpCq().curr_compl_s_U());

          // create a new iterates vector (with allocated space)
          // for the affine scaling step
          SmartPtr<IteratesVector> step_aff = delta_corr->MakeNewIteratesVector(true);

          // Now solve the primal-dual system to get the step
          pd_solver_->Solve(-1.0, 0.0, *rhs_aff, *step_aff, false);

          DBG_PRINT_VECTOR(2, "step_aff", *step_aff);

          IpData().set_delta_aff(step_aff);
          IpData().SetHaveAffineDeltas(true);
        }

        DBG_ASSERT(IpData().HaveAffineDeltas());

        const SmartPtr<const IteratesVector> delta_aff = IpData().delta_aff();

        delta_corr->Copy(*actual_delta);

        // create a rhs vector and allocate entries
        SmartPtr<IteratesVector> rhs = actual_delta->MakeNewIteratesVector(true);

        rhs->x_NonConst()->Set(0.);
        rhs->s_NonConst()->Set(0.);
        rhs->y_c_NonConst()->Set(0.);
        rhs->y_d_NonConst()->Set(0.);
        IpNLP().Px_L()->TransMultVector(-1., *delta_aff->x(), 0., *rhs->z_L_NonConst());
        rhs->z_L_NonConst()->ElementWiseMultiply(*delta_aff->z_L());
        IpNLP().Px_U()->TransMultVector(1., *delta_aff->x(), 0., *rhs->z_U_NonConst());
        rhs->z_U_NonConst()->ElementWiseMultiply(*delta_aff->z_U());
        IpNLP().Pd_L()->TransMultVector(-1., *delta_aff->s(), 0., *rhs->v_L_NonConst());
        rhs->v_L_NonConst()->ElementWiseMultiply(*delta_aff->v_L());
        IpNLP().Pd_U()->TransMultVector(1., *delta_aff->s(), 0., *rhs->v_U_NonConst());
        rhs->v_U_NonConst()->ElementWiseMultiply(*delta_aff->v_U());

        pd_solver_->Solve(1.0, 1.0, *rhs, *delta_corr, true);

        DBG_PRINT_VECTOR(2, "delta_corr", *delta_corr);
      }
      break;
      case PRIMAL_DUAL_CORRECTOR : {
        // 2: Second order correction for primal-dual step to
        // primal-dual mu

        delta_corr->Copy(*actual_delta);

        // allocate space for the rhs
        SmartPtr<IteratesVector> rhs = actual_delta->MakeNewIteratesVector(true);

        rhs->x_NonConst()->Set(0.);
        rhs->s_NonConst()->Set(0.);
        rhs->y_c_NonConst()->Set(0.);
        rhs->y_d_NonConst()->Set(0.);

        Number mu = IpData().curr_mu();
        SmartPtr<Vector> tmp;

        rhs->z_L_NonConst()->Copy(*IpCq().curr_slack_x_L());
        IpNLP().Px_L()->TransMultVector(-1., *actual_delta->x(),
                                        -1., *rhs->z_L_NonConst());
        tmp = actual_delta->z_L()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->z_L(), 1., *actual_delta->z_L(), 0.);
        rhs->z_L_NonConst()->ElementWiseMultiply(*tmp);
        rhs->z_L_NonConst()->AddScalar(mu);

        rhs->z_U_NonConst()->Copy(*IpCq().curr_slack_x_U());
        IpNLP().Px_U()->TransMultVector(1., *actual_delta->x(),
                                        -1., *rhs->z_U_NonConst());
        tmp = actual_delta->z_U()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->z_U(), 1., *actual_delta->z_U(), 0.);
        rhs->z_U_NonConst()->ElementWiseMultiply(*tmp);
        rhs->z_U_NonConst()->AddScalar(mu);

        rhs->v_L_NonConst()->Copy(*IpCq().curr_slack_s_L());
        IpNLP().Pd_L()->TransMultVector(-1., *actual_delta->s(),
                                        -1., *rhs->v_L_NonConst());
        tmp = actual_delta->v_L()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->v_L(), 1., *actual_delta->v_L(), 0.);
        rhs->v_L_NonConst()->ElementWiseMultiply(*tmp);
        rhs->v_L_NonConst()->AddScalar(mu);

        rhs->v_U_NonConst()->Copy(*IpCq().curr_slack_s_U());
        IpNLP().Pd_U()->TransMultVector(1., *actual_delta->s(),
                                        -1., *rhs->v_U_NonConst());
        tmp = actual_delta->v_U()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->v_U(), 1., *actual_delta->v_U(), 0.);
        rhs->v_U_NonConst()->ElementWiseMultiply(*tmp);
        rhs->v_U_NonConst()->AddScalar(mu);

        DBG_PRINT_VECTOR(2, "rhs", *rhs);

        pd_solver_->Solve(1.0, 1.0, *rhs, *delta_corr, true);

        DBG_PRINT_VECTOR(2, "delta_corr", *delta_corr);
      }
      break;
      default:
      DBG_ASSERT("Unknown corrector_type value.");
    }

    // Compute step size
    Number alpha_primal_corr =
      IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                      *delta_corr->x(),
                                      *delta_corr->s());
    // Set the primal trial point
    IpData().SetTrialPrimalVariablesFromStep(alpha_primal_corr, *delta_corr->x(), *delta_corr->s());

    // Check if we want to not even try the filter criterion
    Number alpha_dual_max =
      IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                    *delta_corr->z_L(), *delta_corr->z_U(),
                                    *delta_corr->v_L(), *delta_corr->v_U());

    IpData().SetTrialBoundMultipliersFromStep(alpha_dual_max, *delta_corr->z_L(), *delta_corr->z_U(), *delta_corr->v_L(), *delta_corr->v_U());

    Number trial_avrg_compl = IpCq().trial_avrg_compl();
    Number curr_avrg_compl = IpCq().curr_avrg_compl();
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "avrg_compl(curr) = %e, avrg_compl(trial) = %e\n",
                   curr_avrg_compl, trial_avrg_compl);
    if (corrector_type_==AFFINE_CORRECTOR &&
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
      Jnlst().Printf(J_WARNING, J_MAIN,
                     "Warning: Corrector step rejected due to evaluation error\n");
      accept = false;
    }

    if (accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Corrector step accepted with alpha_primal = %e\n",
                     alpha_primal_corr);
      // Accept all SOC quantities
      alpha_primal = alpha_primal_corr;
      actual_delta = delta_corr;

      if (Jnlst().ProduceOutput(J_MOREVECTOR, J_MAIN)) {
        Jnlst().Printf(J_MOREVECTOR, J_MAIN,
                       "*** Accepted corrector for Iteration: %d\n",
                       IpData().iter_count());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_corr", *delta_corr);
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
    DBG_PRINT_VECTOR(2, "s in", *IpData().trial()->s());
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

    SmartPtr<Vector> delta_s_magic = IpData().trial()->s()->MakeNew();
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
    tmp->Copy(*IpData().trial()->s());
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
      if (delta_s_magic_max > 10*mach_eps*IpData().trial()->s()->Amax()) {
        IpData().Append_info_string("M");
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Magic step with max-norm %.6e taken.\n", delta_s_magic->Amax());
        Jnlst().PrintVector(J_MOREVECTOR, J_LINE_SEARCH,
                            "delta_s_magic", *delta_s_magic);
      }

      // now finally compute the new overall slacks
      delta_s_magic->Axpy(1., *IpData().trial()->s());
      SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
      trial->Set_s(*delta_s_magic);
      IpData().set_trial(trial);
    }

    DBG_PRINT((1,"Outgoing barr = %e and constrviol %e\n", IpCq().trial_barrier_obj(), IpCq().trial_constraint_violation()));
    DBG_PRINT_VECTOR(2, "s out", *IpData().trial()->s());
    DBG_PRINT_VECTOR(2, "d minus s out", *IpCq().trial_d_minus_s());
    DBG_PRINT_VECTOR(2, "slack_s_L out", *IpCq().trial_slack_s_L());
    DBG_PRINT_VECTOR(2, "slack_s_U out", *IpCq().trial_slack_s_U());
  }

  bool
  FilterLineSearch::DetectTinyStep()
  {
    DBG_START_METH("FilterLineSearch::DetectTinyStep",
                   dbg_verbosity);

    Number max_step_x;
    Number max_step_s;

    if (tiny_step_tol_==0.)
      return false;

    // ToDo try to find more efficient implementation

    SmartPtr<Vector> tmp = IpData().curr()->x()->MakeNew();
    tmp->Copy(*IpData().curr()->x());
    tmp->ElementWiseAbs();
    tmp->AddScalar(1.);

    SmartPtr<Vector> tmp2 = IpData().curr()->x()->MakeNew();
    tmp2->Copy(*IpData().delta()->x());
    tmp2->ElementWiseDivide(*tmp);
    max_step_x = tmp2->Amax();
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "Relative step size for delta_x = %e\n",
                   max_step_x);
    if (max_step_x > tiny_step_tol_)
      return false;

    tmp = IpData().curr()->s()->MakeNew();
    tmp->Copy(*IpData().curr()->s());
    tmp->ElementWiseAbs();
    tmp->AddScalar(1.);

    tmp2 = IpData().curr()->s()->MakeNew();
    tmp2->Copy(*IpData().delta()->s());
    tmp2->ElementWiseDivide(*tmp);
    max_step_s = tmp2->Amax();
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "Relative step size for delta_s = %e\n",
                   max_step_s);
    if (max_step_s > tiny_step_tol_)
      return false;

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Tiny step of relative size %e detected.\n",
                   Max(max_step_x, max_step_s));

    return true;
  }

  bool FilterLineSearch::CurrentIsAcceptable()
  {
    return (IsValid(conv_check_) &&
            conv_check_->CurrentIsAcceptable());
  }

  void FilterLineSearch::StoreAcceptablePoint()
  {
    DBG_START_METH("FilterLineSearch::StoreAcceptablePoint",
                   dbg_verbosity);

    acceptable_iterate_ = IpData().curr();
  }

  bool FilterLineSearch::RestoreAcceptablePoint()
  {
    DBG_START_METH("FilterLineSearch::RestoreAcceptablePoint",
                   dbg_verbosity);

    if (!IsValid(acceptable_iterate_)) {
      return false;
    }

    SmartPtr<IteratesVector> prev_iterate = acceptable_iterate_->MakeNewContainer();
    IpData().set_trial(prev_iterate);
    IpData().AcceptTrialPoint();

    return true;
  }


} // namespace Ipopt
