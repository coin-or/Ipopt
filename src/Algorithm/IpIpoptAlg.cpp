// Copyright (C) 2004, 2012 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpIpoptAlg.hpp"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpBacktrackingLineSearch.hpp"

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif

#if IPOPT_CHECKLEVEL > 1 && defined(IPOPT_HAS_FEENABLEEXCEPT)
#include <cfenv>
#endif

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

IpoptAlgorithm::IpoptAlgorithm(
   const SmartPtr<SearchDirectionCalculator>& search_dir_calculator,
   const SmartPtr<LineSearch>&                line_search,
   const SmartPtr<MuUpdate>&                  mu_update,
   const SmartPtr<ConvergenceCheck>&          conv_check,
   const SmartPtr<IterateInitializer>&        iterate_initializer,
   const SmartPtr<IterationOutput>&           iter_output,
   const SmartPtr<HessianUpdater>&            hessian_updater,
   const SmartPtr<EqMultiplierCalculator>&    eq_multiplier_calculator, /* = NULL*/
   const std::string&                         linear_solver_name /* = "" */
)
   : search_dir_calculator_(search_dir_calculator),
     line_search_(line_search),
     mu_update_(mu_update),
     conv_check_(conv_check),
     iterate_initializer_(iterate_initializer),
     iter_output_(iter_output),
     hessian_updater_(hessian_updater),
     eq_multiplier_calculator_(eq_multiplier_calculator),
     linear_solver_name_(linear_solver_name)
{
   DBG_START_METH("IpoptAlgorithm::IpoptAlgorithm",
                  dbg_verbosity);
   DBG_ASSERT(IsValid(search_dir_calculator_));
   DBG_ASSERT(IsValid(line_search_));
   DBG_ASSERT(IsValid(mu_update_));
   DBG_ASSERT(IsValid(conv_check_));
   DBG_ASSERT(IsValid(iterate_initializer_));
   DBG_ASSERT(IsValid(iter_output_));
   DBG_ASSERT(IsValid(hessian_updater_));
}

IpoptAlgorithm::~IpoptAlgorithm()
{
   DBG_START_METH("IpoptAlgorithm::~IpoptAlgorithm()",
                  dbg_verbosity);
}

void IpoptAlgorithm::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->SetRegisteringCategory("Line Search");
   roptions->AddLowerBoundedNumberOption(
      "kappa_sigma",
      "Factor limiting the deviation of dual variables from primal estimates.",
      0., true,
      1e10,
      "If the dual variables deviate from their primal estimates, a correction is performed. "
      "See Eqn. (16) in the implementation paper. "
      "Setting the value to less than 1 disables the correction.",
      true);
   roptions->AddStringOption2(
      "recalc_y",
      "Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates.",
      "no",
      "no", "use the Newton step to update the multipliers",
      "yes", "use least-square multiplier estimates",
      "This asks the algorithm to recompute the multipliers, whenever the current infeasibility is less than recalc_y_feas_tol. "
      "Choosing yes might be helpful in the quasi-Newton option. "
      "However, each recalculation requires an extra factorization of the linear system. "
      "If a limited memory quasi-Newton option is chosen, this is used by default.");
   roptions->AddLowerBoundedNumberOption(
      "recalc_y_feas_tol",
      "Feasibility threshold for recomputation of multipliers.",
      0., true,
      1e-6,
      "If recalc_y is chosen and the current infeasibility is less than this value, then the multipliers are recomputed.");
   roptions->SetRegisteringCategory("Step Calculation");
   roptions->AddBoolOption(
      "mehrotra_algorithm",
      "Indicates whether to do Mehrotra's predictor-corrector algorithm.",
      false,
      "If enabled, line search is disabled and the (unglobalized) adaptive mu strategy is chosen "
      "with the \"probing\" oracle, and \"corrector_type=affine\" is used without any safeguards; "
      "you should not set any of those options explicitly in addition. "
      "Also, unless otherwise specified, the values of \"bound_push\", \"bound_frac\", and "
      "\"bound_mult_init_val\" are set more aggressive, and sets \"alpha_for_y=bound_mult\". "
      "The Mehrotra's predictor-corrector algorithm works usually very well for LPs and convex QPs.");
   roptions->SetRegisteringCategory("Undocumented");
   roptions->AddBoolOption("sb",
                           "whether to skip printing Ipopt copyright banner",
                           false);
   roptions->SetRegisteringCategory("Miscellaneous");
   roptions->AddBoolOption(
      "timing_statistics",
      "Indicates whether to measure time spend in components of Ipopt and NLP evaluation",
      false,
      "The overall algorithm time is unaffected by this option.");
}

static bool copyright_message_printed = false;

bool IpoptAlgorithm::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   DBG_START_METH("IpoptAlgorithm::InitializeImpl",
                  dbg_verbosity);

   // disable detailed timing, if not required
   bool timing_statistics;
   options.GetBoolValue("timing_statistics", timing_statistics, "");
   if( !timing_statistics )
   {
      IpData().TimingStats().DisableTimes();
   }

   SmartPtr<const OptionsList> my_options;
   options.GetBoolValue("mehrotra_algorithm", mehrotra_algorithm_, prefix);
   if( mehrotra_algorithm_ )
   {
      // Verify a few options and set a few new ones.  But we better
      // make a copy of the incoming options.
      SmartPtr<OptionsList> new_options = new OptionsList(options);
      // Check required options are set correctly
      std::string string_option;
      if( new_options->GetStringValue("adaptive_mu_globalization", string_option, prefix) )
      {
         ASSERT_EXCEPTION(string_option == "never-monotone-mode", OPTION_INVALID,
                          "If mehrotra_algorithm=yes, adaptive_mu_globalization must be \"never-monotone-mode\".");
      }
      else
      {
         new_options->SetStringValue("adaptive_mu_globalization", "never-monotone-mode", false);
      }
      // The corrector step is already taken case of in
      // ComputeSearchDirection below
      if( new_options->GetStringValue("corrector_type", string_option, prefix) )
      {
         ASSERT_EXCEPTION(string_option == "none", OPTION_INVALID,
                          "If mehrotra_algorithm=yes, corrector_type must be \"none\".");
      }
      else
      {
         new_options->SetStringValue("corrector_type", "none", false);
      }
      if( new_options->GetStringValue("accept_every_trial_step", string_option, prefix) )
      {
         ASSERT_EXCEPTION(string_option == "yes", OPTION_INVALID,
                          "If mehrotra_algorithm=yes, accept_every_trial_step must be \"yes\".");
      }
      else
      {
         new_options->SetStringValue("accept_every_trial_step", "yes", false);
      }

      // Change some default options
      new_options->SetNumericValueIfUnset("bound_push", 10.);
      new_options->SetNumericValueIfUnset("bound_frac", 0.2);
      new_options->SetNumericValueIfUnset("bound_mult_init_val", 10.);
      new_options->SetNumericValueIfUnset("constr_mult_init_max", 0.);
      new_options->SetStringValueIfUnset("alpha_for_y", "bound_mult");
      new_options->SetStringValueIfUnset("least_square_init_primal", "yes");

      my_options = ConstPtr(new_options);
   }
   else
   {
      my_options = &options;
   }
   bool bval;
   options.GetBoolValue("sb", bval, prefix);
   if( bval )
   {
      copyright_message_printed = true;
   }

   // Read the IpoptAlgorithm options
   // Initialize the Data object
   bool retvalue = IpData().Initialize(Jnlst(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the IpIpoptData object failed to initialize.");

   // Initialize the CQ object
   retvalue = IpCq().Initialize(Jnlst(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the IpIpoptCalculatedQuantities object failed to initialize.");

   // Initialize the NLP object
   retvalue = IpNLP().Initialize(Jnlst(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the IpIpoptNLP object failed to initialize.");

   // Initialize all the strategies
   retvalue = iterate_initializer_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the iterate_initializer strategy failed to initialize.");

   retvalue = mu_update_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the mu_update strategy failed to initialize.");

   retvalue = search_dir_calculator_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the search_direction_calculator strategy failed to initialize.");
   retvalue = line_search_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the line_search strategy failed to initialize.");

   retvalue = conv_check_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the conv_check strategy failed to initialize.");

   retvalue = iter_output_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the iter_output strategy failed to initialize.");

   retvalue = hessian_updater_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
   ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION, "the hessian_updater strategy failed to initialize.");

   my_options->GetNumericValue("kappa_sigma", kappa_sigma_, prefix);
   if( !my_options->GetBoolValue("recalc_y", recalc_y_, prefix) )
   {
      Index enum_int;
      if( my_options->GetEnumValue("hessian_approximation", enum_int, prefix) )
      {
         HessianApproximationType hessian_approximation = HessianApproximationType(enum_int);
         if( hessian_approximation == LIMITED_MEMORY )
         {
            recalc_y_ = true;
         }
      }
   }
   if( recalc_y_ )
   {
      my_options->GetNumericValue("recalc_y_feas_tol", recalc_y_feas_tol_, prefix);
   }

   my_options->GetNumericValue("constr_viol_tol", constr_viol_tol_, prefix);

   if( prefix == "resto." )
   {
      skip_print_problem_stats_ = true;
   }
   else
   {
      skip_print_problem_stats_ = false;
   }

   return true;
}

// class to execute some tasks when object is destructed
// ends a TimedTask, disables floating-point exceptions
class EndTasks
{
private:
   TimedTask& task_;
   int excepts_;
public:
   EndTasks(
      TimedTask& task,
      int        excepts
   )
      : task_(task),
        excepts_(excepts)
   {
      DBG_ASSERT(task_.IsStarted());
   }

   ~EndTasks()
   {
      task_.End();
#if IPOPT_CHECKLEVEL > 1 && defined(IPOPT_HAS_FEENABLEEXCEPT)
      fedisableexcept(excepts_);
#else
      (void)excepts_;
#endif
   }
};

SolverReturn IpoptAlgorithm::Optimize(
   bool isResto /*= false */
)
{
   DBG_START_METH("IpoptAlgorithm::Optimize", dbg_verbosity);

   // Start measuring CPU time
   IpData().TimingStats().OverallAlgorithm().Start();

   int disable_excepts = 0;
#if IPOPT_CHECKLEVEL > 1 && defined(IPOPT_HAS_FEENABLEEXCEPT)
   int orig_excepts = feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
   if( orig_excepts >= 0 )
   {
      // disable exceptions FE_DIVBYZERO, FE_INVALID, FE_OVERFLOW later, if not part of originally set exceptions
      disable_excepts = ~orig_excepts & (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
   }
#endif

   // ensure timed task is ended when Optimize() is left
   // ensure additional floating point exceptions are no longer raised when Optimize() is left
   EndTasks endtasks(IpData().TimingStats().OverallAlgorithm(), disable_excepts);

   if( !copyright_message_printed )
   {
      print_copyright_message(Jnlst());
   }

   if( !isResto )
   {
      Jnlst().Printf(J_ITERSUMMARY, J_MAIN,
                     "This is Ipopt version " IPOPT_VERSION ", running with linear solver %s.\n\n", linear_solver_name_.c_str());
   }

   SolverReturn retval = UNASSIGNED;

   try
   {
      IpData().TimingStats().InitializeIterates().Start();
      // Initialize the iterates
      InitializeIterates();
      IpData().TimingStats().InitializeIterates().End();

      if( !skip_print_problem_stats_ )
      {
         IpData().TimingStats().PrintProblemStatistics().Start();
         PrintProblemStatistics();
         IpData().TimingStats().PrintProblemStatistics().End();
      }

      IpData().TimingStats().CheckConvergence().Start();
      ConvergenceCheck::ConvergenceStatus conv_status = conv_check_->CheckConvergence();
      IpData().TimingStats().CheckConvergence().End();

      // main loop
      while( conv_status == ConvergenceCheck::CONTINUE )
      {
         // Set the Hessian Matrix
         IpData().TimingStats().UpdateHessian().Start();
         UpdateHessian();
         IpData().TimingStats().UpdateHessian().End();

         // do all the output for this iteration
         IpData().TimingStats().OutputIteration().Start();
         OutputIteration();
         IpData().ResetInfo();
         IpData().TimingStats().OutputIteration().End();

         // initialize the flag that is set to true if the algorithm
         // has to continue with an emergency fallback mode.  For
         // example, when no search direction can be computed, continue
         // with the restoration phase
         bool emergency_mode = false;

         // update the barrier parameter
         IpData().TimingStats().UpdateBarrierParameter().Start();
         emergency_mode = !UpdateBarrierParameter();
         IpData().TimingStats().UpdateBarrierParameter().End();

         if( !emergency_mode )
         {
            // solve the primal-dual system to get the full step
            IpData().TimingStats().ComputeSearchDirection().Start();
            emergency_mode = !ComputeSearchDirection();
            IpData().TimingStats().ComputeSearchDirection().End();
         }

         // If we are in the emergency mode, ask the line search object
         // to go to the fallback options.  If that isn't possible,
         // issue error message
         if( emergency_mode )
         {
            if( line_search_->ActivateFallbackMechanism() )
            {
               Jnlst().Printf(J_WARNING, J_MAIN,
                              "WARNING: Problem in step computation; switching to emergency mode.\n");
            }
            else
            {
               Jnlst().Printf(J_ERROR, J_MAIN,
                              "ERROR: Problem in step computation, but emergency mode cannot be activated.\n");
               THROW_EXCEPTION(STEP_COMPUTATION_FAILED, "Step computation failed.");
            }
         }

         // Compute the new iterate
         IpData().TimingStats().ComputeAcceptableTrialPoint().Start();
         ComputeAcceptableTrialPoint();
         IpData().TimingStats().ComputeAcceptableTrialPoint().End();

         // Accept the new iterate
         IpData().TimingStats().AcceptTrialPoint().Start();
         AcceptTrialPoint();
         IpData().TimingStats().AcceptTrialPoint().End();

         IpData().Set_iter_count(IpData().iter_count() + 1);

         if( IpCq().IsSquareProblem() )
         {
            ComputeFeasibilityMultipliers();
         }

         IpData().TimingStats().CheckConvergence().Start();
         conv_status = conv_check_->CheckConvergence();
         IpData().TimingStats().CheckConvergence().End();
      }

      IpData().TimingStats().OutputIteration().Start();
      OutputIteration();
      IpData().TimingStats().OutputIteration().End();

      bool stop_watchdog = false;
      switch( conv_status )
      {
         case ConvergenceCheck::CONVERGED:
            retval = SUCCESS;
            break;
         case ConvergenceCheck::CONVERGED_TO_ACCEPTABLE_POINT:
            retval = STOP_AT_ACCEPTABLE_POINT;
            break;
         case ConvergenceCheck::MAXITER_EXCEEDED:
            retval = MAXITER_EXCEEDED;
            stop_watchdog = true;
            break;
         case ConvergenceCheck::CPUTIME_EXCEEDED:
            retval = CPUTIME_EXCEEDED;
            stop_watchdog = true;
            break;
         case ConvergenceCheck::WALLTIME_EXCEEDED:
            retval = WALLTIME_EXCEEDED;
            stop_watchdog = true;
            break;
         case ConvergenceCheck::DIVERGING:
            retval = DIVERGING_ITERATES;
            break;
         case ConvergenceCheck::USER_STOP:
            retval = USER_REQUESTED_STOP;
            break;
         default:
            retval = INTERNAL_ERROR;
            break;
      }

      // stop watchdog if interrupted due to limit to restore iterate from before watchdog (#289)
      if( stop_watchdog && dynamic_cast<BacktrackingLineSearch*>(GetRawPtr(line_search_)) != NULL )
      {
         static_cast<BacktrackingLineSearch*>(GetRawPtr(line_search_))->StopWatchDog();
      }
   }
   catch( TINY_STEP_DETECTED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().UpdateBarrierParameter().EndIfStarted();
      retval = STOP_AT_TINY_STEP;
   }
   catch( ACCEPTABLE_POINT_REACHED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = STOP_AT_ACCEPTABLE_POINT;
   }
   catch( LOCALLY_INFEASIBLE& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().CheckConvergence().EndIfStarted();
      retval = LOCAL_INFEASIBILITY;
   }
   catch( RESTORATION_CONVERGED_TO_FEASIBLE_POINT& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      if( IpCq().IsSquareProblem() )
      {
         // make sure the multipliers are computed properly
         ComputeFeasibilityMultipliersPostprocess();
         retval = FEASIBLE_POINT_FOUND;
      }
      else
      {
         retval = RESTORATION_FAILURE;
      }
   }
   catch( RESTORATION_FAILED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = RESTORATION_FAILURE;
   }
   catch( RESTORATION_MAXITER_EXCEEDED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = MAXITER_EXCEEDED;
   }
   catch( RESTORATION_CPUTIME_EXCEEDED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = CPUTIME_EXCEEDED;
   }
   catch( RESTORATION_WALLTIME_EXCEEDED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = WALLTIME_EXCEEDED;
   }
   catch( RESTORATION_USER_STOP& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = USER_REQUESTED_STOP;
   }
   catch( STEP_COMPUTATION_FAILED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = ERROR_IN_STEP_COMPUTATION;
   }
   catch( IpoptNLP::Eval_Error& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().InitializeIterates().EndIfStarted();
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = INVALID_NUMBER_DETECTED;
   }
   catch( FEASIBILITY_PROBLEM_SOLVED& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      if( IpCq().IsSquareProblem() )
      {
         // make the sure multipliers are computed properly
         ComputeFeasibilityMultipliersPostprocess();
      }
      retval = FEASIBLE_POINT_FOUND;
   }
   catch( TOO_FEW_DOF& exc )
   {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      retval = TOO_FEW_DEGREES_OF_FREEDOM;
   }
   catch( INTERNAL_ABORT& exc )
   {
      exc.ReportException(Jnlst());
      retval = INTERNAL_ERROR;
   }

   DBG_ASSERT(retval != UNASSIGNED && "Unknown return code in the algorithm");
   return retval;
}

void IpoptAlgorithm::UpdateHessian()
{
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n");
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "*** Update HessianMatrix for Iteration %" IPOPT_INDEX_FORMAT ":", IpData().iter_count());
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n\n");
   hessian_updater_->UpdateHessian();
}

bool IpoptAlgorithm::UpdateBarrierParameter()
{
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n");
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "*** Update Barrier Parameter for Iteration %" IPOPT_INDEX_FORMAT ":", IpData().iter_count());
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n\n");
   bool retval = mu_update_->UpdateBarrierParameter();

   if( retval )
   {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Barrier Parameter: %e\n", IpData().curr_mu());
   }
   else
   {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Barrier parameter could not be updated!\n");
   }

   return retval;
}

bool IpoptAlgorithm::ComputeSearchDirection()
{
   DBG_START_METH("IpoptAlgorithm::ComputeSearchDirection", dbg_verbosity);

   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n");
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "*** Solving the Primal Dual System for Iteration %" IPOPT_INDEX_FORMAT ":", IpData().iter_count());
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n\n");

   bool retval = search_dir_calculator_->ComputeSearchDirection();

   if( retval )
   {
      Jnlst().Printf(J_MOREVECTOR, J_MAIN,
                     "*** Step Calculated for Iteration: %" IPOPT_INDEX_FORMAT "\n", IpData().iter_count());
      IpData().delta()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "delta");
   }
   else
   {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "*** Step could not be computed in iteration %" IPOPT_INDEX_FORMAT "!\n", IpData().iter_count());
   }

   return retval;
}

void IpoptAlgorithm::ComputeAcceptableTrialPoint()
{
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n");
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "*** Finding Acceptable Trial Point for Iteration %" IPOPT_INDEX_FORMAT ":", IpData().iter_count());
   Jnlst().Printf(J_DETAILED, J_MAIN,
                  "\n**************************************************\n\n");
   line_search_->FindAcceptableTrialPoint();
}

void IpoptAlgorithm::OutputIteration()
{
   iter_output_->WriteOutput();
}

void IpoptAlgorithm::InitializeIterates()
{
   DBG_START_METH("IpoptAlgorithm::InitializeIterates", dbg_verbosity);

   bool retval = iterate_initializer_->SetInitialIterates();
   ASSERT_EXCEPTION(retval, FAILED_INITIALIZATION, "Error while obtaining initial iterates.");
}

void IpoptAlgorithm::AcceptTrialPoint()
{
   DBG_START_METH("IpoptAlgorithm::AcceptTrialPoint", dbg_verbosity);
   // If the line search didn't determine a new acceptable trial
   // point, do not accept a new iterate
   if( line_search_->CheckSkippedLineSearch() )
   {
      Jnlst().Printf(J_SUMMARY, J_MAIN,
                     "Line search didn't find acceptable trial point.\n");
      return;
   }

   // Adjust the bounds if necessary
   Index adjusted_slacks = IpCq().AdjustedTrialSlacks();
   DBG_PRINT((1, "adjusted_slacks = %" IPOPT_INDEX_FORMAT "\n", adjusted_slacks));
   if( adjusted_slacks > 0 )
   {
      IpCq().ResetAdjustedTrialSlacks();
      if( adjusted_slacks == 1 )
      {
         Jnlst().Printf(J_WARNING, J_MAIN,
                        "In iteration %" IPOPT_INDEX_FORMAT ", %" IPOPT_INDEX_FORMAT " Slack too small, adjusting variable bound\n", IpData().iter_count(), adjusted_slacks);
      }
      else
      {
         Jnlst().Printf(J_WARNING, J_MAIN,
                        "In iteration %" IPOPT_INDEX_FORMAT ", %" IPOPT_INDEX_FORMAT " Slacks too small, adjusting variable bounds\n", IpData().iter_count(), adjusted_slacks);
      }
      if( Jnlst().ProduceOutput(J_VECTOR, J_MAIN) )
      {
         IpNLP().x_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_x_L");
         IpNLP().x_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_x_U");
         IpNLP().d_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_d_L");
         IpNLP().d_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_d_U");
      }

      SmartPtr<Vector> new_x_l = IpNLP().x_L()->MakeNew();
      IpNLP().Px_L()->TransMultVector(1.0, *IpData().trial()->x(), 0.0, *new_x_l);
      new_x_l->Axpy(-1.0, *IpCq().trial_slack_x_L());

      SmartPtr<Vector> new_x_u = IpNLP().x_U()->MakeNew();
      IpNLP().Px_U()->TransMultVector(1.0, *IpData().trial()->x(), 0.0, *new_x_u);
      new_x_u->Axpy(1.0, *IpCq().trial_slack_x_U());

      SmartPtr<Vector> new_d_l = IpNLP().d_L()->MakeNew();
      IpNLP().Pd_L()->TransMultVector(1.0, *IpData().trial()->s(), 0.0, *new_d_l);
      new_d_l->Axpy(-1.0, *IpCq().trial_slack_s_L());

      SmartPtr<Vector> new_d_u = IpNLP().d_U()->MakeNew();
      IpNLP().Pd_U()->TransMultVector(1.0, *IpData().trial()->s(), 0.0, *new_d_u);
      new_d_u->Axpy(1.0, *IpCq().trial_slack_s_U());

      IpNLP().AdjustVariableBounds(*new_x_l, *new_x_u, *new_d_l, *new_d_u);

      if( Jnlst().ProduceOutput(J_VECTOR, J_MAIN) )
      {
         IpNLP().x_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_x_L");
         IpNLP().x_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_x_U");
         IpNLP().d_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_d_L");
         IpNLP().d_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_d_U");
      }

   }

   // Make sure that bound multipliers are not too far from \mu * S^{-1}
   // (see kappa_sigma in paper)
   bool corrected = false;
   Number max_correction;
   SmartPtr<const Vector> new_z_L;
   max_correction = correct_bound_multiplier(*IpData().trial()->z_L(), *IpCq().trial_slack_x_L(),
                    *IpCq().trial_compl_x_L(), new_z_L);
   if( max_correction > 0. )
   {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in z_L becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
   }
   SmartPtr<const Vector> new_z_U;
   max_correction = correct_bound_multiplier(*IpData().trial()->z_U(), *IpCq().trial_slack_x_U(),
                    *IpCq().trial_compl_x_U(), new_z_U);
   if( max_correction > 0. )
   {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in z_U becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
   }
   SmartPtr<const Vector> new_v_L;
   max_correction = correct_bound_multiplier(*IpData().trial()->v_L(), *IpCq().trial_slack_s_L(),
                    *IpCq().trial_compl_s_L(), new_v_L);
   if( max_correction > 0. )
   {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in v_L becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
   }
   SmartPtr<const Vector> new_v_U;
   max_correction = correct_bound_multiplier(*IpData().trial()->v_U(), *IpCq().trial_slack_s_U(),
                    *IpCq().trial_compl_s_U(), new_v_U);
   if( max_correction > 0. )
   {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in v_U becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
   }
   SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
   trial->Set_bound_mult(*new_z_L, *new_z_U, *new_v_L, *new_v_U);
   IpData().set_trial(trial);

   if( corrected )
   {
      IpData().Append_info_string("z");
   }

   // Accept the step
   IpData().AcceptTrialPoint();

   // If we want to recalculate the multipliers (e.g., as least
   // square estimates), call the calculator for that
   if( recalc_y_ )
   {
      // There is no point in doing this if there are no constraints
      if( IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim() == 0 )
      {
         recalc_y_ = false;
      }
   }
   if( recalc_y_ && IpCq().curr_constraint_violation() < recalc_y_feas_tol_ )
   {
      if( Jnlst().ProduceOutput(J_MOREDETAILED, J_MAIN) )
      {
         Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                        "dual infeasibility before least square multiplier update = %e\n",
                        IpCq().curr_dual_infeasibility(NORM_MAX));
      }
      IpData().Append_info_string("y ");
      DBG_ASSERT(IsValid(eq_multiplier_calculator_));
      if( IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim() > 0 )
      {
         SmartPtr<Vector> y_c = IpData().curr()->y_c()->MakeNew();
         SmartPtr<Vector> y_d = IpData().curr()->y_d()->MakeNew();
         bool retval = eq_multiplier_calculator_->CalculateMultipliers(*y_c, *y_d);
         if( retval )
         {
            SmartPtr<const IteratesVector> curr = IpData().curr();
            SmartPtr<IteratesVector> iterates = curr->MakeNewContainer();
            iterates->Set_x(*curr->x());
            iterates->Set_s(*curr->s());
            iterates->Set_z_L(*curr->z_L());
            iterates->Set_z_U(*curr->z_U());
            iterates->Set_v_L(*curr->v_L());
            iterates->Set_v_U(*curr->v_U());
            iterates->Set_y_c(*y_c);
            iterates->Set_y_d(*y_d);
            IpData().set_trial(iterates);
            IpData().AcceptTrialPoint();
         }
         else
         {
            Jnlst().Printf(J_DETAILED, J_MAIN,
                           "Recalculation of y multipliers skipped because eq_mult_calc returned false.\n");
         }
      }
   }
}

void IpoptAlgorithm::PrintProblemStatistics()
{
   if( !Jnlst().ProduceOutput(J_SUMMARY, J_STATISTICS) )
   {
      // nothing to print
      return;
   }

   Index nx_tot, nx_only_lower, nx_both, nx_only_upper;
   calc_number_of_bounds(*IpData().curr()->x(), *IpNLP().x_L(), *IpNLP().x_U(), *IpNLP().Px_L(), *IpNLP().Px_U(),
                         nx_tot, nx_only_lower, nx_both, nx_only_upper);

   Index ns_tot, ns_only_lower, ns_both, ns_only_upper;
   calc_number_of_bounds(*IpData().curr()->s(), *IpNLP().d_L(), *IpNLP().d_U(), *IpNLP().Pd_L(), *IpNLP().Pd_U(),
                         ns_tot, ns_only_lower, ns_both, ns_only_upper);

   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "Total number of variables............................: %8" IPOPT_INDEX_FORMAT "\n", nx_tot);
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "                     variables with only lower bounds: %8" IPOPT_INDEX_FORMAT "\n", nx_only_lower);
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "                variables with lower and upper bounds: %8" IPOPT_INDEX_FORMAT "\n", nx_both);
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "                     variables with only upper bounds: %8" IPOPT_INDEX_FORMAT "\n", nx_only_upper);
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "Total number of equality constraints.................: %8" IPOPT_INDEX_FORMAT "\n", IpData().curr()->y_c()->Dim());
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "Total number of inequality constraints...............: %8" IPOPT_INDEX_FORMAT "\n", ns_tot);
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "        inequality constraints with only lower bounds: %8" IPOPT_INDEX_FORMAT "\n", ns_only_lower);
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "   inequality constraints with lower and upper bounds: %8" IPOPT_INDEX_FORMAT "\n", ns_both);
   Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                  "        inequality constraints with only upper bounds: %8" IPOPT_INDEX_FORMAT "\n\n", ns_only_upper);
}

void IpoptAlgorithm::ComputeFeasibilityMultipliers()
{
   DBG_START_METH("IpoptAlgorithm::ComputeFeasibilityMultipliers",
                  dbg_verbosity);
   DBG_ASSERT(IpCq().IsSquareProblem());

   // if not primal feasible yet, then do not compute multipliers yet
   Number constr_viol = IpCq().unscaled_curr_nlp_constraint_violation(NORM_MAX);
   if( constr_viol > constr_viol_tol_ )
   {
      return;
   }

   // if we don't have an object for computing least square
   // multipliers we don't compute them
   if( IsNull(eq_multiplier_calculator_) )
   {
      Jnlst().Printf(J_DETAILED, J_SOLUTION,
                     "No eq_mult_calculator object available in IpoptAlgorithm to recompute multipliers at solution for square problem.\n");
      return;
   }

   IpData().TimingStats().CheckConvergence().Start();
   ConvergenceCheck::ConvergenceStatus conv_status = conv_check_->CheckConvergence(false);
   IpData().TimingStats().CheckConvergence().End();

   // if converged or reached some limit, then do not update multipliers
   // status CONTINUE likely means that we are not dual feasible yet, which we try to fix below
   if( conv_status != ConvergenceCheck::CONTINUE )
   {
      return;
   }

   // backup current iterate for case eq_mult_calculator fails
   // TODO we could avoid this backup&restore, if eq_multiplier_calculator_ could be told to calculate the multipliers for the trial instead of the current iterate
   SmartPtr<const IteratesVector> curr_backup = IpData().curr();

   SmartPtr<IteratesVector> iterates = IpData().curr()->MakeNewContainer();
   SmartPtr<Vector> tmp = iterates->z_L()->MakeNew();
   tmp->Set(0.);
   iterates->Set_z_L(*tmp);
   tmp = iterates->z_U()->MakeNew();
   tmp->Set(0.);
   iterates->Set_z_U(*tmp);
   tmp = iterates->v_L()->MakeNew();
   tmp->Set(0.);
   iterates->Set_v_L(*tmp);
   tmp = iterates->v_U()->MakeNew();
   tmp->Set(0.);
   iterates->Set_v_U(*tmp);
   SmartPtr<Vector> y_c = iterates->y_c()->MakeNew();
   SmartPtr<Vector> y_d = iterates->y_d()->MakeNew();
   IpData().set_trial(iterates);
   IpData().AcceptTrialPoint();
   bool retval = eq_multiplier_calculator_->CalculateMultipliers(*y_c, *y_d);
   if( retval )
   {
      //TODO Check if following line is really necessary
      iterates = IpData().curr()->MakeNewContainer();
      iterates->Set_y_c(*y_c);
      iterates->Set_y_d(*y_d);
      IpData().set_trial(iterates);
      IpData().AcceptTrialPoint();

      // check whether the new iterate satisfies convergence criteria now
      // if not, then we better continue with backed up iterate
      IpData().TimingStats().CheckConvergence().Start();
      ConvergenceCheck::ConvergenceStatus conv_status = conv_check_->CheckConvergence(false);
      IpData().TimingStats().CheckConvergence().End();

      if( conv_status == ConvergenceCheck::CONVERGED || conv_status == ConvergenceCheck::CONVERGED_TO_ACCEPTABLE_POINT )
      {
         return;
      }

      Jnlst().Printf(J_DETAILED, J_SOLUTION,
                     "Multipliers for feasibility problem using eq_mult_calculator does not lead to converged status yet.\n");
   }
   else
   {
      Jnlst().Printf(J_DETAILED, J_SOLUTION,
                     "Failed to compute multipliers for feasibility problem using eq_mult_calculator.\n");
   }

   // restore original iterate
   Jnlst().Printf(J_DETAILED, J_SOLUTION,
                  "Restoring iterate from before trying eq_mult_calculator.\n");
   SmartPtr<IteratesVector> orig_iterate = curr_backup->MakeNewContainer();
   IpData().set_trial(orig_iterate);
   IpData().AcceptTrialPoint();
}

void IpoptAlgorithm::ComputeFeasibilityMultipliersPostprocess()
{
   DBG_START_METH("IpoptAlgorithm::ComputeFeasibilityMultipliersPostprocess",
                  dbg_verbosity);
   DBG_ASSERT(IpCq().IsSquareProblem());

   // if we don't have an object for computing least square
   // multipliers we don't compute them
   if( IsNull(eq_multiplier_calculator_) )
   {
      Jnlst().Printf(J_WARNING, J_SOLUTION,
                     "No eq_mult_calculator object available in IpoptAlgorithm to recompute multipliers at solution for square problem.\n");
      return;
   }

   SmartPtr<IteratesVector> iterates = IpData().curr()->MakeNewContainer();
   SmartPtr<Vector> tmp = iterates->z_L()->MakeNew();
   tmp->Set(0.);
   iterates->Set_z_L(*tmp);
   tmp = iterates->z_U()->MakeNew();
   tmp->Set(0.);
   iterates->Set_z_U(*tmp);
   tmp = iterates->v_L()->MakeNew();
   tmp->Set(0.);
   iterates->Set_v_L(*tmp);
   tmp = iterates->v_U()->MakeNew();
   tmp->Set(0.);
   iterates->Set_v_U(*tmp);
   SmartPtr<Vector> y_c = iterates->y_c()->MakeNew();
   SmartPtr<Vector> y_d = iterates->y_d()->MakeNew();
   IpData().set_trial(iterates);
   IpData().AcceptTrialPoint();
   bool retval = eq_multiplier_calculator_->CalculateMultipliers(*y_c, *y_d);
   if( retval )
   {
      //TODO Check if following line is really necessary
      iterates = IpData().curr()->MakeNewContainer();
      iterates->Set_y_c(*y_c);
      iterates->Set_y_d(*y_d);
      IpData().set_trial(iterates);
      IpData().AcceptTrialPoint();
   }
   else
   {
      Jnlst().Printf(J_WARNING, J_SOLUTION,
                     "Failed to compute multipliers for feasibility problem using eq_mult_calculator.\n");
   }
}

void IpoptAlgorithm::calc_number_of_bounds(
   const Vector& x,
   const Vector& x_L,
   const Vector& x_U,
   const Matrix& Px_L,
   const Matrix& Px_U,
   Index&        n_tot,
   Index&        n_only_lower,
   Index&        n_both,
   Index&        n_only_upper
)
{
   DBG_START_METH("IpoptAlgorithm::calc_number_of_bounds",
                  dbg_verbosity);

   n_tot = x.Dim();

   SmartPtr<Vector> tmpx = x.MakeNew();
   SmartPtr<Vector> tmpxL = x_L.MakeNew();
   SmartPtr<Vector> tmpxU = x_U.MakeNew();

   tmpxL->Set(-1.);
   tmpxU->Set(2.);
   Px_L.MultVector(1.0, *tmpxL, 0.0, *tmpx);
   Px_U.MultVector(1.0, *tmpxU, 1.0, *tmpx);
   // Now, x has elements
   //  -1 : if component has only lower bound
   //   0 : if component has no bound
   //   1 : if component has both lower and upper bound
   //   2 : if component has only upper bound
   DBG_PRINT_VECTOR(2, "x-indicator", *tmpx);

   SmartPtr<Vector> tmpx0 = x.MakeNew();
   tmpx0->Set(0.);

   SmartPtr<Vector> tmpx2 = x.MakeNew();
   tmpx2->Set(-1.0);
   tmpx2->Axpy(1.0, *tmpx);
   tmpx2->ElementWiseMax(*tmpx0); // tmpx2 is now 1 in those
   // components with only upper bounds
   n_only_upper = (Index) tmpx2->Asum();

   tmpx->Axpy(-2., *tmpx2);       // now make all those entries for
   // only upper bounds zero in tmpx

   tmpx2->Copy(*tmpx);
   tmpx2->ElementWiseMax(*tmpx0); // tmpx2 is now 1 in those
   // components with both bounds
   n_both = (Index) tmpx2->Asum();

   tmpx->Axpy(-1., *tmpx2);
   tmpx->ElementWiseMin(*tmpx);   // tmpx is now -1 in those with only
   // lower bounds
   n_only_lower = (Index) tmpx->Asum();

}

Number IpoptAlgorithm::correct_bound_multiplier(
   const Vector&           trial_z,
   const Vector&           trial_slack,
   const Vector&           trial_compl,
   SmartPtr<const Vector>& new_trial_z
)
{
   DBG_START_METH("IpoptAlgorithm::CorrectBoundMultiplier",
                  dbg_verbosity);

   if( kappa_sigma_ < 1. || trial_z.Dim() == 0 )
   {
      new_trial_z = &trial_z;
      return 0.;
   }

   // We choose as barrier parameter to be used either the current
   // algorithmic barrier parameter (if we are not in the free mode),
   // or the average complementarity (at the trial point)
   Number mu;
   if( IpData().FreeMuMode() )
   {
      mu = IpCq().trial_avrg_compl();
      mu = Min(mu, Number(1e3));
   }
   else
   {
      mu = IpData().curr_mu();
   }
   DBG_PRINT((1, "mu = %8.2e\n", mu));
   DBG_PRINT_VECTOR(2, "trial_z", trial_z);

   // First check quickly if anything need to be corrected, using the
   // trial complementarity directly.  Here, Amax is the same as Max
   // (and we use Amax because that can be used later)
   if( trial_compl.Amax() <= kappa_sigma_ * mu && trial_compl.Min() >= 1. / kappa_sigma_ * mu )
   {
      new_trial_z = &trial_z;
      return 0.;
   }

   SmartPtr<Vector> one_over_s = trial_z.MakeNew();
   one_over_s->Copy(trial_slack);
   one_over_s->ElementWiseReciprocal();

   SmartPtr<Vector> step_z = trial_z.MakeNew();
   step_z->AddTwoVectors(kappa_sigma_ * mu, *one_over_s, -1., trial_z, 0.);

   DBG_PRINT_VECTOR(2, "step_z", *step_z);

   Number max_correction_up = Max(Number(0.), -step_z->Min());
   if( max_correction_up > 0. )
   {
      SmartPtr<Vector> tmp = trial_z.MakeNew();
      tmp->Set(0.);
      step_z->ElementWiseMin(*tmp);
      tmp->AddTwoVectors(1., trial_z, 1., *step_z, 0.);
      new_trial_z = GetRawPtr(tmp);
   }
   else
   {
      new_trial_z = &trial_z;
   }

   step_z->AddTwoVectors(1. / kappa_sigma_ * mu, *one_over_s, -1., *new_trial_z, 0.);

   Number max_correction_low = Max(Number(0.), step_z->Max());
   if( max_correction_low > 0. )
   {
      SmartPtr<Vector> tmp = trial_z.MakeNew();
      tmp->Set(0.);
      step_z->ElementWiseMax(*tmp);
      tmp->AddTwoVectors(1., *new_trial_z, 1., *step_z, 0.);
      new_trial_z = GetRawPtr(tmp);
   }

   DBG_PRINT_VECTOR(2, "new_trial_z", *new_trial_z);

   return Max(max_correction_up, max_correction_low);
}

void IpoptAlgorithm::print_copyright_message(
   const Journalist& jnlst
)
{
   jnlst.Printf(J_INSUPPRESSIBLE, J_MAIN,
                "\n******************************************************************************\n"
                "This program contains Ipopt, a library for large-scale nonlinear optimization.\n"
                " Ipopt is released as open source code under the Eclipse Public License (EPL).\n"
                "         For more information visit https://github.com/coin-or/Ipopt\n"
                "******************************************************************************\n\n");
   copyright_message_printed = true;
}

} // namespace Ipopt
