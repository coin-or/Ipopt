// Copyright (C) 2004, 2012 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02

#include "IpoptConfig.h"
#include "IpIpoptApplication.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpIpoptAlg.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpAlgBuilder.hpp"
#include "IpSolveStatistics.hpp"
#include "IpLinearSolversRegOp.hpp"
#include "IpInterfacesRegOp.hpp"
#include "IpAlgorithmRegOp.hpp"
#include "IpCGPenaltyRegOp.hpp"
#include "IpNLPBoundsRemover.hpp"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

#ifdef BUILD_INEXACT
# include "IpInexactRegOp.hpp"
# include "IpInexactAlgBuilder.hpp"
#endif

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include <fstream>

// Factory to facilitate creating IpoptApplication objects from within a DLL

Ipopt::IpoptApplication * IpoptApplicationFactory()
{
  return new Ipopt::IpoptApplication;
}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  IpoptApplication::IpoptApplication(bool create_console_out /* = true */,
                                     bool create_empty /* = false */)
      :
      read_params_dat_(true),
      rethrow_nonipoptexception_(false),
      inexact_algorithm_(false),
      replace_bounds_(false)
  {
    options_ = new OptionsList();
    if (create_empty)
      return;

    jnlst_ = new Journalist();
    try {
# if COIN_IPOPT_VERBOSITY > 0
      DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst_));
      SmartPtr<Journal> debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_ITERSUMMARY);
      debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
# endif

      DBG_START_METH("IpoptApplication::IpoptApplication()",
                     dbg_verbosity);

      if (create_console_out) {
        SmartPtr<Journal> stdout_jrnl =
          jnlst_->AddFileJournal("console", "stdout", J_ITERSUMMARY);
        stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
      }

      // Register the valid options
      reg_options_ = new RegisteredOptions();
      RegisterAllIpoptOptions(reg_options_);

      options_->SetJournalist(jnlst_);
      options_->SetRegisteredOptions(reg_options_);
    }
    catch (IpoptException& exc) {
      exc.ReportException(*jnlst_);
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR,
                      "Caught unknown Ipopt exception");
    }
    catch (std::bad_alloc&) {
      jnlst_->Printf(J_ERROR, J_MAIN, "\nEXIT: Not enough memory.\n");
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR,
                      "Not enough memory");
    }
    catch (...) {
      IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
      exc.ReportException(*jnlst_);
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR,
                      "Caught unknown exception");
    }
  }

  IpoptApplication::IpoptApplication(SmartPtr<RegisteredOptions> reg_options,
                                     SmartPtr<OptionsList> options,
                                     SmartPtr<Journalist> jnlst)
      :
      read_params_dat_(true),
      rethrow_nonipoptexception_(false),
      jnlst_(jnlst),
      reg_options_(reg_options),
      options_(options),
      inexact_algorithm_(false),
      replace_bounds_(false)
  {}

  SmartPtr<IpoptApplication> IpoptApplication::clone()
  {
    SmartPtr<IpoptApplication> retval = new IpoptApplication(false, true);
    retval->jnlst_ = Jnlst();
    retval->reg_options_ = RegOptions();
    *retval->options_ = *Options();

    retval->read_params_dat_ = read_params_dat_;
    retval->inexact_algorithm_ = inexact_algorithm_;
    retval->replace_bounds_ = replace_bounds_;
    retval->rethrow_nonipoptexception_ = rethrow_nonipoptexception_;

    return retval;
  }

  ApplicationReturnStatus
  IpoptApplication::Initialize(std::istream& is, bool allow_clobber)
  {
    try {
      // Get the options
      if (is.good()) {
        // stream exists, read the content
        options_->ReadFromStream(*jnlst_, is, allow_clobber);
      }

      bool no_output;
      options_->GetBoolValue("suppress_all_output", no_output, "");

      if (no_output) {
        jnlst_->DeleteAllJournals();
      }
      else {
        Index ivalue;
        options_->GetIntegerValue("print_level", ivalue, "");
        EJournalLevel print_level = (EJournalLevel)ivalue;
        SmartPtr<Journal> stdout_jrnl = jnlst_->GetJournal("console");
        if (IsValid(stdout_jrnl)) {
          // Set printlevel for stdout
          stdout_jrnl->SetAllPrintLevels(print_level);
          stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
        }

        bool option_set;

#if COIN_IPOPT_VERBOSITY > 0
        // Set printlevel for debug
        option_set = options_->GetIntegerValue("debug_print_level",
                                               ivalue, "");
        EJournalLevel debug_print_level;
        if (option_set) {
          debug_print_level = (EJournalLevel)ivalue;
        }
        else {
          debug_print_level = print_level;
        }
        SmartPtr<Journal> debug_jrnl = jnlst_->GetJournal("Debug");
        if (IsNull(debug_jrnl)) {
          debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_ITERSUMMARY);
        }
        debug_jrnl->SetAllPrintLevels(debug_print_level);
        debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
#endif

        // Open an output file if required
        std::string output_filename;
        options_->GetStringValue("output_file", output_filename, "");
        if (output_filename != "") {
          EJournalLevel file_print_level;
          option_set = options_->GetIntegerValue("file_print_level", ivalue, "");
          if (option_set) {
            file_print_level = (EJournalLevel)ivalue;
          }
          else {
            file_print_level = print_level;
          }
          bool openend = OpenOutputFile(output_filename, file_print_level);
          if (!openend) {
            jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                           "Error opening output file \"%s\"\n",
                           output_filename.c_str());
            return Invalid_Option;
          }
        }
      }

      // output a description of all the options
      bool print_options_documentation;
      options_->GetBoolValue("print_options_documentation",
                             print_options_documentation, "");
      if (print_options_documentation) {
        bool latex;
        options_->GetBoolValue("print_options_latex_mode", latex, "");
        if (latex) {
          std::list<std::string> options_to_print;
          options_to_print.push_back("#Output");
          options_to_print.push_back("print_level");
          options_to_print.push_back("print_user_options");
          options_to_print.push_back("print_options_documentation");
          options_to_print.push_back("print_frequency_iter");
          options_to_print.push_back("print_frequency_time");
          options_to_print.push_back("output_file");
          options_to_print.push_back("file_print_level");
          options_to_print.push_back("option_file_name");
          options_to_print.push_back("print_info_string");
          options_to_print.push_back("inf_pr_output");
          options_to_print.push_back("print_timing_statistics");

          options_to_print.push_back("#Termination");
          options_to_print.push_back("tol");
          options_to_print.push_back("max_iter");
          options_to_print.push_back("max_cpu_time");
          options_to_print.push_back("dual_inf_tol");
          options_to_print.push_back("constr_viol_tol");
          options_to_print.push_back("compl_inf_tol");
          options_to_print.push_back("acceptable_tol");
          options_to_print.push_back("acceptable_iter");
          options_to_print.push_back("acceptable_constr_viol_tol");
          options_to_print.push_back("acceptable_dual_inf_tol");
          options_to_print.push_back("acceptable_compl_inf_tol");
          options_to_print.push_back("acceptable_obj_change_tol");
          options_to_print.push_back("diverging_iterates_tol");

          options_to_print.push_back("#NLP Scaling");
          options_to_print.push_back("obj_scaling_factor");
          options_to_print.push_back("nlp_scaling_method");
          options_to_print.push_back("nlp_scaling_max_gradient");
          options_to_print.push_back("nlp_scaling_min_value");

          options_to_print.push_back("#NLP");
          options_to_print.push_back("bound_relax_factor");
          options_to_print.push_back("honor_original_bounds");
          options_to_print.push_back("check_derivatives_for_naninf");
          options_to_print.push_back("nlp_lower_bound_inf");
          options_to_print.push_back("nlp_upper_bound_inf");
          options_to_print.push_back("fixed_variable_treatment");
          options_to_print.push_back("jac_c_constant");
          options_to_print.push_back("jac_d_constant");
          options_to_print.push_back("hessian_constant");

          options_to_print.push_back("#Initialization");
          options_to_print.push_back("bound_frac");
          options_to_print.push_back("bound_push");
          options_to_print.push_back("slack_bound_frac");
          options_to_print.push_back("slack_bound_push");
          options_to_print.push_back("bound_mult_init_val");
          options_to_print.push_back("constr_mult_init_max");
          options_to_print.push_back("bound_mult_init_method");

          options_to_print.push_back("#Barrier Parameter");
          options_to_print.push_back("mehrotra_algorithm");
          options_to_print.push_back("mu_strategy");
          options_to_print.push_back("mu_oracle");
          options_to_print.push_back("quality_function_max_section_steps");
          options_to_print.push_back("fixed_mu_oracle");
          options_to_print.push_back("adaptive_mu_globalization");
          options_to_print.push_back("mu_init");
          options_to_print.push_back("mu_max_fact");
          options_to_print.push_back("mu_max");
          options_to_print.push_back("mu_min");
          options_to_print.push_back("mu_target");
          options_to_print.push_back("barrier_tol_factor");
          options_to_print.push_back("mu_linear_decrease_factor");
          options_to_print.push_back("mu_superlinear_decrease_power");

          options_to_print.push_back("#Multiplier Updates");
          options_to_print.push_back("alpha_for_y");
          options_to_print.push_back("alpha_for_y_tol");
          options_to_print.push_back("recalc_y");
          options_to_print.push_back("recalc_y_feas_tol");

          options_to_print.push_back("#Line Search");
          options_to_print.push_back("max_soc");
          options_to_print.push_back("watchdog_shortened_iter_trigger");
          options_to_print.push_back("watchdog_trial_iter_max");
          options_to_print.push_back("accept_every_trial_step");
          options_to_print.push_back("corrector_type");
          options_to_print.push_back("soc_method");

          options_to_print.push_back("#Warm Start");
          options_to_print.push_back("warm_start_init_point");
          options_to_print.push_back("warm_start_bound_push");
          options_to_print.push_back("warm_start_bound_frac");
          options_to_print.push_back("warm_start_slack_bound_frac");
          options_to_print.push_back("warm_start_slack_bound_push");
          options_to_print.push_back("warm_start_mult_bound_push");
          options_to_print.push_back("warm_start_mult_init_max");

          options_to_print.push_back("#Restoration Phase");
          options_to_print.push_back("expect_infeasible_problem");
          options_to_print.push_back("expect_infeasible_problem_ctol");
          options_to_print.push_back("expect_infeasible_problem_ytol");
          options_to_print.push_back("start_with_resto");
          options_to_print.push_back("soft_resto_pderror_reduction_factor");
          options_to_print.push_back("required_infeasibility_reduction");
          options_to_print.push_back("bound_mult_reset_threshold");
          options_to_print.push_back("constr_mult_reset_threshold");
          options_to_print.push_back("evaluate_orig_obj_at_resto_trial");

          options_to_print.push_back("#Linear Solver");
          options_to_print.push_back("linear_solver");
          options_to_print.push_back("linear_system_scaling");
          options_to_print.push_back("linear_scaling_on_demand");
          options_to_print.push_back("max_refinement_steps");
          options_to_print.push_back("min_refinement_steps");
          options_to_print.push_back("neg_curv_test_reg");
          options_to_print.push_back("neg_curv_test_tol");

          options_to_print.push_back("#Hessian Perturbation");
          options_to_print.push_back("max_hessian_perturbation");
          options_to_print.push_back("min_hessian_perturbation");
          options_to_print.push_back("first_hessian_perturbation");
          options_to_print.push_back("perturb_inc_fact_first");
          options_to_print.push_back("perturb_inc_fact");
          options_to_print.push_back("perturb_dec_fact");
          options_to_print.push_back("jacobian_regularization_value");

          options_to_print.push_back("#Quasi-Newton");
          options_to_print.push_back("hessian_approximation");
          options_to_print.push_back("limited_memory_update_type");
          options_to_print.push_back("limited_memory_max_history");
          options_to_print.push_back("limited_memory_max_skipping");
          options_to_print.push_back("limited_memory_initialization");
          options_to_print.push_back("limited_memory_init_val");
          options_to_print.push_back("limited_memory_init_val_max");
          options_to_print.push_back("limited_memory_init_val_min");
          options_to_print.push_back("limited_memory_special_for_resto");

          options_to_print.push_back("#Derivative Test");
          options_to_print.push_back("derivative_test");
          options_to_print.push_back("derivative_test_perturbation");
          options_to_print.push_back("derivative_test_tol");
          options_to_print.push_back("derivative_test_print_all");
          options_to_print.push_back("derivative_test_first_index");
          options_to_print.push_back("point_perturbation_radius");

          // Special linear solver
#if defined(COINHSL_HAS_MA27) || defined(HAVE_LINEARSOLVERLOADER)

          options_to_print.push_back("#MA27 Linear Solver");
          options_to_print.push_back("ma27_pivtol");
          options_to_print.push_back("ma27_pivtolmax");
          options_to_print.push_back("ma27_liw_init_factor");
          options_to_print.push_back("ma27_la_init_factor");
          options_to_print.push_back("ma27_meminc_factor");
#endif

#if defined(COINHSL_HAS_MA57) || defined(HAVE_LINEARSOLVERLOADER)

          options_to_print.push_back("#MA57 Linear Solver");
          options_to_print.push_back("ma57_pivtol");
          options_to_print.push_back("ma57_pivtolmax");
          options_to_print.push_back("ma57_pre_alloc");
          options_to_print.push_back("ma57_pivot_order");
          options_to_print.push_back("ma57_automatic_scaling");
          options_to_print.push_back("ma57_block_size");
          options_to_print.push_back("ma57_node_amalgamation");
          options_to_print.push_back("ma57_small_pivot_flag");
#endif

#if defined(COINHSL_HAS_MA77) || defined(HAVE_LINEARSOLVERLOADER)

          options_to_print.push_back("#MA77 Linear Solver");
          options_to_print.push_back("ma77_print_level");
          options_to_print.push_back("ma77_buffer_lpage");
          options_to_print.push_back("ma77_buffer_npage");
          options_to_print.push_back("ma77_file_size");
          options_to_print.push_back("ma77_maxstore");
          options_to_print.push_back("ma77_nemin");
#if defined(COINHSL_HSL2013) || defined(HAVE_LINEARSOLVERLOADER)
          options_to_print.push_back("ma77_order");
#endif
          options_to_print.push_back("ma77_small");
          options_to_print.push_back("ma77_static");
          options_to_print.push_back("ma77_u");
          options_to_print.push_back("ma77_umax");
#endif

#if defined(COINHSL_HAS_MA86) || defined(HAVE_LINEARSOLVERLOADER)

          options_to_print.push_back("#MA86 Linear Solver");
          options_to_print.push_back("ma86_print_level");
          options_to_print.push_back("ma86_nemin");
          options_to_print.push_back("ma86_order");
          options_to_print.push_back("ma86_scaling");
          options_to_print.push_back("ma86_small");
          options_to_print.push_back("ma86_static");
          options_to_print.push_back("ma86_u");
          options_to_print.push_back("ma86_umax");
#endif

#if defined(COINHSL_HAS_MA97) || defined(HAVE_LINEARSOLVERLOADER)

          options_to_print.push_back("#MA97 Linear Solver");
          options_to_print.push_back("ma97_print_level");
          options_to_print.push_back("ma97_nemin");
          options_to_print.push_back("ma97_order");
          options_to_print.push_back("ma97_scaling");
          options_to_print.push_back("ma97_scaling1");
          options_to_print.push_back("ma97_scaling2");
          options_to_print.push_back("ma97_scaling3");
          options_to_print.push_back("ma97_small");
          options_to_print.push_back("ma97_solve_blas3");
          options_to_print.push_back("ma97_switch1");
          options_to_print.push_back("ma97_switch2");
          options_to_print.push_back("ma97_switch3");
          options_to_print.push_back("ma97_u");
          options_to_print.push_back("ma97_umax");
#endif

#ifdef COIN_HAS_MUMPS

          options_to_print.push_back("#MUMPS Linear Solver");
          options_to_print.push_back("mumps_pivtol");
          options_to_print.push_back("mumps_pivtolmax");
          options_to_print.push_back("mumps_mem_percent");
          options_to_print.push_back("mumps_permuting_scaling");
          options_to_print.push_back("mumps_pivot_order");
          options_to_print.push_back("mumps_scaling");
#endif

#if defined(HAVE_PARDISO) || defined(HAVE_LINEARSOLVERLOADER)

          options_to_print.push_back("#Pardiso Linear Solver");
          options_to_print.push_back("pardiso_matching_strategy");
          options_to_print.push_back("pardiso_max_iterative_refinement_steps");
          options_to_print.push_back("pardiso_msglvl");
          options_to_print.push_back("pardiso_order");
          //options_to_print.push_back("pardiso_out_of_core_power");
#endif

#ifdef HAVE_WSMP

          options_to_print.push_back("#WSMP Linear Solver");
          options_to_print.push_back("wsmp_num_threads");
          options_to_print.push_back("wsmp_ordering_option");
          options_to_print.push_back("wsmp_pivtol");
          options_to_print.push_back("wsmp_pivtolmax");
          options_to_print.push_back("wsmp_scaling");
          options_to_print.push_back("wsmp_singularity_threshold");
#endif

          reg_options_->OutputLatexOptionDocumentation(*jnlst_, options_to_print);
        }
        else {
          std::list<std::string> categories;
          categories.push_back("Output");
          /*categories.push_back("Main Algorithm");*/
          categories.push_back("Convergence");
          categories.push_back("NLP Scaling");
          categories.push_back("NLP");
          categories.push_back("Initialization");
          categories.push_back("Barrier Parameter Update");
          categories.push_back("Line Search");
          categories.push_back("Warm Start");
          categories.push_back("Linear Solver");
          categories.push_back("Step Calculation");
          categories.push_back("Restoration Phase");
          categories.push_back("Derivative Checker");
          categories.push_back("Hessian Approximation");
          categories.push_back("MA27 Linear Solver");
          categories.push_back("MA57 Linear Solver");
          categories.push_back("Pardiso Linear Solver");
#ifdef HAVE_WSMP

          categories.push_back("WSMP Linear Solver");
#endif
#ifdef COIN_HAS_MUMPS

          categories.push_back("Mumps Linear Solver");
#endif
          categories.push_back("MA28 Linear Solver");

          categories.push_back("Uncategorized");
          //categories.push_back("Undocumented Options");
          reg_options_->OutputOptionDocumentation(*jnlst_, categories);
        }
      }

#ifdef BUILD_INEXACT
      // Check if we are to use the inexact linear solver option
      options_->GetBoolValue("inexact_algorithm", inexact_algorithm_, "");
      // Change the default flags for the inexact algorithm
      if (inexact_algorithm_) {
        AddInexactDefaultOptions(*options_);
      }
#endif

      options_->GetBoolValue("replace_bounds", replace_bounds_, "");
    }
    catch (OPTION_INVALID& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      return Invalid_Option;
    }
    catch (IpoptException& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      return Unrecoverable_Exception;
    }
    catch (std::bad_alloc&) {
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
      return Insufficient_Memory;
    }
    catch (...) {
      if( !rethrow_nonipoptexception_ )
      {
        IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
        exc.ReportException(*jnlst_);
        return NonIpopt_Exception_Thrown;
      }
      else
      {
        throw;
      }
    }
    return Solve_Succeeded;
  }

  ApplicationReturnStatus
  IpoptApplication::Initialize(std::string params_file, bool allow_clobber)
  {
    std::ifstream is;
    if (params_file != "") {
      try {
        is.open(params_file.c_str());
      }
      catch (std::bad_alloc&) {
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
        return Insufficient_Memory;
      }
      catch (...) {
        if( !rethrow_nonipoptexception_ )
        {
          IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
          exc.ReportException(*jnlst_);
          return NonIpopt_Exception_Thrown;
        }
        else
        {
          throw;
        }
      }
    }
    ApplicationReturnStatus retval = Initialize(is, allow_clobber);
    if (is) {
      is.close();
    }
    return retval;
  }

  ApplicationReturnStatus
  IpoptApplication::Initialize(bool allow_clobber)
  {
     std::string option_file_name;
     options_->GetStringValue("option_file_name", option_file_name, "");
     if (option_file_name != "" && option_file_name != "ipopt.opt")
        jnlst_->Printf(J_SUMMARY, J_MAIN, "Using option file \"%s\".\n\n", option_file_name.c_str());

     return Initialize(option_file_name, allow_clobber);
  }

  IpoptApplication::~IpoptApplication()
  {
    DBG_START_METH("IpoptApplication::~IpoptApplication()",
                   dbg_verbosity);
  }

  void IpoptApplication::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Output");
    roptions->AddBoundedIntegerOption(
      "print_level",
      "Output verbosity level.",
      0, J_LAST_LEVEL-1, J_ITERSUMMARY,
      "Sets the default verbosity level for console output. The "
      "larger this value the more detailed is the output.");

    roptions->AddStringOption1(
      "output_file",
      "File name of desired output file (leave unset for no file output).",
      "",
      "*", "Any acceptable standard file name",
      "NOTE: This option only works when read from the ipopt.opt options file! "
      "An output file with this name will be written (leave unset for no "
      "file output).  The verbosity level is by default set to \"print_level\", "
      "but can be overridden with \"file_print_level\".  The file name is "
      "changed to use only small letters.");
    roptions->AddBoundedIntegerOption(
      "file_print_level",
      "Verbosity level for output file.",
      0, J_LAST_LEVEL-1, J_ITERSUMMARY,
      "NOTE: This option only works when read from the ipopt.opt options file! "
      "Determines the verbosity level for the file specified by "
      "\"output_file\".  By default it is the same as \"print_level\".");
    roptions->AddStringOption2(
      "print_user_options",
      "Print all options set by the user.",
      "no",
      "no", "don't print options",
      "yes", "print options",
      "If selected, the algorithm will print the list of all options set by "
      "the user including their values and whether they have been used.  In "
      "some cases this information might be incorrect, due to the internal "
      "program flow.");
    roptions->AddStringOption2(
      "print_options_documentation",
      "Switch to print all algorithmic options.",
      "no",
      "no", "don't print list",
      "yes", "print list",
      "If selected, the algorithm will print the list of all available "
      "algorithmic options with some documentation before solving the "
      "optimization problem.");

#if COIN_IPOPT_VERBOSITY > 0

    roptions->AddBoundedIntegerOption(
      "debug_print_level",
      "Verbosity level for debug file.",
      0, J_LAST_LEVEL-1, J_ITERSUMMARY,
      "This Ipopt library has been compiled in debug mode, and a file "
      "\"debug.out\" is produced for every run.  This option determines "
      "the verbosity level for this file.  By default it is the same as "
      "\"print_level\".");

#endif

    roptions->AddStringOption2(
      "print_timing_statistics",
      "Switch to print timing statistics.",
      "no",
      "no", "don't print statistics",
      "yes", "print all timing statistics",
      "If selected, the program will print the CPU usage (user time) for "
      "selected tasks.");

    roptions->AddStringOption1(
      "option_file_name",
      "File name of options file.",
      "ipopt.opt",
      "*", "Any acceptable standard file name",
      "By default, the name of the Ipopt options file is \"ipopt.opt\" - or "
      "something else if specified in the IpoptApplication::Initialize call. "
      "If this option is set by SetStringValue BEFORE the options file is "
      "read, it specifies the name of the options file.  It does not make any "
      "sense to specify this option within the options file. "
      "Setting this option to an empty string disables reading of an options file.");

    roptions->AddStringOption2(
      "replace_bounds",
      "Indicates if all variable bounds should be replaced by inequality constraints", "no",
      "no", "leave bounds on variables",
      "yes", "replace variable bounds by inequality constraints",
      "This option must be set for the inexact algorithm");
    roptions->AddStringOption2(
      "skip_finalize_solution_call",
      "Indicates if call to NLP::FinalizeSolution after optimization should be suppressed", "no",
      "no", "call FinalizeSolution",
      "yes", "do not call FinalizeSolution",
      "In some Ipopt applications, the user might want to call the "
      "FinalizeSolution method separately.  Setting this option to \"yes\" "
      "will cause the IpoptApplication object to suppress the default call to "
      "that method.");

    roptions->SetRegisteringCategory("Undocumented");
    roptions->AddStringOption2(
      "print_options_latex_mode",
      "Undocumented", "no",
      "no", "Undocumented",
      "yes", "Undocumented",
      "Undocumented");
    roptions->AddStringOption2(
      "suppress_all_output",
      "Undocumented", "no",
      "no", "Undocumented",
      "yes", "Undocumented",
      "Undocumented");
#ifdef BUILD_INEXACT
    roptions->AddStringOption2(
      "inexact_algorithm",
      "Activate the version of Ipopt that allows iterative linear solvers.",
      "no",
      "no", "use default algorithm with direct linear solvers",
      "yes", "use the EXPERIMENTAL iterative linear solver option",
      "");
#endif
  }

  ApplicationReturnStatus
  IpoptApplication::OptimizeTNLP(const SmartPtr<TNLP>& tnlp)
  {
    nlp_adapter_ = new TNLPAdapter(GetRawPtr(tnlp), ConstPtr(jnlst_));
    return OptimizeNLP(nlp_adapter_);
  }

  ApplicationReturnStatus
  IpoptApplication::ReOptimizeTNLP(const SmartPtr<TNLP>& tnlp)
  {
    ASSERT_EXCEPTION(IsValid(nlp_adapter_), INVALID_WARMSTART,
                     "ReOptimizeTNLP called before OptimizeTNLP.");
    TNLPAdapter* adapter =
      static_cast<TNLPAdapter*> (GetRawPtr(nlp_adapter_));
    DBG_ASSERT(dynamic_cast<TNLPAdapter*> (GetRawPtr(nlp_adapter_)));
    ASSERT_EXCEPTION(adapter->tnlp()==tnlp, INVALID_WARMSTART,
                     "ReOptimizeTNLP called for different TNLP.")

    return ReOptimizeNLP(nlp_adapter_);
  }

  ApplicationReturnStatus
  IpoptApplication::OptimizeNLP(const SmartPtr<NLP>& nlp)
  {
    SmartPtr<AlgorithmBuilder> alg_builder = NULL;
    return OptimizeNLP(nlp, alg_builder);
  }

  ApplicationReturnStatus
  IpoptApplication::OptimizeNLP(const SmartPtr<NLP>& nlp, SmartPtr<AlgorithmBuilder>& alg_builder)
  {
    ApplicationReturnStatus retValue = Internal_Error;

    // Prepare internal data structures of the algorithm
    try {

      if (IsNull(alg_builder)) {
#ifdef BUILD_INEXACT
        if (inexact_algorithm_) {
          alg_builder = new InexactAlgorithmBuilder();
        }
        else {
#endif
          alg_builder = new AlgorithmBuilder();
#ifdef BUILD_INEXACT
        }
#endif
      }

      SmartPtr<NLP> use_nlp;
      if (replace_bounds_) {
        use_nlp = new NLPBoundsRemover(*nlp);
      }
      else {
        use_nlp = nlp;
      }
      alg_builder->BuildIpoptObjects(*jnlst_, *options_, "", use_nlp,
                                     ip_nlp_, ip_data_, ip_cq_);

      alg_ = GetRawPtr(alg_builder->BuildBasicAlgorithm(*jnlst_, *options_, ""));

      // finally call the optimization
      retValue = call_optimize();
    }
    catch (OPTION_INVALID& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid option encountered.\n");
      retValue = Invalid_Option;
    }
    catch (IpoptException& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Some uncaught Ipopt exception encountered.\n");
      retValue = Unrecoverable_Exception;
    }
    catch (std::bad_alloc&) {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
    }
    catch (...) {
      if( !rethrow_nonipoptexception_ )
      {
        IpoptException exc("Unknown Exception caught in Ipopt", "Unknown File", -1);
        exc.ReportException(*jnlst_, J_ERROR);
        retValue = NonIpopt_Exception_Thrown;
      }
      else
      {
        throw;
      }
    }

    jnlst_->FlushBuffer();

    return retValue;
  }

  ApplicationReturnStatus
  IpoptApplication::ReOptimizeNLP(const SmartPtr<NLP>& nlp)
  {
    ASSERT_EXCEPTION(IsValid(alg_), INVALID_WARMSTART,
                     "ReOptimizeNLP called before OptimizeNLP.");
    OrigIpoptNLP* orig_nlp =
      static_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_));
    DBG_ASSERT(dynamic_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_)));
    ASSERT_EXCEPTION(orig_nlp->nlp()==nlp, INVALID_WARMSTART,
                     "ReOptimizeTNLP called for different NLP.")

    return call_optimize();
  }


  ApplicationReturnStatus IpoptApplication::call_optimize()
  {
    // Reset the print-level for the screen output
    Index ivalue;
    options_->GetIntegerValue("print_level", ivalue, "");
    EJournalLevel print_level = (EJournalLevel)ivalue;
    SmartPtr<Journal> stdout_jrnl = jnlst_->GetJournal("console");
    if (IsValid(stdout_jrnl)) {
      // Set printlevel for stdout
      stdout_jrnl->SetAllPrintLevels(print_level);
      stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
    }

    statistics_ = NULL; /* delete old statistics */
    // Get the pointers to the real objects (need to do it that
    // awkwardly since otherwise we would have to include so many
    // things in IpoptApplication, which a user would have to
    // include, too (SmartPtr doesn't work otherwise on AIX)
    IpoptAlgorithm* p2alg = static_cast<IpoptAlgorithm*> (GetRawPtr(alg_));
    DBG_ASSERT(dynamic_cast<IpoptAlgorithm*> (GetRawPtr(alg_)));
    IpoptData* p2ip_data = static_cast<IpoptData*> (GetRawPtr(ip_data_));
    DBG_ASSERT(dynamic_cast<IpoptData*> (GetRawPtr(ip_data_)));
    OrigIpoptNLP* p2ip_nlp = static_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_));
    DBG_ASSERT(dynamic_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_)));
    IpoptCalculatedQuantities* p2ip_cq = static_cast<IpoptCalculatedQuantities*> (GetRawPtr(ip_cq_));
    DBG_ASSERT(dynamic_cast<IpoptCalculatedQuantities*> (GetRawPtr(ip_cq_)));

    // Reset Timing statistics
    ip_data_->TimingStats().ResetTimes();
    p2ip_nlp->ResetTimes();

    ApplicationReturnStatus retValue = Internal_Error;
    SolverReturn status = INTERNAL_ERROR;
    /** Flag indicating if the NLP:FinalizeSolution method should not
     *  be called after optimization. */
    bool skip_finalize_solution_call = false;
    try {

      // Set up the algorithm
      p2alg->Initialize(*jnlst_, *p2ip_nlp, *p2ip_data, *p2ip_cq,
                        *options_, "");

      // Process the options used below
      bool print_timing_statistics;
      options_->GetBoolValue("print_timing_statistics",
                             print_timing_statistics, "");

      // If selected, print the user options
      bool print_user_options;
      options_->GetBoolValue("print_user_options", print_user_options, "");
      if (print_user_options) {
        std::string liststr;
        options_->PrintUserOptions(liststr);
        jnlst_->Printf(J_ERROR, J_MAIN, "\nList of user-set options:\n\n%s", liststr.c_str());
      }

      if ( jnlst_->ProduceOutput(J_DETAILED, J_MAIN) ) {
        // Print out the options (including the number of times they were used
        std::string liststr;
        options_->PrintList(liststr);
        jnlst_->Printf(J_DETAILED, J_MAIN, "\nList of options:\n\n%s", liststr.c_str());
      }

      // Run the algorithm
      status = p2alg->Optimize();

      // Since all the output below doesn't make any sense in this
      // case, we rethrow the TOO_FEW_DOF exception here
      ASSERT_EXCEPTION(status != TOO_FEW_DEGREES_OF_FREEDOM, TOO_FEW_DOF,
                       "Too few degrees of freedom (rethrown)!");

      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "\nNumber of Iterations....: %d\n",
                     p2ip_data->iter_count());

      if (status!=INVALID_NUMBER_DETECTED) {
        try {
          jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                         "\n                                   (scaled)                 (unscaled)\n");
          jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                         "Objective...............: %24.16e  %24.16e\n",
                         p2ip_cq->curr_f(),
                         p2ip_cq->unscaled_curr_f());
          jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                         "Dual infeasibility......: %24.16e  %24.16e\n",
                         p2ip_cq->curr_dual_infeasibility(NORM_MAX),
                         p2ip_cq->unscaled_curr_dual_infeasibility(NORM_MAX));
          jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                         "Constraint violation....: %24.16e  %24.16e\n",
                         p2ip_cq->curr_nlp_constraint_violation(NORM_MAX),
                         p2ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX));
          jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                         "Complementarity.........: %24.16e  %24.16e\n",
                         p2ip_cq->curr_complementarity(0., NORM_MAX),
                         p2ip_cq->unscaled_curr_complementarity(0., NORM_MAX));
          jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                         "Overall NLP error.......: %24.16e  %24.16e\n\n",
                         p2ip_cq->curr_nlp_error(),
                         p2ip_cq->unscaled_curr_nlp_error());
        }
        catch (IpoptNLP::Eval_Error& exc) {
          status = INVALID_NUMBER_DETECTED;
          exc.ReportException(*jnlst_, J_ERROR);
        }
      }

      p2ip_data->curr()->x()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "x");
      p2ip_data->curr()->y_c()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_c");
      p2ip_data->curr()->y_d()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_d");
      p2ip_data->curr()->z_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_L");
      p2ip_data->curr()->z_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_U");
      p2ip_data->curr()->v_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_L");
      p2ip_data->curr()->v_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_U");

      if (status==LOCAL_INFEASIBILITY) {
        p2ip_cq->curr_c()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "curr_c");
        p2ip_cq->curr_d_minus_s()->Print(*jnlst_, J_VECTOR, J_SOLUTION,
                                         "curr_d_minus_s");
      }

      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "\nNumber of objective function evaluations             = %d\n",
                     p2ip_nlp->f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of objective gradient evaluations             = %d\n",
                     p2ip_nlp->grad_f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of equality constraint evaluations            = %d\n",
                     p2ip_nlp->c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of inequality constraint evaluations          = %d\n",
                     p2ip_nlp->d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of equality constraint Jacobian evaluations   = %d\n",
                     p2ip_nlp->jac_c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of inequality constraint Jacobian evaluations = %d\n",
                     p2ip_nlp->jac_d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of Lagrangian Hessian evaluations             = %d\n",
                     p2ip_nlp->h_evals());
      Number cpu_time_overall_alg = p2ip_data->TimingStats().OverallAlgorithm().TotalCpuTime();
      Number cpu_time_funcs = p2ip_nlp->TotalFunctionEvaluationCpuTime();
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Total CPU secs in IPOPT (w/o function evaluations)   = %10.3f\n",
                     cpu_time_overall_alg-cpu_time_funcs);
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Total CPU secs in NLP function evaluations           = %10.3f\n",
                     cpu_time_funcs);

      // Write timing statistics information
      if (print_timing_statistics) {
        jnlst_->Printf(J_SUMMARY, J_TIMING_STATISTICS,
                       "\n\nTiming Statistics:\n\n");
        p2ip_data->TimingStats().PrintAllTimingStatistics(*jnlst_, J_SUMMARY,
            J_TIMING_STATISTICS);
        p2ip_nlp->PrintTimingStatistics(*jnlst_, J_SUMMARY,
                                        J_TIMING_STATISTICS);
      }

      // Write EXIT message
      if (status == SUCCESS) {
        retValue = Solve_Succeeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Optimal Solution Found.\n");
      }
      else if (status == MAXITER_EXCEEDED) {
        retValue = Maximum_Iterations_Exceeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum Number of Iterations Exceeded.\n");
      }
      else if (status == CPUTIME_EXCEEDED) {
        retValue = Maximum_CpuTime_Exceeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum CPU time exceeded.\n");
      }
      else if (status == STOP_AT_TINY_STEP) {
        retValue = Search_Direction_Becomes_Too_Small;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Search Direction is becoming Too Small.\n");
      }
      else if (status == STOP_AT_ACCEPTABLE_POINT) {
        retValue = Solved_To_Acceptable_Level;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Solved To Acceptable Level.\n");
      }
      else if (status == FEASIBLE_POINT_FOUND) {
        retValue = Feasible_Point_Found;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Feasible point for square problem found.\n");
      }
      else if (status == DIVERGING_ITERATES) {
        retValue = Diverging_Iterates;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Iterates diverging; problem might be unbounded.\n");
      }
      else if (status == RESTORATION_FAILURE) {
        retValue = Restoration_Failed;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Restoration Failed!\n");
      }
      else if (status == ERROR_IN_STEP_COMPUTATION) {
        retValue = Error_In_Step_Computation;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Error in step computation (regularization becomes too large?)!\n");
      }
      else if (status == LOCAL_INFEASIBILITY) {
        retValue = Infeasible_Problem_Detected;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Converged to a point of local infeasibility. Problem may be infeasible.\n");
      }
      else if (status == USER_REQUESTED_STOP) {
        retValue = User_Requested_Stop;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Stopping optimization at current point as requested by user.\n");
      }
      else if (status == INVALID_NUMBER_DETECTED) {
        retValue = Invalid_Number_Detected;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid number in NLP function or derivative detected.\n");
      }
      else {
        retValue = Internal_Error;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
        return retValue;
      }

      if (status!=INVALID_NUMBER_DETECTED) {
        // Create a SolveStatistics object
        statistics_ = new SolveStatistics(p2ip_nlp, p2ip_data, p2ip_cq);
      }
    }
    catch (TOO_FEW_DOF& exc) {
      exc.ReportException(*jnlst_, J_STRONGWARNING);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Problem has too few degrees of freedom.\n");
      retValue = Not_Enough_Degrees_Of_Freedom;
      status = TOO_FEW_DEGREES_OF_FREEDOM;
    }
    catch (OPTION_INVALID& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid option encountered.\n");
      retValue = Invalid_Option;
      status = INVALID_OPTION;
    }
    catch (NO_FREE_VARIABLES_BUT_FEASIBLE& exc) {
      exc.ReportException(*jnlst_, J_MOREDETAILED);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Optimal Solution Found.\n");
      retValue = Solve_Succeeded;
      status = SUCCESS;
      skip_finalize_solution_call = true; /* has already been called by TNLPAdapter (and we don't know the correct primal solution) */
    }
    catch (NO_FREE_VARIABLES_AND_INFEASIBLE& exc) {
      exc.ReportException(*jnlst_, J_MOREDETAILED);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Problem has only fixed variables and constraints are infeasible.\n");
      retValue = Infeasible_Problem_Detected;
      status = LOCAL_INFEASIBILITY;
      skip_finalize_solution_call = true; /* has already been called by TNLPAdapter (and we don't know the correct primal solution) */
    }
    catch (INCONSISTENT_BOUNDS& exc) {
      exc.ReportException(*jnlst_, J_MOREDETAILED);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Problem has inconsistent variable bounds or constraint sides.\n");
      retValue = Infeasible_Problem_Detected;
      status = LOCAL_INFEASIBILITY;
      skip_finalize_solution_call = true; /* has already been called by TNLPAdapter (and we don't know the correct primal solution) */
    }
    catch (IpoptException& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Some uncaught Ipopt exception encountered.\n");
      retValue = Unrecoverable_Exception;
    }
    catch (std::bad_alloc&) {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
      status = OUT_OF_MEMORY;
    }
    catch (...) {
      if( !rethrow_nonipoptexception_ )
      {
        IpoptException exc("Unknown Exception caught in Ipopt", "Unknown File", -1);
        exc.ReportException(*jnlst_, J_ERROR);
        retValue = NonIpopt_Exception_Thrown;
      }
      else
      {
        jnlst_->FlushBuffer();
        throw;
      }
    }

    if (!skip_finalize_solution_call)
      options_->GetBoolValue("skip_finalize_solution_call", skip_finalize_solution_call, "");

    if (!skip_finalize_solution_call && IsValid(p2ip_data->curr()) && IsValid(p2ip_data->curr()->x())) {
      SmartPtr<const Vector> c;
      SmartPtr<const Vector> d;
      SmartPtr<const Vector> zL;
      SmartPtr<const Vector> zU;
      SmartPtr<const Vector> yc;
      SmartPtr<const Vector> yd;
      Number obj = 0.;

      switch (status) {
      case SUCCESS:
      case MAXITER_EXCEEDED:
      case STOP_AT_TINY_STEP:
      case STOP_AT_ACCEPTABLE_POINT:
      case LOCAL_INFEASIBILITY:
      case USER_REQUESTED_STOP:
      case FEASIBLE_POINT_FOUND:
      case DIVERGING_ITERATES:
      case RESTORATION_FAILURE:
      case ERROR_IN_STEP_COMPUTATION:
        c = p2ip_cq->curr_c();
        d = p2ip_cq->curr_d();
        obj = p2ip_cq->curr_f();
        zL = p2ip_data->curr()->z_L();
        zU = p2ip_data->curr()->z_U();
        yc = p2ip_data->curr()->y_c();
        yd = p2ip_data->curr()->y_d();
        break;
      default: {
          SmartPtr<Vector> tmp = p2ip_data->curr()->y_c()->MakeNew();
          tmp->Set(0.);
          c = ConstPtr(tmp);
          yc = ConstPtr(tmp);
          tmp = p2ip_data->curr()->y_d()->MakeNew();
          tmp->Set(0.);
          d = ConstPtr(tmp);
          yd = ConstPtr(tmp);
          tmp = p2ip_data->curr()->z_L()->MakeNew();
          tmp->Set(0.);
          zL = ConstPtr(tmp);
          tmp = p2ip_data->curr()->z_U()->MakeNew();
          tmp->Set(0.);
          zU = ConstPtr(tmp);
        }
      }

      p2ip_nlp->FinalizeSolution(status,
                                 *p2ip_data->curr()->x(),
                                 *zL, *zU, *c, *d, *yc, *yd,
                                 obj, p2ip_data, p2ip_cq);
    }

    jnlst_->FlushBuffer();

    return retValue;
  }

  bool IpoptApplication::OpenOutputFile(std::string file_name,
                                        EJournalLevel print_level)
  {
    SmartPtr<Journal> file_jrnl = jnlst_->GetJournal("OutputFile:"+file_name);

    if (IsNull(file_jrnl)) {
      file_jrnl = jnlst_->AddFileJournal("OutputFile:"+file_name,
                                         file_name.c_str(),
                                         print_level);
    }

    // Check, if the output file could be created properly
    if (IsNull(file_jrnl)) {
      return false;
    }

    file_jrnl->SetPrintLevel(J_DBG, J_NONE);

    return true;
  }

  void IpoptApplication::RegisterAllIpoptOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    RegisterOptions_Interfaces(roptions);
    RegisterOptions_Algorithm(roptions);
    RegisterOptions_CGPenalty(roptions);
    RegisterOptions_LinearSolvers(roptions);
#ifdef BUILD_INEXACT
    RegisterOptions_Inexact(roptions);
#endif
  }

  SmartPtr<SolveStatistics> IpoptApplication::Statistics()
  {
    return statistics_;
  }

  SmartPtr<IpoptNLP> IpoptApplication::IpoptNLPObject()
  {
    return ip_nlp_;
  }

  SmartPtr<IpoptData> IpoptApplication::IpoptDataObject()
  {
    return ip_data_;
  }

  SmartPtr<IpoptCalculatedQuantities> IpoptApplication::IpoptCQObject()
  {
    return ip_cq_;
  }

  SmartPtr<IpoptAlgorithm> IpoptApplication::AlgorithmObject()
  {
    return alg_;
  }

  void IpoptApplication::PrintCopyrightMessage()
  {
    IpoptAlgorithm::print_copyright_message(*jnlst_);
  }

} // namespace Ipopt
