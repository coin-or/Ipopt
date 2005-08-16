// Copyright (C) 2004, 2005, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02

#include "IpIpoptApplication.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpIpoptAlg.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpAlgBuilder.hpp"
#include "IpUserScaling.hpp"
#include "IpGradientScaling.hpp"

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  DefineOptionsRegistrar(IpoptApplication);

  IpoptApplication::IpoptApplication(bool read_params_dat,
                                     bool create_console_out)
      :
      jnlst_(new Journalist()),
      options_(new OptionsList()),
      statistics_(NULL)
  {
# ifdef IP_DEBUG
    DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst_));
    SmartPtr<Journal> debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_SUMMARY);
    debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
# endif

    DBG_START_METH("IpoptApplication::IpoptApplication()",
                   dbg_verbosity);


    SmartPtr<Journal> stdout_jrnl = NULL;
    if (create_console_out) {
      stdout_jrnl =
        jnlst_->AddFileJournal("console", "stdout", J_SUMMARY);
      stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
    }

    // Register the valid options
    SmartPtr<RegisteredOptions> reg_options = new RegisteredOptions();
    OptionsRegistrar::RegisterAllOptions(reg_options);

    options_->SetJournalist(jnlst_);
    options_->SetRegisteredOptions(reg_options);

    // Get the options
    if (read_params_dat) {
      FILE* fp_options = fopen("PARAMS.DAT", "r");
      if (fp_options) {
        // PARAMS.DAT exists, read the content
        options_->ReadFromFile(*jnlst_, fp_options);
        fclose(fp_options);
        fp_options=NULL;
      }
    }

    Index ivalue;
    options_->GetIntegerValue("print_level", ivalue, "");
    EJournalLevel print_level = (EJournalLevel)ivalue;
    if (create_console_out) {
      // Set printlevel for stdout
      stdout_jrnl->SetAllPrintLevels(print_level);
      stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
    }

    bool option_set;

#ifdef IP_DEBUG
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
    debug_jrnl->SetAllPrintLevels(debug_print_level);
    debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
#endif

    // Open an output file if required
    std::string output_filename;
    options_->GetValue("output_file", output_filename, "");
    if (output_filename != "") {
      EJournalLevel file_print_level;
      option_set = options_->GetIntegerValue("file_print_level", ivalue, "");
      if (option_set) {
        file_print_level = (EJournalLevel)ivalue;
      }
      else {
        file_print_level = print_level;
      }
      OpenOutputFile(output_filename, file_print_level);
    }

    // output a description of all the options
    bool print_options_documentation;
    options_->GetBoolValue("print_options_documentation",
                           print_options_documentation, "");
    if (print_options_documentation) {
      reg_options->OutputOptionDocumentation(*jnlst_);
    }
  }

  IpoptApplication::~IpoptApplication()
  {
    DBG_START_METH("IpoptApplication::~IpoptApplication()",
                   dbg_verbosity);
  }

  void IpoptApplication::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedIntegerOption(
      "print_level",
      "Output verbosity level.",
      0, J_LAST_LEVEL-1, J_SUMMARY,
      "Sets the default verbosity level, in particular to screen.  The "
      "larger this value the more detailed is the output.");

#if IP_DEBUG

    roptions->AddBoundedIntegerOption(
      "debug_print_level",
      "Verbosity level for debug file.",
      0, J_LAST_LEVEL-1, J_SUMMARY,
      "This Ipopt library has been compiled in debug mode, and a file "
      "\"debug.out\" is produced for every run.  This option determines "
      "the verbosity level for this file.  By default it is the same as "
      "\"print_level\".");
#endif

    roptions->AddStringOption1(
      "output_file",
      "File name of an output file (leave unset for no file output)",
      "",
      "*", "Any acceptable standard file name",
      "NOTE: This option only works when read from PARAMS.DAT options file! "
      "An output file with this name will be written (leave unset for no "
      "file output).  The verbosity level is either \"print_level\", or "
      "\"file_print_level\" if that is given.");
    roptions->AddBoundedIntegerOption(
      "file_print_level",
      "Verbosity level output file.",
      0, J_LAST_LEVEL-1, J_SUMMARY,
      "NOTE: This option only works when read from PARAMS.DAT options file! "
      "Determines the verbosity level for the file specified by "
      "\"output_file\".  By defauly it is the same as \"print_level\".");
    roptions->AddStringOption2(
      "print_options_documentation",
      "Switch to print list all algorithmic options",
      "no",
      "no", "don't print list",
      "yes", "print list",
      "If selected, the algorithm will print the list of all available "
      "algorithmic options with some documentation before solving the "
      "optimization problem.");

    roptions->AddStringOption3(
      "nlp_scaling_method",
      "Select the technique used for scaling the NLP", "gradient_based",
      "none", "no problem scaling will be performed",
      "user_scaling", "scaling parameters will come from the user",
      "gradient_based", "scale the problem so the maximum gradient at the starting point is scaling_max_gradient",
      "Selects the technique used for scaling the problem before it is solved."
      " For user-scaling, the parameters come from the NLP. If you are using "
      "AMPL, they can be specified through suffixes (scaling_factor)");
  }

  ApplicationReturnStatus IpoptApplication::OptimizeTNLP(const SmartPtr<TNLP>& nlp)
  {
    SmartPtr<NLP> nlp_adapter =
      new TNLPAdapter(GetRawPtr(nlp));

    return OptimizeNLP(nlp_adapter);
  }

  ApplicationReturnStatus
  IpoptApplication::OptimizeNLP(const SmartPtr<NLP>& nlp)
  {
    statistics_ = NULL; /* delete old statistics */
    ApplicationReturnStatus retValue = Internal_Error;
    SmartPtr<IpoptData> ip_data;
    SmartPtr<IpoptCalculatedQuantities> ip_cq;
    SmartPtr<IpoptNLP> ip_nlp;
    try {

      SmartPtr<NLPScalingObject> nlp_scaling;
      std::string nlp_scaling_method;
      options_->GetValue("nlp_scaling_method", nlp_scaling_method, "");
      if (nlp_scaling_method == "user_scaling") {
        nlp_scaling = new UserScaling(ConstPtr(nlp));
      }
      else if (nlp_scaling_method == "gradient_based") {
        nlp_scaling = new GradientScaling(nlp);
      }
      else {
        nlp_scaling = new NoNLPScalingObject();
      }

      ip_nlp =  new OrigIpoptNLP(ConstPtr(jnlst_), GetRawPtr(nlp), nlp_scaling);

      // Create the IpoptData
      if (IsNull(ip_data)) {
        ip_data = new IpoptData();
      }

      // Create the IpoptCalculators
      if (IsNull(ip_cq)) {
        ip_cq = new IpoptCalculatedQuantities(ip_nlp, ip_data);
      }

      // Create the Algorithm object
      SmartPtr<IpoptAlgorithm> alg
      = AlgorithmBuilder::BuildBasicAlgorithm(*jnlst_, *options_, "");

      // Set up the algorithm
      alg->Initialize(*jnlst_, *ip_nlp, *ip_data, *ip_cq, *options_, "");

      if( jnlst_->ProduceOutput(J_DETAILED, J_MAIN) ) {
        // Print out the options (including the number of times they were used
        std::string liststr;
        options_->PrintList(liststr);
        jnlst_->Printf(J_DETAILED, J_MAIN, "\nList of options:\n\n%s", liststr.c_str());
      }

      // Run the algorithm
      SolverReturn status = alg->Optimize();

      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "\nNumber of Iterations....: %d\n",
                     ip_data->iter_count());

      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "\n                                   (scaled)                 (unscaled)\n");
      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "Objective...............: %24.16e  %24.16e\n", ip_cq->curr_f(), ip_cq->unscaled_curr_f());
      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "Dual infeasibility......: %24.16e  %24.16e\n", ip_cq->curr_dual_infeasibility(NORM_MAX), ip_cq->unscaled_curr_dual_infeasibility(NORM_MAX));
      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "Constraint violation....: %24.16e  %24.16e\n", ip_cq->curr_nlp_constraint_violation(NORM_MAX), ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX));
      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "Complementarity.........: %24.16e  %24.16e\n", ip_cq->curr_complementarity(0., NORM_MAX), ip_cq->unscaled_curr_complementarity(0., NORM_MAX));
      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "Overall NLP error.......: %24.16e  %24.16e\n\n", ip_cq->curr_nlp_error(), ip_cq->unscaled_curr_nlp_error());

      ip_data->curr()->x()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "x");
      ip_data->curr()->y_c()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_c");
      ip_data->curr()->y_d()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_d");
      ip_data->curr()->z_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_L");
      ip_data->curr()->z_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_U");
      ip_data->curr()->v_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_L");
      ip_data->curr()->v_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_U");

      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "\nNumber of objective function evaluations             = %d\n",
                     ip_nlp->f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of objective gradient evaluations             = %d\n",
                     ip_nlp->grad_f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of equality constraint evaluations            = %d\n",
                     ip_nlp->c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of inequality constraint evaluations          = %d\n",
                     ip_nlp->d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of equality constraint Jacobian evaluations   = %d\n",
                     ip_nlp->jac_c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of inequality constraint Jacobian evaluations = %d\n",
                     ip_nlp->jac_d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of Lagrangian Hessian evaluations             = %d\n",
                     ip_nlp->h_evals());

      // Write EXIT message
      if (status == SUCCESS) {
        retValue = Solve_Succeeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Optimal Solution Found.\n");
      }
      else if (status == MAXITER_EXCEEDED) {
        retValue = Maximum_Iterations_Exceeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum Number of Iterations Exceeded.\n");
      }
      else if (status == STOP_AT_TINY_STEP) {
        retValue = Solved_To_Best_Possible_Precision;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Solved To Best Possible Precision.\n");
      }
      else if (status == STOP_AT_ACCEPTABLE_POINT) {
        retValue = Solved_To_Acceptable_Level;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Solved To Acceptable Level.\n");
      }
      else if (status == RESTORATION_FAILURE) {
        retValue = Restoration_Failed;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Restoration Failed!\n");
      }
      else if (status == LOCAL_INFEASIBILITY) {
        retValue = Infeasible_Problem_Detected;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Converged to a point of local infeasibility. Problem may be infeasible.\n");
      }
      else {
        retValue = Internal_Error;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
        return retValue;
      }
      ip_nlp->FinalizeSolution(status,
                               *ip_data->curr()->x(), *ip_data->curr()->z_L(), *ip_data->curr()->z_U(),
                               *ip_cq->curr_c(), *ip_cq->curr_d(), *ip_data->curr()->y_c(), *ip_data->curr()->y_d(),
                               ip_cq->curr_f());
      // Create a SolveStatistics object
      statistics_ = new SolveStatistics(ip_nlp, ip_data, ip_cq);
    }
    catch(TOO_FEW_DOF& exc) {
      exc.ReportException(*jnlst_);
      retValue = Not_Enough_Degrees_Of_Freedom;
    }
    catch(OPTION_INVALID& exc) {
      exc.ReportException(*jnlst_);
      retValue = Invalid_Option;
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst_);
      retValue = Unrecoverable_Exception;
    }
    catch(std::bad_alloc& exc) {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
    }
    catch(...) {
      IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
      exc.ReportException(*jnlst_);
      retValue = NonIpopt_Exception_Thrown;
    }

    jnlst_->FlushBuffer();

    return retValue;

  }

  bool IpoptApplication::OpenOutputFile(std::string file_name,
                                        EJournalLevel print_level)
  {
    SmartPtr<Journal> file_jrnl = jnlst_->AddFileJournal("OutputFile:"+file_name,
                                  file_name.c_str(),
                                  print_level);
    file_jrnl->SetPrintLevel(J_DBG, J_NONE);

    return true;
  }

} // namespace Ipopt
