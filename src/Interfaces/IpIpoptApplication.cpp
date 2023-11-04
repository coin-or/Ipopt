// Copyright (C) 2004, 2012 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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
#include "IpLibraryLoader.hpp"
#include "IpLinearSolvers.h"

#ifdef BUILD_INEXACT
# include "IpInexactRegOp.hpp"
# include "IpInexactAlgBuilder.hpp"
#endif

#include <cassert>
#include <cmath>
#include <fstream>

// Factory to facilitate creating IpoptApplication objects from within a DLL

Ipopt::IpoptApplication* IpoptApplicationFactory()
{
   return new Ipopt::IpoptApplication;
}

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
// Kluge: Add reference counter for DebugJournalistWrapper::jrnl
static SmartPtr<Journalist> smart_jnlst(NULL);
#endif

IpoptApplication::IpoptApplication(
   bool create_console_out /* = true */,
   bool create_empty /* = false */
)
   : read_params_dat_(true),
     rethrow_nonipoptexception_(false),
     options_(new OptionsList()),
     inexact_algorithm_(false),
     replace_bounds_(false)
{
   if( create_empty )
   {
      return;
   }

   jnlst_ = new Journalist();
   try
   {
# if IPOPT_VERBOSITY > 0
      // Kludge: If this is the first IpoptApplication, then store jnlst_ in smart_jnlst, too, so that it doesn't
      // get freed when the IpoptApplication is freed and DebugJournalistWrapper::jrnl becomes a dangling pointer.
      // Also add the Debug journal that writes to debug.out.
      if( IsNull(smart_jnlst) )
      {
         smart_jnlst = jnlst_;
         DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst_));
         SmartPtr<Journal> debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_ITERSUMMARY);
         debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
      }
# endif

      DBG_START_METH("IpoptApplication::IpoptApplication()",
                     dbg_verbosity);

      if( create_console_out )
      {
         SmartPtr<Journal> stdout_jrnl = jnlst_->AddFileJournal("console", "stdout", J_ITERSUMMARY);
         stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
      }

      // Register the valid options
      reg_options_ = new RegisteredOptions();
      RegisterAllIpoptOptions(reg_options_);

      options_->SetJournalist(jnlst_);
      options_->SetRegisteredOptions(reg_options_);
   }
   catch( IpoptException& exc )
   {
      exc.ReportException(*jnlst_);
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR, "Caught unknown Ipopt exception");
   }
   catch( std::bad_alloc& )
   {
      jnlst_->Printf(J_ERROR, J_MAIN, "\nEXIT: Not enough memory.\n");
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR, "Not enough memory");
   }
   catch( std::overflow_error& )
   {
      jnlst_->Printf(J_ERROR, J_MAIN, "\nEXIT: Integer type too small for required memory.\n");
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR, "Not enough memory");
   }
   catch( ... )
   {
      IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
      exc.ReportException(*jnlst_);
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR, "Caught unknown exception");
   }
}

IpoptApplication::IpoptApplication(
   SmartPtr<RegisteredOptions> reg_options,
   SmartPtr<OptionsList>       options,
   SmartPtr<Journalist>        jnlst
)
   : read_params_dat_(true),
     rethrow_nonipoptexception_(false),
     jnlst_(jnlst),
     reg_options_(reg_options),
     options_(options),
     inexact_algorithm_(false),
     replace_bounds_(false)
{
#if IPOPT_VERBOSITY > 0
   // Kludge: If this is the first IpoptApplication, then store jnlst_ in smart_jnlst, too, so that it doesn't
   // get freed when the IpoptApplication is freed and DebugJournalistWrapper::jrnl becomes a dangling pointer.
   // Also add the Debug journal that writes to debug.out.
   if( IsNull(smart_jnlst) )
   {
      smart_jnlst = jnlst_;
      DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst_));
      SmartPtr<Journal> debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_ITERSUMMARY);
      debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
   }
#endif
}

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

ApplicationReturnStatus IpoptApplication::Initialize(
   std::istream& is,
   bool          allow_clobber
)
{
   try
   {
      // Get the options
      if( is.good() )
      {
         // stream exists, read the content
         options_->ReadFromStream(*jnlst_, is, allow_clobber);
      }

      bool no_output;
      options_->GetBoolValue("suppress_all_output", no_output, "");

      if( no_output )
      {
         jnlst_->DeleteAllJournals();
      }
      else
      {
         Index ivalue;
         options_->GetIntegerValue("print_level", ivalue, "");
         EJournalLevel print_level = (EJournalLevel) ivalue;
         SmartPtr<Journal> stdout_jrnl = jnlst_->GetJournal("console");
         if( IsValid(stdout_jrnl) )
         {
            // Set printlevel for stdout
            stdout_jrnl->SetAllPrintLevels(print_level);
            stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
         }

#if IPOPT_VERBOSITY > 0
         // Set printlevel for debug
         EJournalLevel debug_print_level;
         if( options_->GetIntegerValue("debug_print_level", ivalue, "") )
         {
            debug_print_level = (EJournalLevel)ivalue;
         }
         else
         {
            debug_print_level = print_level;
         }
         assert(IsValid(smart_jnlst));  /* should have been created in constructor */
         SmartPtr<Journal> debug_jrnl = smart_jnlst->GetJournal("Debug");
         assert(IsValid(debug_jrnl));  /* should have been added in constructor */
         debug_jrnl->SetAllPrintLevels(debug_print_level);
         debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
         if( IsNull(jnlst_->GetJournal("Debug")) )
         {
            jnlst_->AddJournal(debug_jrnl);
         }
#endif

         // Open an output file if required
         std::string output_filename;
         options_->GetStringValue("output_file", output_filename, "");
         if( output_filename != "" )
         {
            EJournalLevel file_print_level;
            if( options_->GetIntegerValue("file_print_level", ivalue, "") )
            {
               file_print_level = (EJournalLevel) ivalue;
            }
            else
            {
               file_print_level = print_level;
            }
            bool file_append;
            options_->GetBoolValue("file_append", file_append, "");
            bool openend = OpenOutputFile(output_filename, file_print_level, file_append);
            if( !openend )
            {
               jnlst_->Printf(J_ERROR, J_INITIALIZATION, "Error opening output file \"%s\"\n", output_filename.c_str());
               return Invalid_Option;
            }
         }
      }

      // output a description of all the options
      bool print_options_documentation;
      options_->GetBoolValue("print_options_documentation", print_options_documentation, "");
      if( print_options_documentation )
      {
         reg_options_->OutputOptionDocumentation(*jnlst_, options_);
      }

#ifdef BUILD_INEXACT
      // Check if we are to use the inexact linear solver option
      options_->GetBoolValue("inexact_algorithm", inexact_algorithm_, "");
      // Change the default flags for the inexact algorithm
      if (inexact_algorithm_)
      {
         AddInexactDefaultOptions(*options_);
      }
#endif

      options_->GetBoolValue("replace_bounds", replace_bounds_, "");
   }
   catch( OPTION_INVALID& exc )
   {
      exc.ReportException(*jnlst_, J_ERROR);
      return Invalid_Option;
   }
   catch( IpoptException& exc )
   {
      exc.ReportException(*jnlst_, J_ERROR);
      return Unrecoverable_Exception;
   }
   catch( std::bad_alloc& )
   {
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
      return Insufficient_Memory;
   }
   catch( std::overflow_error& )
   {
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Integer type too small for required memory.\n");
      return Insufficient_Memory;
   }
   catch( ... )
   {
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

ApplicationReturnStatus IpoptApplication::Initialize(
   std::string params_file,
   bool        allow_clobber
)
{
   std::ifstream is;
   if( params_file != "" )
   {
      try
      {
         is.open(params_file.c_str());
      }
      catch( std::bad_alloc& )
      {
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
         return Insufficient_Memory;
      }
      catch( std::overflow_error& )
      {
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Integer type too small for required memory.\n");
         return Insufficient_Memory;
      }
      catch( ... )
      {
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
   if( is )
   {
      is.close();
   }
   return retval;
}

ApplicationReturnStatus IpoptApplication::Initialize(
   bool allow_clobber
)
{
   std::string option_file_name;
   options_->GetStringValue("option_file_name", option_file_name, "");
   if( option_file_name != "" && option_file_name != "ipopt.opt" )
   {
      jnlst_->Printf(J_SUMMARY, J_MAIN, "Using option file \"%s\".\n\n", option_file_name.c_str());
   }

   return Initialize(option_file_name, allow_clobber);
}

IpoptApplication::~IpoptApplication()
{
   DBG_START_METH("IpoptApplication::~IpoptApplication()",
                  dbg_verbosity);
}

void IpoptApplication::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->SetRegisteringCategory("Output");
   roptions->AddBoundedIntegerOption(
      "print_level",
      "Output verbosity level.",
      0, J_LAST_LEVEL - 1,
      J_ITERSUMMARY,
      "Sets the default verbosity level for console output. "
      "The larger this value the more detailed is the output.");

   roptions->AddStringOption1(
      "output_file",
      "File name of desired output file (leave unset for no file output).",
      "",
      "*", "Any acceptable standard file name",
      "NOTE: This option only works when read from the ipopt.opt options file! "
      "An output file with this name will be written (leave unset for no file output). "
      "The verbosity level is by default set to \"print_level\", but can be overridden with \"file_print_level\". "
      "The file name is changed to use only small letters.");
   roptions->AddBoundedIntegerOption(
      "file_print_level",
      "Verbosity level for output file.",
      0, J_LAST_LEVEL - 1,
      J_ITERSUMMARY,
      "NOTE: This option only works when read from the ipopt.opt options file! "
      "Determines the verbosity level for the file specified by \"output_file\". "
      "By default it is the same as \"print_level\".");
   roptions->AddBoolOption(
      "file_append",
      "Whether to append to output file, if set, instead of truncating.",
      false,
      "NOTE: This option only works when read from the ipopt.opt options file!");
   roptions->AddBoolOption(
      "print_user_options",
      "Print all options set by the user.",
      false,
      "If selected, the algorithm will print the list of all options set by the user including their values and whether they have been used. "
      "In some cases this information might be incorrect, due to the internal program flow.");
   roptions->AddBoolOption(
      "print_options_documentation",
      "Switch to print all algorithmic options with some documentation before solving the optimization problem.",
      false);

#if IPOPT_VERBOSITY > 0
   roptions->AddBoundedIntegerOption(
      "debug_print_level",
      "Verbosity level for debug file.",
      0, J_LAST_LEVEL - 1,
      J_ITERSUMMARY,
      "This Ipopt library has been compiled in debug mode, and a file \"debug.out\" is produced for every run. "
      "This option determines the verbosity level for this file. "
      "By default it is the same as \"print_level\".");
#endif

   roptions->AddBoolOption(
      "print_timing_statistics",
      "Switch to print timing statistics.",
      false,
      "If selected, the program will print the time spend for selected tasks. "
      "This implies timing_statistics=yes.");

   roptions->SetRegisteringCategory("Miscellaneous");
   roptions->AddStringOption1(
      "option_file_name",
      "File name of options file.",
      "ipopt.opt",
      "*", "Any acceptable standard file name",
      "By default, the name of the Ipopt options file is \"ipopt.opt\" - "
      "or something else if specified in the IpoptApplication::Initialize call. "
      "If this option is set by SetStringValue BEFORE the options file is read, it specifies the name of the options file. "
      "It does not make any sense to specify this option within the options file. "
      "Setting this option to an empty string disables reading of an options file.");

   roptions->AddBoolOption(
      "replace_bounds",
      "Whether all variable bounds should be replaced by inequality constraints",
      false,
      "This option must be set for the inexact algorithm.",
      true);
   roptions->AddBoolOption(
      "skip_finalize_solution_call",
      "Whether a call to NLP::FinalizeSolution after optimization should be suppressed",
      false,
      "In some Ipopt applications, the user might want to call the FinalizeSolution method separately. "
      "Setting this option to \"yes\" will cause the IpoptApplication object to suppress the default call to that method.",
      true);

   roptions->SetRegisteringCategory("Undocumented");
   roptions->AddBoolOption(
      "suppress_all_output",
      "",
      false,
      "",
      true);
#ifdef BUILD_INEXACT
   roptions->AddBoolOption(
      "inexact_algorithm",
      "Whether to activate the version of Ipopt that allows iterative linear solvers.",
      false,
      "EXPERIMENTAL",
      true);
#endif
}

ApplicationReturnStatus IpoptApplication::OptimizeTNLP(
   const SmartPtr<TNLP>& tnlp
)
{
   nlp_adapter_ = new TNLPAdapter(GetRawPtr(tnlp), ConstPtr(jnlst_));
   return OptimizeNLP(nlp_adapter_);
}

ApplicationReturnStatus IpoptApplication::ReOptimizeTNLP(
   const SmartPtr<TNLP>& tnlp
)
{
   ASSERT_EXCEPTION(IsValid(nlp_adapter_), INVALID_WARMSTART, "ReOptimizeTNLP called before OptimizeTNLP.");
   TNLPAdapter* adapter = static_cast<TNLPAdapter*>(GetRawPtr(nlp_adapter_));
   DBG_ASSERT(dynamic_cast<TNLPAdapter*> (GetRawPtr(nlp_adapter_)));
   ASSERT_EXCEPTION(adapter->tnlp() == tnlp, INVALID_WARMSTART, "ReOptimizeTNLP called for different TNLP.")

   return ReOptimizeNLP(nlp_adapter_);
}

ApplicationReturnStatus IpoptApplication::OptimizeNLP(
   const SmartPtr<NLP>& nlp
)
{
   SmartPtr<AlgorithmBuilder> alg_builder = NULL;
   return OptimizeNLP(nlp, alg_builder);
}

ApplicationReturnStatus IpoptApplication::OptimizeNLP(
   const SmartPtr<NLP>&        nlp,
   SmartPtr<AlgorithmBuilder>& alg_builder
)
{
   ApplicationReturnStatus retValue = Internal_Error;

   // Prepare internal data structures of the algorithm
   try
   {

      if( IsNull(alg_builder) )
      {
#ifdef BUILD_INEXACT
         if (inexact_algorithm_)
         {
            alg_builder = new InexactAlgorithmBuilder();
         }
         else
         {
#endif
            alg_builder = new AlgorithmBuilder();
#ifdef BUILD_INEXACT
         }
#endif
      }

      SmartPtr<NLP> use_nlp;
      if( replace_bounds_ )
      {
         use_nlp = new NLPBoundsRemover(*nlp);
      }
      else
      {
         use_nlp = nlp;
      }
      alg_builder->BuildIpoptObjects(*jnlst_, *options_, "", use_nlp, ip_nlp_, ip_data_, ip_cq_);

      alg_ = GetRawPtr(alg_builder->BuildBasicAlgorithm(*jnlst_, *options_, ""));

      // finally call the optimization
      retValue = call_optimize();
   }
   catch( OPTION_INVALID& exc )
   {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid option encountered.\n");
      retValue = Invalid_Option;
   }
   catch( IpoptException& exc )
   {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Some uncaught Ipopt exception encountered.\n");
      retValue = Unrecoverable_Exception;
   }
   catch( std::bad_alloc& )
   {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
   }
   catch( std::overflow_error& )
   {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Integer type too small for required memory.\n");
   }
   catch( ... )
   {
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

ApplicationReturnStatus IpoptApplication::ReOptimizeNLP(
   const SmartPtr<NLP>& nlp
)
{
   ASSERT_EXCEPTION(IsValid(alg_), INVALID_WARMSTART, "ReOptimizeNLP called before OptimizeNLP.");
   OrigIpoptNLP* orig_nlp = static_cast<OrigIpoptNLP*>(GetRawPtr(ip_nlp_));
   DBG_ASSERT(dynamic_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_)));
   ASSERT_EXCEPTION(orig_nlp->nlp() == nlp, INVALID_WARMSTART, "ReOptimizeTNLP called for different NLP.")

   return call_optimize();
}

ApplicationReturnStatus IpoptApplication::call_optimize()
{
   // Reset the print-level for the screen output
   Index ivalue;
   options_->GetIntegerValue("print_level", ivalue, "");
   EJournalLevel print_level = (EJournalLevel) ivalue;
   SmartPtr<Journal> stdout_jrnl = jnlst_->GetJournal("console");
   if( IsValid(stdout_jrnl) )
   {
      // Set printlevel for stdout
      stdout_jrnl->SetAllPrintLevels(print_level);
      stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
   }

   statistics_ = NULL; /* delete old statistics */
   // Get the pointers to the real objects (need to do it that
   // awkwardly since otherwise we would have to include so many
   // things in IpoptApplication, which a user would have to
   // include, too (SmartPtr doesn't work otherwise on AIX)
   IpoptAlgorithm* p2alg = static_cast<IpoptAlgorithm*>(GetRawPtr(alg_));
   DBG_ASSERT(dynamic_cast<IpoptAlgorithm*> (GetRawPtr(alg_)));
   IpoptData* p2ip_data = static_cast<IpoptData*>(GetRawPtr(ip_data_));
   DBG_ASSERT(dynamic_cast<IpoptData*> (GetRawPtr(ip_data_)));
   OrigIpoptNLP* p2ip_nlp = static_cast<OrigIpoptNLP*>(GetRawPtr(ip_nlp_));
   DBG_ASSERT(dynamic_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_)));
   IpoptCalculatedQuantities* p2ip_cq = static_cast<IpoptCalculatedQuantities*>(GetRawPtr(ip_cq_));
   DBG_ASSERT(dynamic_cast<IpoptCalculatedQuantities*> (GetRawPtr(ip_cq_)));

   // Reset Timing statistics
   ip_data_->TimingStats().ResetTimes();

   ApplicationReturnStatus retValue = Internal_Error;
   SolverReturn status = INTERNAL_ERROR;
   try
   {
      // check whether timing statistics need to be printed
      bool print_timing_statistics;
      options_->GetBoolValue("print_timing_statistics", print_timing_statistics, "");
      // enable collecting timing statistics if they need to be printed later
      if( print_timing_statistics )
      {
         options_->SetStringValue("timing_statistics", "yes", true, true);
      }

      // Set up the algorithm
      p2alg->Initialize(*jnlst_, *p2ip_nlp, *p2ip_data, *p2ip_cq, *options_, "");

      // If selected, print the user options
      bool print_user_options;
      options_->GetBoolValue("print_user_options", print_user_options, "");
      if( print_user_options )
      {
         std::string liststr;
         options_->PrintUserOptions(liststr);
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nList of user-set options:\n\n%s", liststr.c_str());
      }

      if( jnlst_->ProduceOutput(J_DETAILED, J_MAIN) )
      {
         // Print out the options (including the number of times they were used
         std::string liststr;
         options_->PrintList(liststr);
         jnlst_->Printf(J_DETAILED, J_MAIN, "\nList of options:\n\n%s", liststr.c_str());
      }

      // Run the algorithm
      status = p2alg->Optimize();

      // Since all the output below doesn't make any sense in this
      // case, we rethrow the TOO_FEW_DOF exception here
      ASSERT_EXCEPTION(status != TOO_FEW_DEGREES_OF_FREEDOM, TOO_FEW_DOF, "Too few degrees of freedom (rethrown)!");

      jnlst_->Printf(J_SUMMARY, J_SOLUTION, "\nNumber of Iterations....: %" IPOPT_INDEX_FORMAT "\n", p2ip_data->iter_count());

      if( status != INVALID_NUMBER_DETECTED )
      {
         try
         {
            jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                           "\n                                   (scaled)                 (unscaled)\n");
            jnlst_->Printf(J_SUMMARY, J_SOLUTION, "Objective...............: %24.16e  %24.16e\n", p2ip_cq->curr_f(),
                           p2ip_cq->unscaled_curr_f());
            jnlst_->Printf(J_SUMMARY, J_SOLUTION, "Dual infeasibility......: %24.16e  %24.16e\n",
                           p2ip_cq->curr_dual_infeasibility(NORM_MAX), p2ip_cq->unscaled_curr_dual_infeasibility(NORM_MAX));
            jnlst_->Printf(J_SUMMARY, J_SOLUTION, "Constraint violation....: %24.16e  %24.16e\n",
                           p2ip_cq->curr_nlp_constraint_violation(NORM_MAX),
                           p2ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX));
            jnlst_->Printf(J_SUMMARY, J_SOLUTION, "Variable bound violation: %24.16e  %24.16e\n",
                           p2ip_cq->curr_orig_bounds_violation(NORM_MAX),
                           p2ip_cq->unscaled_curr_orig_bounds_violation(NORM_MAX));
            jnlst_->Printf(J_SUMMARY, J_SOLUTION, "Complementarity.........: %24.16e  %24.16e\n",
                           p2ip_cq->curr_complementarity(0., NORM_MAX), p2ip_cq->unscaled_curr_complementarity(0., NORM_MAX));
            jnlst_->Printf(J_SUMMARY, J_SOLUTION, "Overall NLP error.......: %24.16e  %24.16e\n\n",
                           p2ip_cq->curr_nlp_error(), p2ip_cq->unscaled_curr_nlp_error());
         }
         catch( IpoptNLP::Eval_Error& exc )
         {
            // this can happen if the final point was accepted because functions can be evaluated,
            // but functions are not differentiable, so dual infeasibility cannot be computed
            status = INVALID_NUMBER_DETECTED;
            exc.ReportException(*jnlst_, J_STRONGWARNING);
         }
      }

      p2ip_data->curr()->x()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "x");
      p2ip_data->curr()->y_c()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_c");
      p2ip_data->curr()->y_d()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_d");
      p2ip_data->curr()->z_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_L");
      p2ip_data->curr()->z_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_U");
      p2ip_data->curr()->v_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_L");
      p2ip_data->curr()->v_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_U");

      if( status == LOCAL_INFEASIBILITY )
      {
         p2ip_cq->curr_c()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "curr_c");
         p2ip_cq->curr_d_minus_s()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "curr_d_minus_s");
      }

      jnlst_->Printf(J_SUMMARY, J_STATISTICS, "\nNumber of objective function evaluations             = %" IPOPT_INDEX_FORMAT "\n",
                     p2ip_nlp->f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Number of objective gradient evaluations             = %" IPOPT_INDEX_FORMAT "\n",
                     p2ip_nlp->grad_f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Number of equality constraint evaluations            = %" IPOPT_INDEX_FORMAT "\n",
                     p2ip_nlp->c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Number of inequality constraint evaluations          = %" IPOPT_INDEX_FORMAT "\n",
                     p2ip_nlp->d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Number of equality constraint Jacobian evaluations   = %" IPOPT_INDEX_FORMAT "\n",
                     p2ip_nlp->jac_c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Number of inequality constraint Jacobian evaluations = %" IPOPT_INDEX_FORMAT "\n",
                     p2ip_nlp->jac_d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Number of Lagrangian Hessian evaluations             = %" IPOPT_INDEX_FORMAT "\n",
                     p2ip_nlp->h_evals());
      Number wall_time_overall_alg = p2ip_data->TimingStats().OverallAlgorithm().TotalWallclockTime();
      if( p2ip_data->TimingStats().IsFunctionEvaluationTimeEnabled() )
      {
         Number wall_time_funcs = p2ip_data->TimingStats().TotalFunctionEvaluationWallclockTime();
         jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Total seconds in IPOPT (w/o function evaluations)    = %10.3f\n",
                        wall_time_overall_alg - wall_time_funcs);
         jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Total seconds in NLP function evaluations            = %10.3f\n",
                        wall_time_funcs);
      }
      else
      {
         jnlst_->Printf(J_SUMMARY, J_STATISTICS, "Total seconds in IPOPT                               = %.3f\n",
                        wall_time_overall_alg);
      }

      // Write timing statistics information
      if( print_timing_statistics )
      {
         jnlst_->Printf(J_SUMMARY, J_TIMING_STATISTICS, "\n\nTiming Statistics:\n\n");
         p2ip_data->TimingStats().PrintAllTimingStatistics(*jnlst_, J_SUMMARY, J_TIMING_STATISTICS);
      }

      // Write EXIT message
      if( status == SUCCESS )
      {
         retValue = Solve_Succeeded;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Optimal Solution Found.\n");
      }
      else if( status == MAXITER_EXCEEDED )
      {
         retValue = Maximum_Iterations_Exceeded;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum Number of Iterations Exceeded.\n");
      }
      else if( status == CPUTIME_EXCEEDED )
      {
         retValue = Maximum_CpuTime_Exceeded;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum CPU time exceeded.\n");
      }
      else if( status == WALLTIME_EXCEEDED )
      {
         retValue = Maximum_WallTime_Exceeded;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum wallclock time exceeded.\n");
      }
      else if( status == STOP_AT_TINY_STEP )
      {
         retValue = Search_Direction_Becomes_Too_Small;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Search Direction is becoming Too Small.\n");
      }
      else if( status == STOP_AT_ACCEPTABLE_POINT )
      {
         retValue = Solved_To_Acceptable_Level;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Solved To Acceptable Level.\n");
      }
      else if( status == FEASIBLE_POINT_FOUND )
      {
         retValue = Feasible_Point_Found;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Feasible point for square problem found.\n");
      }
      else if( status == DIVERGING_ITERATES )
      {
         retValue = Diverging_Iterates;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Iterates diverging; problem might be unbounded.\n");
      }
      else if( status == RESTORATION_FAILURE )
      {
         retValue = Restoration_Failed;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Restoration Failed!\n");
      }
      else if( status == ERROR_IN_STEP_COMPUTATION )
      {
         retValue = Error_In_Step_Computation;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Error in step computation!\n");
      }
      else if( status == LOCAL_INFEASIBILITY )
      {
         retValue = Infeasible_Problem_Detected;
         jnlst_->Printf(J_SUMMARY, J_MAIN,
                        "\nEXIT: Converged to a point of local infeasibility. Problem may be infeasible.\n");
      }
      else if( status == USER_REQUESTED_STOP )
      {
         retValue = User_Requested_Stop;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Stopping optimization at current point as requested by user.\n");
      }
      else if( status == INVALID_NUMBER_DETECTED )
      {
         retValue = Invalid_Number_Detected;
         jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid number in NLP function or derivative detected.\n");
      }
      else
      {
         retValue = Internal_Error;
         jnlst_->Printf(J_SUMMARY, J_MAIN,
                        "\nEXIT: INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
         return retValue;
      }

      if( status != INVALID_NUMBER_DETECTED )
      {
         // Create a SolveStatistics object
         statistics_ = new SolveStatistics(p2ip_nlp, p2ip_data, p2ip_cq);
      }
   }
   catch( TOO_FEW_DOF& exc )
   {
      exc.ReportException(*jnlst_, J_STRONGWARNING);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Problem has too few degrees of freedom.\n");
      retValue = Not_Enough_Degrees_Of_Freedom;
      status = TOO_FEW_DEGREES_OF_FREEDOM;
   }
   catch( OPTION_INVALID& exc )
   {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid option encountered.\n");
      retValue = Invalid_Option;
      status = INVALID_OPTION;
   }
   catch( DYNAMIC_LIBRARY_FAILURE& exc )
   {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Library loading failure.\n");
      retValue = Invalid_Option;
   }
   catch( INCONSISTENT_BOUNDS& exc )
   {
      exc.ReportException(*jnlst_, J_MOREDETAILED);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Problem has inconsistent variable bounds or constraint sides.\n");
      retValue = Invalid_Problem_Definition;
      status = LOCAL_INFEASIBILITY;
   }
   catch( IpoptException& exc )
   {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Some uncaught Ipopt exception encountered.\n");
      retValue = Unrecoverable_Exception;
   }
   catch( std::bad_alloc& )
   {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
      status = OUT_OF_MEMORY;
   }
   catch( std::overflow_error& )
   {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Integer type too small for required memory.\n");
      status = OUT_OF_MEMORY;
   }
   catch( ... )
   {
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

   /** Flag indicating if the NLP:FinalizeSolution method should not
    *  be called after optimization. */
   bool skip_finalize_solution_call;
   options_->GetBoolValue("skip_finalize_solution_call", skip_finalize_solution_call, "");

   if( !skip_finalize_solution_call && IsValid(p2ip_data->curr()) && IsValid(p2ip_data->curr()->x()) )
   {
      SmartPtr<const Vector> c;
      SmartPtr<const Vector> d;
      SmartPtr<const Vector> zL;
      SmartPtr<const Vector> zU;
      SmartPtr<const Vector> yc;
      SmartPtr<const Vector> yd;
      Number obj = 0.;

      switch( status )
      {
         case SUCCESS:
         case MAXITER_EXCEEDED:
         case CPUTIME_EXCEEDED:
         case WALLTIME_EXCEEDED:
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
         default:
         {
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

      p2ip_nlp->FinalizeSolution(status, *p2ip_data->curr()->x(), *zL, *zU, *c, *d, *yc, *yd, obj, p2ip_data, p2ip_cq);
   }

   jnlst_->FlushBuffer();

   return retValue;
}

bool IpoptApplication::OpenOutputFile(
   std::string   file_name,
   EJournalLevel print_level,
   bool          file_append
)
{
   SmartPtr<Journal> file_jrnl = jnlst_->GetJournal("OutputFile:" + file_name);

   if( IsNull(file_jrnl) )
   {
      file_jrnl = jnlst_->AddFileJournal("OutputFile:" + file_name, file_name.c_str(), print_level, file_append);
   }

   // Check, if the output file could be created properly
   if( IsNull(file_jrnl) )
   {
      return false;
   }

   file_jrnl->SetPrintLevel(J_DBG, J_NONE);

   return true;
}

void IpoptApplication::RegisterAllIpoptOptions(
   const SmartPtr<RegisteredOptions>& roptions
)
{
   // create Ipopt categories here to have place where to specify its priorities
   roptions->SetRegisteringCategory("Termination", 600000);
   roptions->SetRegisteringCategory("Output", 500000);
   roptions->SetRegisteringCategory("NLP", 480000);
   roptions->SetRegisteringCategory("NLP Scaling", 470000);
   roptions->SetRegisteringCategory("Initialization", 460000);
   roptions->SetRegisteringCategory("Warm Start", 450000);
   roptions->SetRegisteringCategory("Miscellaneous", 400000);
   roptions->SetRegisteringCategory("Barrier Parameter Update", 390000);
   roptions->SetRegisteringCategory("Line Search", 380000);
   roptions->SetRegisteringCategory("Linear Solver", 360000);
   roptions->SetRegisteringCategory("Step Calculation", 350000);
   roptions->SetRegisteringCategory("Restoration Phase", 340000);
   roptions->SetRegisteringCategory("Hessian Approximation", 290000);
   roptions->SetRegisteringCategory("Derivative Checker", 280000);
   roptions->SetRegisteringCategory("MA27 Linear Solver", 199000);
   roptions->SetRegisteringCategory("MA57 Linear Solver", 198000);
   roptions->SetRegisteringCategory("MA77 Linear Solver", 197000);
   roptions->SetRegisteringCategory("MA86 Linear Solver", 196000);
   roptions->SetRegisteringCategory("MA97 Linear Solver", 195000);
   roptions->SetRegisteringCategory("Pardiso (pardiso-project.org) Linear Solver", 190000);
   roptions->SetRegisteringCategory("Pardiso (MKL) Linear Solver", 189000);
   roptions->SetRegisteringCategory("SPRAL Linear Solver", 180000);
   roptions->SetRegisteringCategory("WSMP Linear Solver", 170000);
   roptions->SetRegisteringCategory("Mumps Linear Solver", 160000);
   roptions->SetRegisteringCategory("MA28 Linear Solver", 150000);

   roptions->SetRegisteringCategory("CG Penalty", -400000);
   roptions->SetRegisteringCategory("Inexact Step Computation", -900000);
   roptions->SetRegisteringCategory("Undocumented", -1000000);

   RegisterOptions_Interfaces(roptions);
   RegisterOptions_Algorithm(roptions);
   RegisterOptions_CGPenalty(roptions);
   RegisterOptions_LinearSolvers(roptions);
#ifdef BUILD_INEXACT
   RegisterOptions_Inexact(roptions);
#endif
   roptions->SetRegisteringCategory("");
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
