// Copyright (C) 2004, International Business Machines and others.
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

namespace Ipopt
{
   IpoptApplication::IpoptApplication()
      :
   output_file_(""),
   output_file_print_level_(J_NONE),
   read_params_dat_(true),
   report_solve_status_(true),
   report_solution_(true),
   force_report_solution_to_console_(false),
   report_statistics_(true)
   {
   }

   IpoptApplication::~IpoptApplication()
   {
   }

   ApplicationReturnStatus IpoptApplication::OptimizeTNLP(const SmartPtr<TNLP>& nlp,
      const SmartPtr<OptionsList> additional_options)
   {
    SmartPtr<NLP> nlp_adapter =
      new TNLPAdapter(GetRawPtr(nlp));

    return OptimizeNLP(nlp_adapter, additional_options);
   }

   ApplicationReturnStatus IpoptApplication::OptimizeNLP(const SmartPtr<NLP>& nlp,
      const SmartPtr<OptionsList> additional_options)
   {
      ApplicationReturnStatus retValue = Solve_Succeeded;
  SmartPtr<Journalist> jnlst = new Journalist();

# ifdef IP_DEBUG
  DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst));
# endif

  Journal* jrnl = jnlst->AddJournal("ConsoleStdOut", "stdout", J_SUMMARY);
  jrnl->SetPrintLevel(J_DBG, J_NONE);

  if (output_file_ != "" && output_file_print_level_ != J_NONE) {
     jnlst->AddJournal("IPOPT_OUT", output_file_, output_file_print_level_);
  }

# ifdef IP_DEBUG
  jrnl = jnlst->AddJournal("Debug", "debug.out", J_DETAILED);
  jrnl->SetPrintLevel(J_DBG, J_ALL);
# endif


  // Get the options
  SmartPtr<OptionsList> options;
  if (IsValid(additional_options)){
     options = new OptionsList(*additional_options);
  }
  else {
     options = new OptionsList();
  }

  if (read_params_dat_) {
    FILE* fp_options = fopen("PARAMS.DAT", "r");
    if (fp_options) {
      // PARAMS.DAT exists, read the content
      options->ReadFromFile(*jnlst, fp_options);
      fclose(fp_options);
      fp_options=NULL;
    }
  }

  try {
    SmartPtr<IpoptNLP> ip_nlp =
      new OrigIpoptNLP(ConstPtr(jnlst), GetRawPtr(nlp));

    // Create the IpoptData
    SmartPtr<IpoptData> ip_data = new IpoptData();

    // Create the IpoptCalculators
    SmartPtr<IpoptCalculatedQuantities> ip_cq
    = new IpoptCalculatedQuantities(ip_nlp, ip_data);

    // Create the Algorithm object
    SmartPtr<IpoptAlgorithm> alg 
       = AlgorithmBuilder::BuildBasicAlgorithm(*jnlst, *options, "");

    // Set up the algorithm
    alg->Initialize(*jnlst, *ip_nlp, *ip_data, *ip_cq, *options, "");

    // Run the algorithm
    IpoptAlgorithm::SolverReturn status = alg->Optimize();

   if (status == IpoptAlgorithm::SUCCESS) {
      if (report_solve_status_) {
      jnlst->Printf(J_SUMMARY, J_MAIN, "\n\nIpopt Optimize Successful...\n");
      }
   }
   else if (status == IpoptAlgorithm::MAXITER_EXCEEDED) {
      retValue = Maximum_Iterations_Exceeded;
      if (report_solve_status_) {
         jnlst->Printf(J_SUMMARY, J_MAIN, "\n\nERROR: Maximum Number of Iterations Exceeded.\n");
      }
   }
   else if (status == IpoptAlgorithm::FAILED) {
      retValue = Solve_Failed;
      if (report_solve_status_) {
      jnlst->Printf(J_SUMMARY, J_MAIN, "\n\nERROR: Algorithm Failed - Check detailed output.\n");
      }
   }
   else {
      retValue = Internal_Error;
      if (report_solve_status_) {
         jnlst->Printf(J_SUMMARY, J_MAIN, "\n\nINTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
      }
   }

    EJournalLevel vector_report_level = J_VECTOR;
    if (status == IpoptAlgorithm::SUCCESS && (report_solution_ || force_report_solution_to_console_)) {
       if (force_report_solution_to_console_) {
          vector_report_level = J_SUMMARY;
       }
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "\n\nOptimal solution found! \n");
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "Optimal Objective Value = %.16E\n", ip_cq->curr_f());
      jnlst->PrintVector(vector_report_level, J_SOLUTION, "x", *ip_data->curr_x());
      jnlst->PrintVector(vector_report_level, J_SOLUTION, "y_c", *ip_data->curr_y_c());
      jnlst->PrintVector(vector_report_level, J_SOLUTION, "y_d", *ip_data->curr_y_d());
      jnlst->PrintVector(vector_report_level, J_SOLUTION, "z_L", *ip_data->curr_z_L());
      jnlst->PrintVector(vector_report_level, J_SOLUTION, "z_U", *ip_data->curr_z_U());
      jnlst->PrintVector(vector_report_level, J_SOLUTION, "v_L", *ip_data->curr_v_L());
      jnlst->PrintVector(vector_report_level, J_SOLUTION, "v_U", *ip_data->curr_v_U());
    }

    if (report_statistics_) {
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "Number of Iterations = %d\n", ip_data->iter_count());
    }
  }
  catch(IpoptException& exc) {
    exc.ReportException(*jnlst);
    retValue = Solve_Failed;
  }
  catch(...) {
    IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
    exc.ReportException(*jnlst);
    retValue = NonIpopt_Exception_Thrown;
  }
  return retValue;

   }
} // namespace Ipopt



