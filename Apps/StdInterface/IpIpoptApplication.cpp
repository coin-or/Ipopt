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
      read_params_dat_(true),
      report_solution_(true),
      force_report_solution_to_console_(false),
      report_statistics_(true),
      jnlst_(new Journalist()),
      options_(new OptionsList())
  {}

  IpoptApplication::~IpoptApplication()
  {}

  ApplicationReturnStatus IpoptApplication::OptimizeTNLP(const SmartPtr<TNLP>& nlp)
  {
    SmartPtr<NLP> nlp_adapter =
      new TNLPAdapter(GetRawPtr(nlp));

    return OptimizeNLP(nlp_adapter);
  }

  ApplicationReturnStatus 
  IpoptApplication::OptimizeTNLP(const SmartPtr<TNLP>& nlp, 
				 SmartPtr<IpoptData>& ip_data,
				 SmartPtr<IpoptCalculatedQuantities>& ip_cq)
  {
    SmartPtr<NLP> nlp_adapter =
      new TNLPAdapter(GetRawPtr(nlp));
    
    return OptimizeNLP(nlp_adapter, ip_data, ip_cq);
  }

  ApplicationReturnStatus IpoptApplication::OptimizeNLP(const SmartPtr<NLP>& nlp)
  {
    SmartPtr<IpoptData> ip_data = NULL;
    SmartPtr<IpoptCalculatedQuantities> ip_cq = NULL;

    return OptimizeNLP(nlp, ip_data, ip_cq);

  }

  ApplicationReturnStatus 
  IpoptApplication::OptimizeNLP(const SmartPtr<NLP>& nlp,
				SmartPtr<IpoptData>& ip_data,
				SmartPtr<IpoptCalculatedQuantities>& ip_cq)
  {
    ApplicationReturnStatus retValue = Solve_Succeeded;

# ifdef IP_DEBUG

    DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst_));
# endif

    Journal* stdout_jrnl =
      jnlst_->AddJournal("ConsoleStdOut", "stdout", J_SUMMARY);
    stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);

# ifdef IP_DEBUG

    Journal* dbg_jrnl = jnlst_->AddJournal("Debug", "debug.out", J_DETAILED);
    dbg_jrnl->SetPrintLevel(J_DBG, J_ALL);
# endif


    // Get the options
    if (read_params_dat_) {
      FILE* fp_options = fopen("PARAMS.DAT", "r");
      if (fp_options) {
        // PARAMS.DAT exists, read the content
        options_->ReadFromFile(*jnlst_, fp_options);
        fclose(fp_options);
        fp_options=NULL;
      }
    }

    // Set printlevel for stdout
    Index ivalue;
    EJournalLevel print_level = J_SUMMARY;
    if (options_->GetIntegerValue("print_level", ivalue, "")) {
      ivalue = Max(0, Min(ivalue, ((Index)J_LAST_LEVEL)-1));
      print_level = (EJournalLevel)ivalue;
      stdout_jrnl->SetAllPrintLevels(print_level);
      stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
    }

    // Open an output file if required
    std::string output_filename;
    if (options_->GetValue("output_file", output_filename, "")) {
      EJournalLevel file_print_level = J_SUMMARY;
      if (options_->GetIntegerValue("file_print_level", ivalue, "")) {
	ivalue = Max(0, Min(ivalue, ((Index)J_LAST_LEVEL)-1));
        file_print_level = (EJournalLevel)ivalue;
      }
      else {
	file_print_level = print_level;
      }
      Journal* file_jrnl = jnlst_->AddJournal("OutputFile", output_filename.c_str(), file_print_level);
      file_jrnl->SetPrintLevel(J_DBG, J_NONE);
    }

    try {
      SmartPtr<IpoptNLP> ip_nlp =
        new OrigIpoptNLP(ConstPtr(jnlst_), GetRawPtr(nlp));

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
      IpoptAlgorithm::SolverReturn status = alg->Optimize();

      EJournalLevel vector_report_level = J_VECTOR;
      if (report_solution_ || force_report_solution_to_console_) {
        if (force_report_solution_to_console_) {
          vector_report_level = J_SUMMARY;
        }
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "\nNumber of Iterations    = %d\n",
                       ip_data->iter_count());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Objective Value         = %23.16e\n",
                       ip_cq->curr_f());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Primal Infeasibility    = %23.16e\n",
                       ip_cq->curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_MAX));
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Dual Infeasibility      = %23.16e\n",
                       ip_cq->curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_MAX));
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Complementarity         = %23.16e\n",
                       ip_cq->curr_complementarity(0., IpoptCalculatedQuantities::NORM_MAX));

        jnlst_->PrintVector(vector_report_level, J_SOLUTION, "x", *ip_data->curr_x());
        jnlst_->PrintVector(vector_report_level, J_SOLUTION, "y_c", *ip_data->curr_y_c());
        jnlst_->PrintVector(vector_report_level, J_SOLUTION, "y_d", *ip_data->curr_y_d());
        jnlst_->PrintVector(vector_report_level, J_SOLUTION, "z_L", *ip_data->curr_z_L());
        jnlst_->PrintVector(vector_report_level, J_SOLUTION, "z_U", *ip_data->curr_z_U());
        jnlst_->PrintVector(vector_report_level, J_SOLUTION, "v_L", *ip_data->curr_v_L());
        jnlst_->PrintVector(vector_report_level, J_SOLUTION, "v_U", *ip_data->curr_v_U());
      }

      if (report_statistics_) {
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "\nNumber of objective function evaluations             = %d\n",
                       ip_nlp->f_evals());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Number of equality constraint evaluations            = %d\n",
                       ip_nlp->c_evals());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Number of inequality constraint evaluations          = %d\n",
                       ip_nlp->d_evals());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Number of equality constraint Jacobian evaluations   = %d\n",
                       ip_nlp->jac_c_evals());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Number of inequality constraint Jacobian evaluations = %d\n",
                       ip_nlp->jac_d_evals());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Number of Lagrangian Hessian evaluations             = %d\n",
                       ip_nlp->h_evals());
      }

      // Write EXIT message
      if (status == IpoptAlgorithm::SUCCESS) {
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Optimal Solution Found.\n");
      }
      else if (status == IpoptAlgorithm::MAXITER_EXCEEDED) {
        retValue = Maximum_Iterations_Exceeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum Number of Iterations Exceeded.\n");
      }
      else if (status == IpoptAlgorithm::FAILED) {
        retValue = Solve_Failed;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Algorithm Failed - Check detailed output.\n");
      }
      else {
        retValue = Internal_Error;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
      }

      nlp->FinalizeSolution(retValue,
                            *ip_data->curr_x(), *ip_data->curr_z_L(), *ip_data->curr_z_U(),
                            *ip_cq->curr_c(), *ip_cq->curr_d(), *ip_data->curr_y_c(), *ip_data->curr_y_d(),
                            ip_cq->curr_f());
    }
    catch(LOCALLY_INFEASIBILE& exc) {
      exc.ReportException(*jnlst_);

      retValue = Infeasible_Problem_Detected;
    }
    catch(TOO_FEW_DOF& exc) {
      exc.ReportException(*jnlst_);

      retValue = Not_Enough_Degrees_Of_Freedom;
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst_);
      
      retValue = Solve_Failed;
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

    return retValue;

  }
} // namespace Ipopt



