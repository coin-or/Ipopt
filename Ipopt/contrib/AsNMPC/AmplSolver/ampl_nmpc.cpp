// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-11

#include "AsAmplNmpcTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "AsNMPCApplication.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "AsAsNMPCRegOp.hpp"


int main(int argv, char**argc)
{
  using namespace Ipopt;

  SmartPtr<IpoptApplication> app_ipopt = new IpoptApplication();

  SmartPtr<NmpcApplication> app_nmpc = new NmpcApplication(app_ipopt->Jnlst(),
							   app_ipopt->Options(),
							   app_ipopt->RegOptions());

  // Register AsNMPC options
  RegisterOptions_AsNMPC(app_ipopt->RegOptions());
  app_ipopt->Options()->SetRegisteredOptions(app_ipopt->RegOptions());

  // Call Initialize the first time to create a journalist, but ignore
  // any options file
  ApplicationReturnStatus retval;
  retval = app_ipopt->Initialize("");
  if (retval != Solve_Succeeded) {
    //printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
    exit(-100);
  }
    
  app_ipopt->Initialize();

  // prepare suffixes, or metadata ...
  SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
  // Modified for warm-start from AMPL
  suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

  suffix_handler->AddAvailableSuffix("parameter", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("nominal_value", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("perturbed_value", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

  // Suffixes for NMPC
  suffix_handler->AddAvailableSuffix("nmpc_init_constr", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);

  int n_nmpc_steps = 0;
  app_ipopt->Options()->GetIntegerValue("n_nmpc_steps",n_nmpc_steps,"");
  std::string state;
  std::string state_value;
  std::string sol_state;
  std::string sol_state_zL;
  std::string sol_state_zU;
  for (int k=0; k<n_nmpc_steps; ++k) {
    state        = "nmpc_state_";
    state_value  = "nmpc_state_value_";
    sol_state    = "nmpc_sol_state_";
    sol_state_zL = sol_state;
    sol_state_zU = sol_state;
    append_Index(state,k);
    append_Index(state_value,k);
    append_Index(sol_state,k);
    append_Index(sol_state_zL,k);
    append_Index(sol_state_zU,k);
    sol_state_zL += "z_L";
    sol_state_zU += "z_U";
    suffix_handler->AddAvailableSuffix(state, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix(state_value, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state, AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state_zL, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state_zU, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    // what's this good for??  suffix_handler->AddAvailableSuffix("nmpc_state_value_1_z_L", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  }

  // for reduced hessian computation
  suffix_handler->AddAvailableSuffix("red_hessian", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);

  // Create AmplOptionsList for AsNMPC AMPL options
  SmartPtr<AmplOptionsList> ampl_options_list = new AmplOptionsList();

  ampl_options_list->AddAmplOption("run_nmpc", "run_nmpc",
				   AmplOptionsList::String_Option,
				   "Set to yes if nmpc algorithm should be run.");
  ampl_options_list->AddAmplOption("compute_red_hessian", "compute_red_hessian", 
				   AmplOptionsList::String_Option,
				   "Set to yes if reduced hessian should be computed.");

  // create AmplSensTNLP from argc.
  SmartPtr<TNLP> nmpc_tnlp = new AmplNmpcTNLP(ConstPtr(app_ipopt->Jnlst()),
					      app_ipopt->Options(),
					      argc, suffix_handler, false,
					      ampl_options_list);

  app_nmpc->Initialize();

  const int n_loops = 1; // make larger for profiling
  for (Index i=0; i<n_loops; i++) {
    retval = app_ipopt->OptimizeTNLP(nmpc_tnlp);
  }

  /* give pointers to Ipopt algorithm objects to NMPC Application */
  app_nmpc->SetIpoptAlgorithmObjects(app_ipopt, retval);

  app_nmpc->Run();


  return 0;
}
