// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
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
  suffix_handler->AddAvailableSuffix("nmpc_state_0", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("nmpc_state_1", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("nmpc_state_2", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("nmpc_state_3", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("nmpc_state_4", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("nmpc_state_5", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("nmpc_state_6", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);

  suffix_handler->AddAvailableSuffix("nmpc_init_constr", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);

  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_0", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_1", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_2", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_3", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_4", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_5", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_6", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_1", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_2", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_3", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_4", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_5", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_6", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_1", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_2", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_3", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_4", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_5", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_6", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);

  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_1_z_L", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_1_z_L", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_sol_state_1_z_U", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->
    AddAvailableSuffix("nmpc_state_value_1_z_U", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

  // for reduced hessian computation
  suffix_handler->AddAvailableSuffix("red_hessian", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);


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

  // Create AmplOptionsList for AsNMPC AMPL options
  SmartPtr<AmplOptionsList> ampl_options_list = new AmplOptionsList();

  ampl_options_list->AddAmplOption("run_nmpc", "run_nmpc",
				   AmplOptionsList::String_Option,
				   "Set to yes if nmpc algorithm should be run.");

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
