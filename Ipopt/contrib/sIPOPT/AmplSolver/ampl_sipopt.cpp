// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-11

#include "SensAmplTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "SensApplication.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "SensRegOp.hpp"


int main(int argv, char**argc)
{
  using namespace Ipopt;

  SmartPtr<IpoptApplication> app_ipopt = new IpoptApplication();

  SmartPtr<SensApplication> app_sens = new SensApplication(app_ipopt->Jnlst(),
							   app_ipopt->Options(),
							   app_ipopt->RegOptions());

  // Register sIPOPT options
  RegisterOptions_sIPOPT(app_ipopt->RegOptions());
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

  // Suffixes for sIPOPT
  suffix_handler->AddAvailableSuffix("sens_init_constr", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);

  int n_sens_steps = 0;
  app_ipopt->Options()->GetIntegerValue("n_sens_steps",n_sens_steps,"");
  std::string state;
  std::string state_value;
  std::string state_value_zL;
  std::string state_value_zU;
  std::string sol_state;
  std::string sol_state_zL;
  std::string sol_state_zU;
  for (int k=0; k<n_sens_steps+1; ++k) {
    state          = "sens_state_";
    state_value    = "sens_state_value_";
    sol_state      = "sens_sol_state_";
    state_value_zL = state_value;
    state_value_zU = state_value;
    sol_state_zL   = sol_state;
    sol_state_zU   = sol_state;
    append_Index(state,k);
    append_Index(state_value,k);
    append_Index(state_value_zL,k);
    append_Index(state_value_zU,k);
    append_Index(sol_state,k);
    append_Index(sol_state_zL,k);
    append_Index(sol_state_zU,k);
    sol_state_zL += "_z_L";
    sol_state_zU += "_z_U";
    state_value_zL += "_z_L";
    state_value_zU += "_z_U";
    suffix_handler->AddAvailableSuffix(state, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix(state_value, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(state_value_zL, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(state_value_zU, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state, AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state_zL, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix(sol_state_zU, AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  }

  // for reduced hessian computation
  suffix_handler->AddAvailableSuffix("red_hessian", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);

  // Create AmplOptionsList for sIPOPT AMPL options
  SmartPtr<AmplOptionsList> ampl_options_list = new AmplOptionsList();

  ampl_options_list->AddAmplOption("run_sens", "run_sens",
				   AmplOptionsList::String_Option,
				   "Set to yes if sens algorithm should be run.");
  ampl_options_list->AddAmplOption("compute_red_hessian", "compute_red_hessian",
				   AmplOptionsList::String_Option,
				   "Set to yes if reduced hessian should be computed.");
  ampl_options_list->AddAmplOption("sens_boundcheck", "sens_boundcheck",
				   AmplOptionsList::String_Option,
				   "Set to yes to enable the fix-relax QP adaption to a possible bound check. This feature is experimental.");
  ampl_options_list->AddAmplOption("n_sens_steps", "n_sens_steps",
				   AmplOptionsList::Integer_Option,
				   "Number of sensitivity steps");

  // create AmplSensTNLP from argc.
  SmartPtr<TNLP> sens_tnlp = new SensAmplTNLP(ConstPtr(app_ipopt->Jnlst()),
					      app_ipopt->Options(),
					      argc, suffix_handler, false,
					      ampl_options_list);

  app_sens->Initialize();

  const int n_loops = 1; // make larger for profiling
  for (Index i=0; i<n_loops; i++) {
    retval = app_ipopt->OptimizeTNLP(sens_tnlp);
  }

  /* give pointers to Ipopt algorithm objects to Sens Application */
  app_sens->SetIpoptAlgorithmObjects(app_ipopt, retval);

  app_sens->Run();


  return 0;
}
