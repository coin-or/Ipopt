// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpUtils.hpp"
#include "AmplTNLP.hpp"
#include "IpIpoptApplication.hpp"

int main(int argv, char**argc)
{
  using namespace Ipopt;

  SmartPtr<IpoptApplication> app = new IpoptApplication();
  app->Jnlst()->Printf(J_ERROR, J_MAIN, "\n\n\n*************************************************************\n");
  app->Jnlst()->Printf(J_ERROR, J_MAIN, "*** Running Ipopt with AMPL Model ***************************\n");
  app->Jnlst()->Printf(J_ERROR, J_MAIN, "*************************************************************\n\n\n");

  // Get the options
  // TODO: Need method in AmplNLP that can fill options with the  AMPL user options
  //...  get options from AMPL
  //... app->Options()->Add...

  // Add the suffix handler for scaling
  SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);

  SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()), argc, suffix_handler);

  ApplicationReturnStatus retval;
  const int n_loops = 1; // make larger for profiling
  for (Index i=0; i<n_loops; i++) {
    retval = app->OptimizeTNLP(ampl_tnlp);
  }

  // finalize_solution method in AmplTNLP writes the solution file

  return (int)retval;
}


