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

  SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(app->Jnlst(), argc);
  ApplicationReturnStatus retval = app->OptimizeTNLP(ampl_tnlp);

  // finalize_solution method in AmplTNLP writes the solution file

  return (int)retval;
}


