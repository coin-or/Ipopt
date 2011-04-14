// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-10-04

#include "MySensTNLP.hpp"

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

  // create AmplSensTNLP from argc. This is an nlp because we are using our own TNLP Adapter
  SmartPtr<TNLP> sens_tnlp = new MySensTNLP();

  app_ipopt->Options()->SetStringValueIfUnset("compute_red_hessian", "yes");

  app_sens->Initialize();

  retval = app_ipopt->OptimizeTNLP(sens_tnlp);

  /* give pointers to Ipopt algorithm objects to Sens Application */
  app_sens->SetIpoptAlgorithmObjects(app_ipopt, retval);

  app_sens->Run();

  return 0;

}
