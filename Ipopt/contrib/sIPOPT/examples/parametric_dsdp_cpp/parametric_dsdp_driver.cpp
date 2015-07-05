// Copyright 2010, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-01-05

#include "parametricTNLP.hpp"

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
  SmartPtr<TNLP> sens_tnlp = new ParametricTNLP();

  app_ipopt->Options()->SetStringValueIfUnset("run_sens", "yes");
  app_ipopt->Options()->SetIntegerValueIfUnset("n_sens_steps", 1);

  app_sens->Initialize();

  retval = app_ipopt->OptimizeTNLP(sens_tnlp);

  /* give pointers to Ipopt algorithm objects to Sens Application */
  app_sens->SetIpoptAlgorithmObjects(app_ipopt, retval);

  printf("\n");
  printf("#-------------------------------------------\n");
  printf("# Sensitivity without bound checking\n");
  printf("#-------------------------------------------\n");

  app_sens->Run();
  
  Index m = app_sens->nl() ;
  Index n = app_sens->nx() ;
  Index nzl = app_sens->nzl() ;
  Index nzu = app_sens->nzu() ;
  Number *SX = new Number[n] ;
  Number *SL = new Number[m] ;
  Number *SZL = new Number[nzl] ;
  Number *SZU = new Number[nzu] ;

  app_sens->GetSensitivities(SX,SL,SZL,SZU) ;

  printf("\n** Sensitivity vector (Eq. 9 of implementation paper) ** \n") ;
  for (int i=0; i<n; ++i) printf("* ds/dp(x)[%i] = %.14g\n",i+1,SX[i]) ;
  for (int i=0; i<m; ++i) printf("* ds/dp(l)[%i] = %.14g\n",i+1,SL[i]) ;
  for (int i=0; i<nzl; ++i) printf("* ds/dp(zl)[%i] = %.14g\n",i+1,SZL[i]) ;
  for (int i=0; i<nzu; ++i) printf("* ds/dp(zu)[%i] = %.14g\n",i+1,SZU[i]) ;

    
  
  printf("\n");
  printf("#-------------------------------------------\n");
  printf("# Sensitivity with bound checking\n");
  printf("#-------------------------------------------\n");
  app_ipopt->Options()->SetStringValue("sens_boundcheck", "yes");
  app_sens->Run();

  app_sens->GetSensitivities(SX,SL,SZL,SZU) ;

  printf("\n** Sensitivity vector (Eq. 9 of implementation paper) ** \n") ;
  for (int i=0; i<n; ++i) printf("* ds/dp(x)[%i] = %.14g\n",i+1,SX[i]) ;
  for (int i=0; i<m; ++i) printf("* ds/dp(l)[%i] = %.14g\n",i+1,SL[i]) ;
  for (int i=0; i<nzl; ++i) printf("* ds/dp(zl)[%i] = %.14g\n",i+1,SZL[i]) ;
  for (int i=0; i<nzu; ++i) printf("* ds/dp(zu)[%i] = %.14g\n",i+1,SZU[i]) ;


  delete [] SX ;
  delete [] SL ;
  delete [] SZL ;
  delete [] SZU ;

}
