// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpUtils.hpp"
#include "IpAlgBuilder.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpStdInterfaceTNLP.hpp"
#include "IpTNLPAdapter.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

OptionsPtr Ipopt_NewOptions()
{
  Ipopt::OptionsList* options = new Ipopt::OptionsList();
  return (OptionsPtr) options;
}

void Ipopt_DeleteOptions(OptionsPtr options)
{
  Ipopt::OptionsList* optlist = (Ipopt::OptionsList*)options;
  delete optlist;
}

Bool Ipopt_AddOption(OptionsPtr options, char* keyword, char* val)
{
  Ipopt::OptionsList* optlist = (Ipopt::OptionsList*)options;
  std::string tag(keyword);
  std::string value(val);

  optlist->SetValue(tag, value);

  return (Bool)true;
}

Bool Ipopt_AddNumOption(OptionsPtr options, char* keyword, Number val)
{
  Ipopt::OptionsList* optlist = (Ipopt::OptionsList*)options;
  std::string tag(keyword);
  Ipopt::Number value=val;

  optlist->SetNumericValue(tag, value);

  return (Bool)true;
}

Bool Ipopt_AddIntOption(OptionsPtr options, char* keyword, Int val)
{
  Ipopt::OptionsList* optlist = (Ipopt::OptionsList*)options;
  std::string tag(keyword);
  Ipopt::Index value=val;

  optlist->SetIntegerValue(tag, value);

  return (Bool)true;
}

Int IpoptSolve(Index n,
               Number* x_,
               Number* x_L_,
               Number* x_U_,
               Index m,
               Number* g,
               Number* g_L_,
               Number* g_U_,
               Index nele_jac,
               Index nele_hess,
               Number* obj_val,
               Number* mult_g_,
               Number* mult_x_L_,
               Number* mult_x_U_,
               Eval_F_CB eval_f,
               Eval_G_CB eval_g,
               Eval_Grad_F_CB eval_grad_f,
               Eval_Jac_G_CB eval_jac_g,
               Eval_H_CB eval_h,
               OptionsPtr options_,
               UserDataPtr user_data)
{
  using namespace Ipopt;

  // Since those arrays are assumed not to change during the optimization
  // make sure we are not accidentially changing them.
  const ::Number* x = x_;
  const ::Number* x_L = x_L_;
  const ::Number* x_U = x_U_;
  const ::Number* g_L = g_L_;
  const ::Number* g_U = g_U_;
  const ::Number* mult_g = mult_g_;
  const ::Number* mult_x_L = mult_x_L_;
  const ::Number* mult_x_U = mult_x_U_;

  // Get the journalist
  ::Int retValue = 0;
  SmartPtr<Journalist> jnlst = new Journalist();

# ifdef IP_DEBUG

  DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst));
# endif

  Journal* jrnl = jnlst->AddJournal("ConsoleStdOut", "stdout", J_SUMMARY);
  jrnl->SetPrintLevel(J_DBG, J_NONE);

# ifdef IP_DEBUG

  jrnl = jnlst->AddJournal("Debug", "debug.out", J_DETAILED);
  jrnl->SetPrintLevel(J_DBG, J_ALL);
# endif

  //  jrnl = jnlst->AddJournal("All", "all.out", J_ALL);

  // Get the options
  OptionsList* options = (OptionsList*) options_;
  bool created_options=false;
  if (!options) {
    created_options=true;
    options = new OptionsList();
    FILE* fp_options = fopen("PARAMS.DAT", "r");
    if (fp_options) {
      // PARAMS.DAT exists, read the content
      options->ReadFromFile(*jnlst, fp_options);
      fclose(fp_options);
      fp_options=NULL;
    }
  }

  try {
    // Create the original nlp and the IpoptNLP
    SmartPtr<StdInterfaceTNLP> stdinterface_nlp =
      new StdInterfaceTNLP(ConstPtr(jnlst), n, x_L, x_U, m, g_L, g_U, nele_jac,
                           nele_hess, x, mult_g, mult_x_L, mult_x_U,
                           eval_f, eval_g, eval_grad_f, eval_jac_g,
                           eval_h, user_data);

    SmartPtr<TNLPAdapter> tnlpadapter =
      new TNLPAdapter(GetRawPtr(stdinterface_nlp));
    SmartPtr<IpoptNLP> ip_nlp =
      new OrigIpoptNLP(ConstPtr(jnlst), GetRawPtr(tnlpadapter));
    //      new OrigIpoptNLP(ConstPtr(jnlst), GetRawPtr(tnlpadapter));

    // Create the IpoptData
    SmartPtr<IpoptData> ip_data = new IpoptData();

    // Create the IpoptCalculators
    SmartPtr<IpoptCalculatedQuantities> ip_cq
    = new IpoptCalculatedQuantities(ip_nlp, ip_data);

    // Create the Algorithm object
    SmartPtr<IpoptAlgorithm> alg = AlgorithmBuilder(*jnlst, *options, "");

    // Set up the algorithm
    alg->Initialize(*jnlst, *ip_nlp, *ip_data, *ip_cq, *options, "");

    // Run the algorithm
    IpoptAlgorithm::SolverReturn status = alg->Optimize();

    if (status == IpoptAlgorithm::SUCCESS) {
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "\n\nOptimal solution found! \n");
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "Optimal Objective Value = %.16E\n", ip_cq->curr_f());
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "Number of Iterations = %d\n", ip_data->iter_count());
      jnlst->PrintVector(J_VECTOR, J_SOLUTION, "x", *ip_data->curr_x());
      jnlst->PrintVector(J_VECTOR, J_SOLUTION, "y_c", *ip_data->curr_y_c());
      jnlst->PrintVector(J_VECTOR, J_SOLUTION, "y_d", *ip_data->curr_y_d());
      jnlst->PrintVector(J_VECTOR, J_SOLUTION, "z_L", *ip_data->curr_z_L());
      jnlst->PrintVector(J_VECTOR, J_SOLUTION, "z_U", *ip_data->curr_z_U());
      jnlst->PrintVector(J_VECTOR, J_SOLUTION, "v_L", *ip_data->curr_v_L());
      jnlst->PrintVector(J_VECTOR, J_SOLUTION, "v_U", *ip_data->curr_v_U());

      retValue = 0;
    }
    else {
      jnlst->Printf(J_ERROR, J_MAIN, "Sorry, things failed ?!?!\n");
      retValue = 1;
    }

    // Copy values into user-provided output arrays
    tnlpadapter->ResortX(*ip_data->curr_x(), x_);
    if (g) {
      tnlpadapter->ResortG(*ip_cq->curr_c(), *ip_cq->curr_d(), g);
    }
    if (obj_val) {
      *obj_val = ip_cq->curr_f();
    }
    if (mult_g_) {
      tnlpadapter->ResortG(*ip_data->curr_y_c(), *ip_data->curr_y_d(), mult_g_);
    }
    if (mult_x_L_ || mult_x_U_) {
      tnlpadapter->ResortBnds(*ip_data->curr_z_L(), mult_x_L_,
                              *ip_data->curr_z_U(), mult_x_U_);
    }
  }
  catch(IpoptException& exc) {
    exc.ReportException(*jnlst);
    retValue = 1;
  }
  catch(...) {
    IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
    exc.ReportException(*jnlst);
    retValue = 2;
  }

  if (created_options) {
    delete options;
  }

  return retValue;
}

