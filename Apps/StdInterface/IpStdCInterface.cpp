// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpUtils.hpp"
#include "IpStdCInterface.h"
#include "IpStdInterfaceTNLP.hpp"
#include "IpOptionsList.hpp"
#include "IpIpoptApplication.hpp"

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

  // Create an IpoptApplication
  SmartPtr<IpoptApplication> app =
    new IpoptApplication();

  // Get the options
  OptionsList* options = (OptionsList*) options_;
  if (options) {
    // Copy the options given by the C interface user to those in the
    // Ipopt Application
    *app->Options() = *options;
  }

  // Create the original nlp
  SmartPtr<TNLP> tnlp =
    new StdInterfaceTNLP(n, x_L, x_U, m, g_L, g_U, nele_jac,
                         nele_hess, x, mult_g, mult_x_L, mult_x_U,
                         eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h,
                         x_, mult_x_L_, mult_x_U_, g, mult_g_, obj_val, user_data);

  ApplicationReturnStatus status = app->OptimizeTNLP(tnlp);

  return (::Int) status;
}

