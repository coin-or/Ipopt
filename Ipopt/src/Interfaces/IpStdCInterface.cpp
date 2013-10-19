// Copyright (C) 2004, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpStdCInterface.h"
#include "IpStdInterfaceTNLP.hpp"
#include "IpOptionsList.hpp"
#include "IpIpoptApplication.hpp"

struct IpoptProblemInfo
{
  Index n;
  Number* x_L;
  Number* x_U;
  Index m;
  Number* g_L;
  Number* g_U;
  Index nele_jac;
  Index nele_hess;
  Index index_style;
  Eval_F_CB eval_f;
  Eval_G_CB eval_g;
  Eval_Grad_F_CB eval_grad_f;
  Eval_Jac_G_CB eval_jac_g;
  Eval_H_CB eval_h;
  Intermediate_CB intermediate_cb;
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
  Number obj_scaling;
  Number* x_scaling;
  Number* g_scaling;
};

IpoptProblem CreateIpoptProblem(
  Index n,
  Number* x_L,
  Number* x_U,
  Index m,
  Number* g_L,
  Number* g_U,
  Index nele_jac,
  Index nele_hess,
  Index index_style,
  Eval_F_CB eval_f,
  Eval_G_CB eval_g,
  Eval_Grad_F_CB eval_grad_f,
  Eval_Jac_G_CB eval_jac_g,
  Eval_H_CB eval_h)
{
  // make sure input is Ok
  if (n<1 || m<0 || !x_L || !x_U || (m>0 && (!g_L || !g_U)) ||
      (m==0 && nele_jac != 0) || (m>0 && nele_jac < 1) || nele_hess < 0 ||
      !eval_f || !eval_grad_f || (m>0 && (!eval_g || !eval_jac_g))) {
    return NULL;
  }

  IpoptProblem retval = new IpoptProblemInfo;

  retval->n = n;
  retval->x_L = new Number[n];
  for (Index i=0; i<n; i++) {
    retval->x_L[i] = x_L[i];
  }
  retval->x_U = new Number[n];
  for (Index i=0; i<n; i++) {
    retval->x_U[i] = x_U[i];
  }

  retval->m = m;
  if (m>0) {
    retval->g_L = new Number[m];
    for (Index i=0; i<m; i++) {
      retval->g_L[i] = g_L[i];
    }
    retval->g_U = new Number[m];
    for (Index i=0; i<m; i++) {
      retval->g_U[i] = g_U[i];
    }
  }
  else {
    retval->g_L = NULL;
    retval->g_U = NULL;
  }

  retval->nele_jac = nele_jac;
  retval->nele_hess = nele_hess;
  retval->index_style = index_style;
  retval->eval_f = eval_f;
  retval->eval_g = eval_g;
  retval->eval_grad_f = eval_grad_f;
  retval->eval_jac_g = eval_jac_g;
  retval->eval_h = eval_h;
  retval->intermediate_cb = NULL;

  retval->app = new Ipopt::IpoptApplication();

  retval->obj_scaling = 1;
  retval->x_scaling = NULL;
  retval->g_scaling = NULL;

  retval->app->RethrowNonIpoptException(false);

  return retval;
}

void FreeIpoptProblem(IpoptProblem ipopt_problem)
{
  delete [] ipopt_problem->x_L;
  delete [] ipopt_problem->x_U;
  if (ipopt_problem->m>0) {
    delete [] ipopt_problem->g_L;
    delete [] ipopt_problem->g_U;
  }

  ipopt_problem->app = NULL;

  delete [] ipopt_problem->x_scaling;
  delete [] ipopt_problem->g_scaling;

  delete ipopt_problem;

}


Bool AddIpoptStrOption(IpoptProblem ipopt_problem, char* keyword, char* val)
{
  std::string tag(keyword);
  std::string value(val);
  return (Bool) ipopt_problem->app->Options()->SetStringValue(tag, value);
}

Bool AddIpoptNumOption(IpoptProblem ipopt_problem, char* keyword, Number val)
{
  std::string tag(keyword);
  Ipopt::Number value=val;
  return (Bool) ipopt_problem->app->Options()->SetNumericValue(tag, value);
}

Bool AddIpoptIntOption(IpoptProblem ipopt_problem, char* keyword, Int val)
{
  std::string tag(keyword);
  Ipopt::Index value=val;
  return (Bool) ipopt_problem->app->Options()->SetIntegerValue(tag, value);
}

Bool OpenIpoptOutputFile(IpoptProblem ipopt_problem, char* file_name,
                         Int print_level)
{
  std::string name(file_name);
  Ipopt::EJournalLevel level = Ipopt::EJournalLevel(print_level);
  return (Bool) ipopt_problem->app->OpenOutputFile(name, level);
}

Bool SetIpoptProblemScaling(IpoptProblem ipopt_problem,
                            Number obj_scaling,
                            Number* x_scaling,
                            Number* g_scaling)
{
  ipopt_problem->obj_scaling = obj_scaling;
  if (x_scaling) {
    if (!ipopt_problem->x_scaling) {
      ipopt_problem->x_scaling = new Number[ipopt_problem->n];
    }
    for (::Index i=0; i<ipopt_problem->n; i++) {
      ipopt_problem->x_scaling[i] = x_scaling[i];
    }
  }
  else {
    delete [] ipopt_problem->x_scaling;
    ipopt_problem->x_scaling = NULL;
  }
  if (g_scaling) {
    if (!ipopt_problem->g_scaling) {
      ipopt_problem->g_scaling = new Number[ipopt_problem->m];
    }
    for (::Index i=0; i<ipopt_problem->m; i++) {
      ipopt_problem->g_scaling[i] = g_scaling[i];
    }
  }
  else {
    delete [] ipopt_problem->g_scaling;
    ipopt_problem->g_scaling = NULL;
  }

  return (Bool)true;
}

Bool SetIntermediateCallback(IpoptProblem ipopt_problem,
                             Intermediate_CB intermediate_cb)
{
  ipopt_problem->intermediate_cb = intermediate_cb;
  return (Bool)true;
}


enum ApplicationReturnStatus IpoptSolve(
  IpoptProblem ipopt_problem,
  Number* x,
  Number* g,
  Number* obj_val,
  Number* mult_g,
  Number* mult_x_L,
  Number* mult_x_U,
  UserDataPtr user_data)
{
  using namespace Ipopt;

  // Initialize and process options
  Ipopt::ApplicationReturnStatus retval = ipopt_problem->app->Initialize();
  if (retval!=Ipopt::Solve_Succeeded) {
    return (::ApplicationReturnStatus) retval;
  }

  if (!x) {
    ipopt_problem->app->Jnlst()->Printf(J_ERROR, J_MAIN,
                                        "Error: Array x with starting point information is NULL.");
    return (::ApplicationReturnStatus) Ipopt::Invalid_Problem_Definition;
  }

  // Copy the starting point information
  ::Number* start_x = new ::Number[ipopt_problem->n];
  for (::Index i=0; i<ipopt_problem->n; i++) {
    start_x[i] = x[i];
  }
  ::Number* start_lam = NULL;
  if (mult_g) {
    start_lam = new ::Number[ipopt_problem->m];
    for (::Index i=0; i<ipopt_problem->m; i++) {
      start_lam[i] = mult_g[i];
    }
  }
  ::Number* start_z_L = NULL;
  if (mult_x_L) {
    start_z_L = new ::Number[ipopt_problem->n];
    for (::Index i=0; i<ipopt_problem->n; i++) {
      start_z_L[i] = mult_x_L[i];
    }
  }
  ::Number* start_z_U = NULL;
  if (mult_x_U) {
    start_z_U = new ::Number[ipopt_problem->n];
    for (::Index i=0; i<ipopt_problem->n; i++) {
      start_z_U[i] = mult_x_U[i];
    }
  }

  // Create the original nlp
  SmartPtr<TNLP> tnlp;

  Ipopt::ApplicationReturnStatus status;
  try {
    tnlp = new StdInterfaceTNLP(ipopt_problem->n, ipopt_problem->x_L,
                                ipopt_problem->x_U, ipopt_problem->m,
                                ipopt_problem->g_L, ipopt_problem->g_U,
                                ipopt_problem->nele_jac,
                                ipopt_problem->nele_hess,
                                ipopt_problem->index_style,
                                start_x, start_lam, start_z_L, start_z_U,
                                ipopt_problem->eval_f, ipopt_problem->eval_g,
                                ipopt_problem->eval_grad_f,
                                ipopt_problem->eval_jac_g,
                                ipopt_problem->eval_h,
                                ipopt_problem->intermediate_cb,
                                x, mult_x_L, mult_x_U, g, mult_g,
                                obj_val, user_data,
                                ipopt_problem->obj_scaling,
                                ipopt_problem->x_scaling,
                                ipopt_problem->g_scaling);
    status = ipopt_problem->app->OptimizeTNLP(tnlp);
  }
  catch (INVALID_STDINTERFACE_NLP& exc) {
    exc.ReportException(*ipopt_problem->app->Jnlst(), J_ERROR);
    status = Ipopt::Invalid_Problem_Definition;
  }
  catch( IpoptException& exc ) {
    exc.ReportException(*ipopt_problem->app->Jnlst(), J_ERROR);
    status = Ipopt::Unrecoverable_Exception;
  }

  delete [] start_x;
  delete [] start_lam;
  delete [] start_z_L;
  delete [] start_z_U;

  return (::ApplicationReturnStatus) status;
}

