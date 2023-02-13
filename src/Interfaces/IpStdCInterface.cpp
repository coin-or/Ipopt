// Copyright (C) 2004, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpStdCInterface.h"
#include "IpStdInterfaceTNLP.hpp"
#include "IpOptionsList.hpp"
#include "IpIpoptApplication.hpp"
#include "IpBlas.hpp"
#include "IpSmartPtr.hpp"

struct IpoptProblemInfo
{
   Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
   Ipopt::SmartPtr<Ipopt::StdInterfaceTNLP> tnlp;
   ipindex         n;
   ipnumber*       x_L;
   ipnumber*       x_U;
   ipindex         m;
   ipnumber*       g_L;
   ipnumber*       g_U;
   ipindex         nele_jac;
   ipindex         nele_hess;
   ipindex         index_style;
   Eval_F_CB       eval_f;
   Eval_G_CB       eval_g;
   Eval_Grad_F_CB  eval_grad_f;
   Eval_Jac_G_CB   eval_jac_g;
   Eval_H_CB       eval_h;
   Intermediate_CB intermediate_cb;
   ipnumber        obj_scaling;
   ipnumber*       x_scaling;
   ipnumber*       g_scaling;
};

IpoptProblem CreateIpoptProblem(
   ipindex        n,
   ipnumber*      x_L,
   ipnumber*      x_U,
   ipindex        m,
   ipnumber*      g_L,
   ipnumber*      g_U,
   ipindex        nele_jac,
   ipindex        nele_hess,
   ipindex        index_style,
   Eval_F_CB      eval_f,
   Eval_G_CB      eval_g,
   Eval_Grad_F_CB eval_grad_f,
   Eval_Jac_G_CB  eval_jac_g,
   Eval_H_CB      eval_h
)
{
   // make sure input is Ok
   if( n < 1 || m < 0 || !x_L || !x_U || (m > 0 && (!g_L || !g_U)) || (m == 0 && nele_jac != 0)
       || (m > 0 && nele_jac < 1) || nele_hess < 0 || !eval_f || !eval_grad_f || (m > 0 && (!eval_g || !eval_jac_g)) )
   {
      return NULL;
   }

   IpoptProblem retval = new IpoptProblemInfo;

   retval->tnlp = NULL;

   retval->n   = n;
   retval->x_L = new ipnumber[n];
   Ipopt::IpBlasCopy(n, x_L, 1, retval->x_L, 1);

   retval->x_U = new ipnumber[n];
   Ipopt::IpBlasCopy(n, x_U, 1, retval->x_U, 1);

   retval->m = m;
   if( m > 0 )
   {
      retval->g_L = new ipnumber[m];
      Ipopt::IpBlasCopy(m, g_L, 1, retval->g_L, 1);

      retval->g_U = new ipnumber[m];
      Ipopt::IpBlasCopy(m, g_U, 1, retval->g_U, 1);
   }
   else
   {
      retval->g_L = NULL;
      retval->g_U = NULL;
   }

   retval->app = new Ipopt::IpoptApplication();

   retval->nele_jac = nele_jac;
   retval->nele_hess = nele_hess;
   retval->index_style = index_style;
   retval->eval_f = eval_f;
   retval->eval_g = eval_g;
   retval->eval_grad_f = eval_grad_f;
   retval->eval_jac_g = eval_jac_g;
   retval->eval_h = eval_h;
   retval->intermediate_cb = NULL;
   retval->obj_scaling = 1;
   retval->x_scaling = NULL;
   retval->g_scaling = NULL;

   retval->app->RethrowNonIpoptException(false);

   return retval;
}

void FreeIpoptProblem(
   IpoptProblem ipopt_problem
)
{
   ipopt_problem->app = NULL;

   delete[] ipopt_problem->x_L;
   delete[] ipopt_problem->x_U;
   delete[] ipopt_problem->g_L;
   delete[] ipopt_problem->g_U;
   delete[] ipopt_problem->x_scaling;
   delete[] ipopt_problem->g_scaling;

   delete ipopt_problem;
}

bool AddIpoptStrOption(
   IpoptProblem ipopt_problem,
   char*        keyword,
   char*        val
)
{
   return ipopt_problem->app->Options()->SetStringValue(keyword, val);
}

bool AddIpoptNumOption(
   IpoptProblem ipopt_problem,
   char*        keyword,
   ipnumber     val
)
{
   return ipopt_problem->app->Options()->SetNumericValue(keyword, val);
}

bool AddIpoptIntOption(
   IpoptProblem ipopt_problem,
   char*        keyword,
   ipindex      val
)
{
   return ipopt_problem->app->Options()->SetIntegerValue(keyword, val);
}

bool OpenIpoptOutputFile(
   IpoptProblem ipopt_problem,
   char*        file_name,
   int          print_level
)
{
   return ipopt_problem->app->OpenOutputFile(file_name, Ipopt::EJournalLevel(print_level));
}

bool SetIpoptProblemScaling(
   IpoptProblem ipopt_problem,
   ipnumber     obj_scaling,
   ipnumber*    x_scaling,
   ipnumber*    g_scaling
)
{
   ipopt_problem->obj_scaling = obj_scaling;

   if( x_scaling )
   {
      if( !ipopt_problem->x_scaling )
      {
         ipopt_problem->x_scaling = new ipnumber[ipopt_problem->n];
      }
      Ipopt::IpBlasCopy(ipopt_problem->n, x_scaling, 1, ipopt_problem->x_scaling, 1);
   }
   else
   {
      delete[] ipopt_problem->x_scaling;
      ipopt_problem->x_scaling = NULL;
   }

   if( g_scaling )
   {
      if( !ipopt_problem->g_scaling )
      {
         ipopt_problem->g_scaling = new ipnumber[ipopt_problem->m];
      }
      Ipopt::IpBlasCopy(ipopt_problem->m, g_scaling, 1, ipopt_problem->g_scaling, 1);
   }
   else
   {
      delete[] ipopt_problem->g_scaling;
      ipopt_problem->g_scaling = NULL;
   }

   return true;
}

bool SetIntermediateCallback(
   IpoptProblem    ipopt_problem,
   Intermediate_CB intermediate_cb
)
{
   ipopt_problem->intermediate_cb = intermediate_cb;

   return true;
}

enum ApplicationReturnStatus IpoptSolve(
   IpoptProblem ipopt_problem,
   ipnumber*    x,
   ipnumber*    g,
   ipnumber*    obj_val,
   ipnumber*    mult_g,
   ipnumber*    mult_x_L,
   ipnumber*    mult_x_U,
   UserDataPtr  user_data
)
{
   // Initialize and process options
   Ipopt::ApplicationReturnStatus retval = ipopt_problem->app->Initialize();
   if( retval != Ipopt::Solve_Succeeded )
   {
      return ApplicationReturnStatus(retval);
   }

   if( !x )
   {
      ipopt_problem->app->Jnlst()->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "Error: Array x with starting point information is NULL.");
      return ApplicationReturnStatus(Ipopt::Invalid_Problem_Definition);
   }

   // Copy the starting point information
   ipnumber* start_x = new ipnumber[ipopt_problem->n];
   Ipopt::IpBlasCopy(ipopt_problem->n, x, 1, start_x, 1);

   ipnumber* start_lam = NULL;
   if( mult_g )
   {
      start_lam = new ipnumber[ipopt_problem->m];
      Ipopt::IpBlasCopy(ipopt_problem->m, mult_g, 1, start_lam, 1);
   }

   ipnumber* start_z_L = NULL;
   if( mult_x_L )
   {
      start_z_L = new ipnumber[ipopt_problem->n];
      Ipopt::IpBlasCopy(ipopt_problem->n, mult_x_L, 1, start_z_L, 1);
   }

   ipnumber* start_z_U = NULL;
   if( mult_x_U )
   {
      start_z_U = new ipnumber[ipopt_problem->n];
      Ipopt::IpBlasCopy(ipopt_problem->n, mult_x_U, 1, start_z_U, 1);
   }

   Ipopt::ApplicationReturnStatus status;
   try
   {
      // Create the original nlp
      ipopt_problem->tnlp = new Ipopt::StdInterfaceTNLP(ipopt_problem->n, ipopt_problem->x_L, ipopt_problem->x_U,
            ipopt_problem->m, ipopt_problem->g_L, ipopt_problem->g_U,
            ipopt_problem->nele_jac, ipopt_problem->nele_hess, ipopt_problem->index_style,
            start_x, start_lam, start_z_L, start_z_U,
            ipopt_problem->eval_f, ipopt_problem->eval_g, ipopt_problem->eval_grad_f, ipopt_problem->eval_jac_g, ipopt_problem->eval_h,
            ipopt_problem->intermediate_cb,
            x, mult_x_L, mult_x_U, g, mult_g, obj_val, user_data,
            ipopt_problem->obj_scaling, ipopt_problem->x_scaling, ipopt_problem->g_scaling);
      status = ipopt_problem->app->OptimizeTNLP(ipopt_problem->tnlp);
   }
   catch( Ipopt::INVALID_STDINTERFACE_NLP& exc )
   {
      exc.ReportException(*ipopt_problem->app->Jnlst(), Ipopt::J_ERROR);
      status = Ipopt::Invalid_Problem_Definition;
   }
   catch( Ipopt::IpoptException& exc )
   {
      exc.ReportException(*ipopt_problem->app->Jnlst(), Ipopt::J_ERROR);
      status = Ipopt::Unrecoverable_Exception;
   }
   ipopt_problem->tnlp = NULL;

   delete[] start_x;
   delete[] start_lam;
   delete[] start_z_L;
   delete[] start_z_U;

   return ApplicationReturnStatus(status);
}

bool GetIpoptCurrentIterate(
   IpoptProblem    ipopt_problem,
   bool            scaled,
   ipindex         n,
   ipnumber*       x,
   ipnumber*       z_L,
   ipnumber*       z_U,
   ipindex         m,
   ipnumber*       g,
   ipnumber*       lambda
)
{
   if( IsNull(ipopt_problem->tnlp) )
   {
      return false;
   }

   return ipopt_problem->tnlp->get_curr_iterate(scaled, n, x, z_L, z_U, m, g, lambda);
}

bool GetIpoptCurrentViolations(
   IpoptProblem  ipopt_problem,
   bool          scaled,
   ipindex       n,
   ipnumber*     x_L_violation,
   ipnumber*     x_U_violation,
   ipnumber*     compl_x_L,
   ipnumber*     compl_x_U,
   ipnumber*     grad_lag_x,
   ipindex       m,
   ipnumber*     nlp_constraint_violation,
   ipnumber*     compl_g
)
{
   if( IsNull(ipopt_problem->tnlp) )
   {
      return false;
   }

   return ipopt_problem->tnlp->get_curr_violations(scaled != 0, n, x_L_violation, x_U_violation, compl_x_L, compl_x_U, grad_lag_x, m, nlp_constraint_violation, compl_g);
}
