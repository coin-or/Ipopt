// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2004-09-23

#include "IpRestoIterationOutput.hpp"
#include "IpRestoIpoptNLP.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  RestoIterationOutput::RestoIterationOutput(const SmartPtr<OrigIterationOutput>& resto_orig_iteration_output)
      :
      resto_orig_iteration_output_(resto_orig_iteration_output)
  {}

  RestoIterationOutput::~RestoIterationOutput()
  {}

  bool RestoIterationOutput::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    bool retval = true;
    if (IsValid(resto_orig_iteration_output_)) {
      retval = resto_orig_iteration_output_->Initialize(Jnlst(), IpNLP(),
               IpData(), IpCq(),
               options, prefix);
    }
    return retval;
  }

  void RestoIterationOutput::WriteOutput()
  {
    // Get pointers to the Original NLP objects
    const RestoIpoptNLP* resto_ipopt_nlp =
      dynamic_cast<const RestoIpoptNLP*>(&IpNLP());
    DBG_ASSERT(resto_ipopt_nlp);

    SmartPtr<IpoptData> orig_ip_data = &resto_ipopt_nlp->OrigIpData();
    SmartPtr<IpoptNLP> orig_ip_nlp = &resto_ipopt_nlp->OrigIpNLP();
    SmartPtr<IpoptCalculatedQuantities> orig_ip_cq =
      &resto_ipopt_nlp->OrigIpCq();

    // Set the iteration counter for the original NLP to the current value
    Index iter = IpData().iter_count();
    orig_ip_data->Set_iter_count(iter);

    // If a resto_orig_iteration_output object was given, first do the
    // WriteOutput method with that one
    if (IsValid(resto_orig_iteration_output_)) {
      resto_orig_iteration_output_->WriteOutput();
    }

    //////////////////////////////////////////////////////////////////////
    //         First print the summary line for the iteration           //
    //////////////////////////////////////////////////////////////////////

    std::string header =
      " iter     objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls\n";
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Summary of Iteration %d for original NLP:", IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");
    if (iter%10 == 0 && !IsValid(resto_orig_iteration_output_)) {
      // output the header
      Jnlst().Printf(J_SUMMARY, J_MAIN, header.c_str());
    }
    else {
      Jnlst().Printf(J_DETAILED, J_MAIN, header.c_str());
    }

    // For now, just print the total NLP error for the restoration
    // phase problem in the dual infeasibility column
    Number inf_du =
      IpCq().curr_dual_infeasibility(NORM_MAX);

    Number mu = IpData().curr_mu();
    Number dnrm = 0.;
    if (IsValid(IpData().delta_x())) {
      dnrm = Max(IpData().delta_x()->Amax(), IpData().delta_s()->Amax());
    }

    // Set  the trial  values  for  the original  Data  object to  the
    // current restoration phase values
    SmartPtr<const Vector> x = IpData().curr_x();
    const CompoundVector* cx =
      dynamic_cast<const CompoundVector*>(GetRawPtr(x));
    DBG_ASSERT(cx);

    SmartPtr<const Vector> x_only = cx->GetComp(0);
    orig_ip_data->SetTrialPrimalVariablesFromPtr(x_only, IpData().curr_s());

    // Compute primal infeasibility
    Number inf_pr = orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
    Number f = orig_ip_cq->trial_f();

    // Retrieve some information set in the different parts of the algorithm
    char info_iter='r';

    Number alpha_primal = IpData().info_alpha_primal();
    char alpha_primal_char = IpData().info_alpha_primal_char();
    Number alpha_dual = IpData().info_alpha_dual();
    Number regu_x = IpData().info_regu_x();
    char regu_x_buf[8];
    char dashes[]="   - ";
    char *regu_x_ptr;
    if (regu_x==.0) {
      regu_x_ptr = dashes;
    }
    else {
      sprintf(regu_x_buf, "%5.1f", log10(regu_x));
      regu_x_ptr = regu_x_buf;
    }
    Index ls_count = IpData().info_ls_count();
    const std::string info_string = IpData().info_string();

    Jnlst().Printf(J_SUMMARY, J_MAIN,
                   "%5d%c %14.7e %7.2e %7.2e %5.1f %7.2e %5s %7.2e %7.2e%c%3d %s\n",
                   iter, info_iter, f, inf_pr, inf_du, log10(mu), dnrm, regu_x_ptr,
                   alpha_dual, alpha_primal, alpha_primal_char,
                   ls_count, info_string.c_str());


    //////////////////////////////////////////////////////////////////////
    //           Now if desired more detail on the iterates             //
    //////////////////////////////////////////////////////////////////////

    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Beginning Iteration %d from the following point:",
                   IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");

    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "Primal infeasibility for restoration phase problem = %.16e\n",
                   IpCq().curr_primal_infeasibility(NORM_MAX));
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "Dual infeasibility for restoration phase problem   = %.16e\n",
                   IpCq().curr_dual_infeasibility(NORM_MAX));

    if (Jnlst().ProduceOutput(J_DETAILED, J_MAIN)) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_x||_inf   = %.16e\n", IpData().curr_x()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_s||_inf   = %.16e\n", IpData().curr_s()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_y_c||_inf = %.16e\n", IpData().curr_y_c()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_y_d||_inf = %.16e\n", IpData().curr_y_d()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_z_L||_inf = %.16e\n", IpData().curr_z_L()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_z_U||_inf = %.16e\n", IpData().curr_z_U()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_v_L||_inf = %.16e\n", IpData().curr_v_L()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_v_U||_inf = %.16e\n", IpData().curr_v_U()->Amax());
    }
    if (Jnlst().ProduceOutput(J_MOREDETAILED, J_MAIN)) {
      if (IsValid(IpData().delta_x())) {
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "\n||delta_x||_inf   = %.16e\n", IpData().delta_x()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_s||_inf   = %.16e\n", IpData().delta_s()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_y_c||_inf = %.16e\n", IpData().delta_y_c()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_y_d||_inf = %.16e\n", IpData().delta_y_d()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_z_L||_inf = %.16e\n", IpData().delta_z_L()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_z_U||_inf = %.16e\n", IpData().delta_z_U()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_v_L||_inf = %.16e\n", IpData().delta_v_L()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_v_U||_inf = %.16e\n", IpData().delta_v_U()->Amax());
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "\nNo search direction has been computed yet.\n");
      }
    }
    if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_x", *IpData().curr_x());
      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_s", *IpData().curr_s());

      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_y_c", *IpData().curr_y_c());
      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_y_d", *IpData().curr_y_d());

      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_slack_x_L", *IpCq().curr_slack_x_L());
      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_slack_x_U", *IpCq().curr_slack_x_U());

      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_z_L", *IpData().curr_z_L());
      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_z_U", *IpData().curr_z_U());

      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_slack_s_L", *IpCq().curr_slack_s_L());
      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_slack_s_U", *IpCq().curr_slack_s_U());

      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_v_L", *IpData().curr_v_L());
      Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_v_U", *IpData().curr_v_U());
    }
    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_MAIN)) {
      Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "curr_grad_lag_x", *IpCq().curr_grad_lag_x());
      Jnlst().PrintVector(J_MOREVECTOR, J_MAIN, "curr_grad_lag_s", *IpCq().curr_grad_lag_s());
      if (IsValid(IpData().delta_x())) {
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_x", *IpData().delta_x());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_s", *IpData().delta_s());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_y_c", *IpData().delta_y_c());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_y_d", *IpData().delta_y_d());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_z_L", *IpData().delta_z_L());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_z_U", *IpData().delta_z_U());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_v_L", *IpData().delta_v_L());
        Jnlst().PrintVector(J_MOREVECTOR, J_MAIN,
                            "delta_v_U", *IpData().delta_v_U());
      }
    }

    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n\n***Current NLP Values for Iteration (Restoration phase problem) %d:\n",
                   IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN, "Objective = %.16e\n", IpCq().curr_f());
    Jnlst().PrintVector(J_VECTOR, J_MAIN, "grad_f", *IpCq().curr_grad_f());
    Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_c", *IpCq().curr_c());
    Jnlst().PrintVector(J_VECTOR, J_MAIN, "curr_d", *IpCq().curr_d());
    Jnlst().PrintVector(J_VECTOR, J_MAIN,
                        "curr_d - curr_s", *IpCq().curr_d_minus_s());

    Jnlst().PrintMatrix(J_MATRIX, J_MAIN, "jac_c", *IpCq().curr_jac_c());
    Jnlst().PrintMatrix(J_MATRIX, J_MAIN, "jac_d", *IpCq().curr_jac_d());
    Jnlst().PrintMatrix(J_MATRIX, J_MAIN, "h", *IpCq().curr_exact_hessian());
    Jnlst().Printf(J_DETAILED, J_MAIN, "\n\n");

  }

} // namespace Ipopt
