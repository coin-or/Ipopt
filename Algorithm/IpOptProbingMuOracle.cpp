// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2004-11-12

#include "IpOptProbingMuOracle.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  OptProbingMuOracle::OptProbingMuOracle(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      MuOracle(),
      pd_solver_(pd_solver)
  {
    DBG_ASSERT(IsValid(pd_solver_));
  }

  OptProbingMuOracle::~OptProbingMuOracle()
  {}

  bool OptProbingMuOracle::InitializeImpl(const OptionsList& options,
                                          const std::string& prefix)
  {
    // Check for the algorithm options
    Number value;
    Index ivalue;
    if (options.GetNumericValue("sigma_max", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"sigma_max\": This value must be positive.");
      sigma_max_ = value;
    }
    else {
      sigma_max_ = 1e2;
    }

    if (options.GetIntegerValue("quality_function_norm", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue>=1 && ivalue<=4, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"quality_function_norm\": This value must be between 1 and 3.");
      quality_function_norm_ = ivalue;
    }
    else {
      quality_function_norm_ = 2;
    }

    if (options.GetIntegerValue("quality_function_normalized", ivalue, prefix)) {
      quality_function_normalized_ = !(ivalue==0);
    }
    else {
      quality_function_normalized_ = false;
    }

    if (options.GetIntegerValue("quality_function_centrality", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue>=0 && ivalue<=3, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"quality_function_centrality\": This value must be 0, 1, 2 or 3.");
      quality_function_centrality_ = ivalue;
    }
    else {
      quality_function_centrality_ = 1;
    }

    if (options.GetIntegerValue("quality_function_dual_inf", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue>=1 && ivalue<=2, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"quality_function_dual_inf\": This value must be 1 or 2.");
      quality_function_dual_inf_ = ivalue;
    }
    else {
      quality_function_dual_inf_ = 1;
    }

    if (options.GetIntegerValue("max_bisection_steps", ivalue, prefix)) {
      //      ASSERT_EXCEPTION(ivalue>0, OptionsList::OPTION_OUT_OF_RANGE,
      //                       "Option \"max_bisection_steps\": This value must be positive.");
      max_bisection_steps_ = ivalue;
    }
    else {
      max_bisection_steps_ = 4;
    }

    if (options.GetNumericValue("bisection_tol", value, prefix)) {
      ASSERT_EXCEPTION(value > 0 && value <1., OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"bisection_tol\": This value must be between 0 and 1.");
      bisection_tol_ = value;
    }
    else {
      bisection_tol_ = 1e-3;
    }

    if (quality_function_normalized_) {
      // Set the scaling values to a negative value to indicate that
      // they still have to be determined
      dual_inf_scal_ = -1.;
      primal_inf_scal_ = -1.;
      compl_inf_scal_ = -1.;
    }
    else {
      dual_inf_scal_ = 1.;
      primal_inf_scal_ = 1.;
      compl_inf_scal_ = 1.;
    }

    if (options.GetIntegerValue("dual_alpha_for_y", ivalue, prefix)) {
      dual_alpha_for_y_ = (ivalue != 0);
    }
    else {
      dual_alpha_for_y_ = false;
    }

    // The following line is only here so that
    // IpoptCalculatedQuantities::CalculateSafeSlack and the first
    // output line have something to work with
    IpData().Set_mu(1.);

    return true;
  }

  Number OptProbingMuOracle::CalculateMu()
  {
    DBG_START_METH("OptProbingMuOracle::CalculateMu",
                   dbg_verbosity);

    /////////////////////////////////////
    // Compute the affine scaling step //
    /////////////////////////////////////

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the affine step\n");
    // First get the right hand side
    SmartPtr<const Vector> rhs_aff_x = IpCq().curr_grad_lag_x();
    SmartPtr<const Vector> rhs_aff_s = IpCq().curr_grad_lag_s();
    SmartPtr<const Vector> rhs_aff_c = IpCq().curr_c();
    SmartPtr<const Vector> rhs_aff_d = IpCq().curr_d_minus_s();
    SmartPtr<const Vector> rhs_aff_x_L = IpCq().curr_compl_x_L();
    SmartPtr<const Vector> rhs_aff_x_U = IpCq().curr_compl_x_U();
    SmartPtr<const Vector> rhs_aff_s_L = IpCq().curr_compl_s_L();
    SmartPtr<const Vector> rhs_aff_s_U = IpCq().curr_compl_s_U();

    // Get space for the affine scaling step
    SmartPtr<Vector> step_aff_x = rhs_aff_x->MakeNew();
    SmartPtr<Vector> step_aff_s = rhs_aff_s->MakeNew();
    SmartPtr<Vector> step_aff_y_c = rhs_aff_c->MakeNew();
    SmartPtr<Vector> step_aff_y_d = rhs_aff_d->MakeNew();
    SmartPtr<Vector> step_aff_z_L = rhs_aff_x_L->MakeNew();
    SmartPtr<Vector> step_aff_z_U = rhs_aff_x_U->MakeNew();
    SmartPtr<Vector> step_aff_v_L = rhs_aff_s_L->MakeNew();
    SmartPtr<Vector> step_aff_v_U = rhs_aff_s_U->MakeNew();

    // Now solve the primal-dual system to get the step
    pd_solver_->Solve(-1.0, 0.0,
                      *rhs_aff_x,
                      *rhs_aff_s,
                      *rhs_aff_c,
                      *rhs_aff_d,
                      *rhs_aff_x_L,
                      *rhs_aff_x_U,
                      *rhs_aff_s_L,
                      *rhs_aff_s_U,
                      *step_aff_x,
                      *step_aff_s,
                      *step_aff_y_c,
                      *step_aff_y_d,
                      *step_aff_z_L,
                      *step_aff_z_U,
                      *step_aff_v_L,
                      *step_aff_v_U,
                      false           // want accurate solution here
                      // because we can use it to
                      // compute the overall search
                      // direction
                     );

    DBG_PRINT_VECTOR(2, "step_aff_x", *step_aff_x);
    DBG_PRINT_VECTOR(2, "step_aff_s", *step_aff_s);
    DBG_PRINT_VECTOR(2, "step_aff_y_c", *step_aff_y_c);
    DBG_PRINT_VECTOR(2, "step_aff_y_d", *step_aff_y_d);
    DBG_PRINT_VECTOR(2, "step_aff_z_L", *step_aff_z_L);
    DBG_PRINT_VECTOR(2, "step_aff_z_U", *step_aff_z_U);
    DBG_PRINT_VECTOR(2, "step_aff_v_L", *step_aff_v_L);
    DBG_PRINT_VECTOR(2, "step_aff_v_U", *step_aff_v_U);

    /////////////////////////////////////
    // Compute the pure centering step //
    /////////////////////////////////////

    Number avrg_compl = IpCq().curr_avrg_compl();

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the centering step\n");
    // First get the right hand side
    SmartPtr<Vector> rhs_cen_x = rhs_aff_x->MakeNew();
    SmartPtr<Vector> rhs_cen_s = rhs_aff_s->MakeNew();
    SmartPtr<Vector> rhs_cen_c = rhs_aff_c->MakeNew();
    SmartPtr<Vector> rhs_cen_d = rhs_aff_d->MakeNew();
    SmartPtr<Vector> rhs_cen_x_L = rhs_aff_x_L->MakeNew();
    SmartPtr<Vector> rhs_cen_x_U = rhs_aff_x_U->MakeNew();
    SmartPtr<Vector> rhs_cen_s_L = rhs_aff_s_L->MakeNew();
    SmartPtr<Vector> rhs_cen_s_U = rhs_aff_s_U->MakeNew();

    rhs_cen_x->Set(0.);
    rhs_cen_s->Set(0.);
    rhs_cen_c->Set(0.);
    rhs_cen_d->Set(0.);
    rhs_cen_x_L->Set(avrg_compl);
    rhs_cen_x_U->Set(avrg_compl);
    rhs_cen_s_L->Set(avrg_compl);
    rhs_cen_s_U->Set(avrg_compl);

    // Get space for the affine scaling step
    SmartPtr<Vector> step_cen_x = rhs_aff_x->MakeNew();
    SmartPtr<Vector> step_cen_s = rhs_aff_s->MakeNew();
    SmartPtr<Vector> step_cen_y_c = rhs_aff_c->MakeNew();
    SmartPtr<Vector> step_cen_y_d = rhs_aff_d->MakeNew();
    SmartPtr<Vector> step_cen_z_L = rhs_aff_x_L->MakeNew();
    SmartPtr<Vector> step_cen_z_U = rhs_aff_x_U->MakeNew();
    SmartPtr<Vector> step_cen_v_L = rhs_aff_s_L->MakeNew();
    SmartPtr<Vector> step_cen_v_U = rhs_aff_s_U->MakeNew();

    // Now solve the primal-dual system to get the step
    pd_solver_->Solve(1.0, 0.0,
                      *rhs_cen_x,
                      *rhs_cen_s,
                      *rhs_cen_c,
                      *rhs_cen_d,
                      *rhs_cen_x_L,
                      *rhs_cen_x_U,
                      *rhs_cen_s_L,
                      *rhs_cen_s_U,
                      *step_cen_x,
                      *step_cen_s,
                      *step_cen_y_c,
                      *step_cen_y_d,
                      *step_cen_z_L,
                      *step_cen_z_U,
                      *step_cen_v_L,
                      *step_cen_v_U,
                      false           // want accurate solution here
                      // because we can use it to
                      // compute the overall search
                      // direction
                     );

    DBG_PRINT_VECTOR(2, "step_cen_x", *step_cen_x);
    DBG_PRINT_VECTOR(2, "step_cen_s", *step_cen_s);
    DBG_PRINT_VECTOR(2, "step_cen_y_c", *step_cen_y_c);
    DBG_PRINT_VECTOR(2, "step_cen_y_d", *step_cen_y_d);
    DBG_PRINT_VECTOR(2, "step_cen_z_L", *step_cen_z_L);
    DBG_PRINT_VECTOR(2, "step_cen_z_U", *step_cen_z_U);
    DBG_PRINT_VECTOR(2, "step_cen_v_L", *step_cen_v_L);
    DBG_PRINT_VECTOR(2, "step_cen_v_U", *step_cen_v_U);

    // We now compute the step for the slack variables.  This safes
    // time, because we then don't have to do this any more for each
    // evaluation of the quality function
    SmartPtr<Vector> step_aff_x_L = step_aff_z_L->MakeNew();
    SmartPtr<Vector> step_aff_x_U = step_aff_z_U->MakeNew();
    SmartPtr<Vector> step_aff_s_L = step_aff_v_L->MakeNew();
    SmartPtr<Vector> step_aff_s_U = step_aff_v_U->MakeNew();
    IpNLP().Px_L()->TransMultVector(1., *step_aff_x, 0., *step_aff_x_L);
    IpNLP().Px_U()->TransMultVector(-1., *step_aff_x, 0., *step_aff_x_U);
    IpNLP().Pd_L()->TransMultVector(1., *step_aff_s, 0., *step_aff_s_L);
    IpNLP().Pd_U()->TransMultVector(-1., *step_aff_s, 0., *step_aff_s_U);
    SmartPtr<Vector> step_cen_x_L = step_cen_z_L->MakeNew();
    SmartPtr<Vector> step_cen_x_U = step_cen_z_U->MakeNew();
    SmartPtr<Vector> step_cen_s_L = step_cen_v_L->MakeNew();
    SmartPtr<Vector> step_cen_s_U = step_cen_v_U->MakeNew();
    IpNLP().Px_L()->TransMultVector(1., *step_cen_x, 0., *step_cen_x_L);
    IpNLP().Px_U()->TransMultVector(-1., *step_cen_x, 0., *step_cen_x_U);
    IpNLP().Pd_L()->TransMultVector(1., *step_cen_s, 0., *step_cen_s_L);
    IpNLP().Pd_U()->TransMultVector(-1., *step_cen_s, 0., *step_cen_s_U);

    // If necessary, compute some products with the constraint Jacobian
    SmartPtr<const Vector> jac_cT_times_step_aff_y_c;
    SmartPtr<const Vector> jac_dT_times_step_aff_y_d;
    SmartPtr<const Vector> jac_cT_times_step_cen_y_c;
    SmartPtr<const Vector> jac_dT_times_step_cen_y_d;
    if (quality_function_dual_inf_==2 && dual_alpha_for_y_) {
      jac_cT_times_step_aff_y_c = IpCq().curr_jac_cT_times_vec(*step_aff_y_c);
      jac_dT_times_step_aff_y_d = IpCq().curr_jac_dT_times_vec(*step_aff_y_d);
      jac_cT_times_step_cen_y_c = IpCq().curr_jac_cT_times_vec(*step_cen_y_c);
      jac_dT_times_step_cen_y_d = IpCq().curr_jac_dT_times_vec(*step_cen_y_d);
    }

    Number sigma;
    if (max_bisection_steps_>0) {
      // Now we do an search for the best centering parameter, that
      // gives us the lower value of a quality function, using golden
      // bisection

      Number sigma_up = Min(1., sigma_max_);
      Number sigma_lo = 1e-9/avrg_compl;
      sigma = PerformGoldenBisection(sigma_up, sigma_lo, bisection_tol_,
                                     *step_aff_x_L,
                                     *step_aff_x_U,
                                     *step_aff_s_L,
                                     *step_aff_s_U,
                                     *step_aff_y_c,
                                     *step_aff_y_d,
				     *step_aff_z_L,
                                     *step_aff_z_U,
                                     *step_aff_v_L,
                                     *step_aff_v_U,
                                     *step_cen_x_L,
                                     *step_cen_x_U,
                                     *step_cen_s_L,
                                     *step_cen_s_U,
                                     *step_cen_y_c,
                                     *step_cen_y_d,
                                     *step_cen_z_L,
                                     *step_cen_z_U,
                                     *step_cen_v_L,
                                     *step_cen_v_U,
				     jac_cT_times_step_aff_y_c,
				     jac_dT_times_step_aff_y_d,
				     jac_cT_times_step_cen_y_c,
				     jac_dT_times_step_cen_y_d);

      if (sigma_max_ > 1. && sigma >= 1.-2*bisection_tol_) {
        // It seems that the optimal value might be larger than one.
        sigma_up = sigma_max_;
        sigma_lo = sigma;
        sigma = PerformGoldenBisection(sigma_up, sigma_lo, bisection_tol_,
                                       *step_aff_x_L,
                                       *step_aff_x_U,
                                       *step_aff_s_L,
                                       *step_aff_s_U,
				       *step_aff_y_c,
				       *step_aff_y_d,
                                       *step_aff_z_L,
                                       *step_aff_z_U,
                                       *step_aff_v_L,
                                       *step_aff_v_U,
                                       *step_cen_x_L,
                                       *step_cen_x_U,
                                       *step_cen_s_L,
                                       *step_cen_s_U,
				       *step_cen_y_c,
				       *step_cen_y_d,
                                       *step_cen_z_L,
                                       *step_cen_z_U,
                                       *step_cen_v_L,
                                       *step_cen_v_U,
				       jac_cT_times_step_aff_y_c,
				       jac_dT_times_step_aff_y_d,
				       jac_cT_times_step_cen_y_c,
				       jac_dT_times_step_cen_y_d);
      }

      //#define tracequalityfunction
#ifdef tracequalityfunction
      //DELETEME
      Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                     "\n");
      Number base = 1.2;
      for (Index l=30;l>=(Index)trunc(-(log(avrg_compl)-log(1e-9))/log(base));l--) {
        Number sig = pow(base, l);
        CalculateQualityFunction(sig,
                                 *step_aff_x_L,
                                 *step_aff_x_U,
                                 *step_aff_s_L,
                                 *step_aff_s_U,
				 *step_aff_y_c,
				 *step_aff_y_d,
                                 *step_aff_z_L,
                                 *step_aff_z_U,
                                 *step_aff_v_L,
                                 *step_aff_v_U,
                                 *step_cen_x_L,
                                 *step_cen_x_U,
                                 *step_cen_s_L,
                                 *step_cen_s_U,
				 *step_cen_y_c,
				 *step_cen_y_d,
                                 *step_cen_z_L,
                                 *step_cen_z_U,
                                 *step_cen_v_L,
                                 *step_cen_v_U,
				 jac_cT_times_step_aff_y_c,
				 jac_dT_times_step_aff_y_d,
				 jac_cT_times_step_cen_y_c,
				 jac_dT_times_step_cen_y_d);
      }
#endif

    }
    else {
      Index l;
      Index l_best;
      Number q_best;

      Number base = 1.2;
      l = 20;
      l_best = l;
      sigma = pow(base, l);
      q_best = CalculateQualityFunction(sigma,
                                        *step_aff_x_L,
                                        *step_aff_x_U,
                                        *step_aff_s_L,
                                        *step_aff_s_U,
					*step_aff_y_c,
					*step_aff_y_d,
                                        *step_aff_z_L,
                                        *step_aff_z_U,
                                        *step_aff_v_L,
                                        *step_aff_v_U,
                                        *step_cen_x_L,
                                        *step_cen_x_U,
                                        *step_cen_s_L,
                                        *step_cen_s_U,
					*step_cen_y_c,
					*step_cen_y_d,
                                        *step_cen_z_L,
                                        *step_cen_z_U,
                                        *step_cen_v_L,
                                        *step_cen_v_U,
					jac_cT_times_step_aff_y_c,
					jac_dT_times_step_aff_y_d,
					jac_cT_times_step_cen_y_c,
					jac_dT_times_step_cen_y_d);
      Index l_min = (Index)trunc(-(log(avrg_compl)-log(1e-9))/log(base))-1;
      for (; l>=l_min; l--) {
        sigma = pow(base, l);
        Number q = CalculateQualityFunction(sigma,
                                            *step_aff_x_L,
                                            *step_aff_x_U,
                                            *step_aff_s_L,
                                            *step_aff_s_U,
					    *step_aff_y_c,
					    *step_aff_y_d,
					    *step_aff_z_L,
                                            *step_aff_z_U,
                                            *step_aff_v_L,
                                            *step_aff_v_U,
                                            *step_cen_x_L,
                                            *step_cen_x_U,
                                            *step_cen_s_L,
                                            *step_cen_s_U,
					    *step_cen_y_c,
					    *step_cen_y_d,
                                            *step_cen_z_L,
                                            *step_cen_z_U,
                                            *step_cen_v_L,
                                            *step_cen_v_U,
					    jac_cT_times_step_aff_y_c,
					    jac_dT_times_step_aff_y_d,
					    jac_cT_times_step_cen_y_c,
					    jac_dT_times_step_cen_y_d);
        if (q<=q_best) {
          q_best = q;
          l_best = l;
        }
      }

      sigma = pow(base, l_best);
    }

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Sigma = %e\n", sigma);
    Number mu = sigma*avrg_compl;

    // Now construct the overall search direction here
    step_aff_x->Axpy(sigma, *step_cen_x);
    step_aff_s->Axpy(sigma, *step_cen_s);
    step_aff_y_c->Axpy(sigma, *step_cen_y_c);
    step_aff_y_d->Axpy(sigma, *step_cen_y_d);
    step_aff_z_L->Axpy(sigma, *step_cen_z_L);
    step_aff_z_U->Axpy(sigma, *step_cen_z_U);
    step_aff_v_L->Axpy(sigma, *step_cen_v_L);
    step_aff_v_U->Axpy(sigma, *step_cen_v_U);
    IpData().SetFromPtr_delta_x(step_aff_x);
    IpData().SetFromPtr_delta_s(step_aff_s);
    IpData().SetFromPtr_delta_y_c(step_aff_y_c);
    IpData().SetFromPtr_delta_y_d(step_aff_y_d);
    IpData().SetFromPtr_delta_z_L(step_aff_z_L);
    IpData().SetFromPtr_delta_z_U(step_aff_z_U);
    IpData().SetFromPtr_delta_v_L(step_aff_v_L);
    IpData().SetFromPtr_delta_v_U(step_aff_v_U);
    IpData().SetHaveDeltas(true);

    // DELETEME
    char ssigma[40];
    sprintf(ssigma, " sigma=%8.2e", sigma);
    IpData().Append_info_string(ssigma);
    sprintf(ssigma, " xi=%8.2e ", IpCq().curr_centrality_measure());
    IpData().Append_info_string(ssigma);
    if (sigma>1.) {
      IpData().Append_info_string("LARGESIGMA");
    }

    return mu;
  }

  Number OptProbingMuOracle::CalculateQualityFunction
  (Number sigma,
   const Vector& step_aff_x_L,
   const Vector& step_aff_x_U,
   const Vector& step_aff_s_L,
   const Vector& step_aff_s_U,
   const Vector& step_aff_y_c,
   const Vector& step_aff_y_d,
   const Vector& step_aff_z_L,
   const Vector& step_aff_z_U,
   const Vector& step_aff_v_L,
   const Vector& step_aff_v_U,
   const Vector& step_cen_x_L,
   const Vector& step_cen_x_U,
   const Vector& step_cen_s_L,
   const Vector& step_cen_s_U,
   const Vector& step_cen_y_c,
   const Vector& step_cen_y_d,
   const Vector& step_cen_z_L,
   const Vector& step_cen_z_U,
   const Vector& step_cen_v_L,
   const Vector& step_cen_v_U,
   SmartPtr<const Vector> jac_cT_times_step_aff_y_c,
   SmartPtr<const Vector> jac_dT_times_step_aff_y_d,
   SmartPtr<const Vector> jac_cT_times_step_cen_y_c,
   SmartPtr<const Vector> jac_dT_times_step_cen_y_d
  )
  {
    DBG_START_METH("OptProbingMuOracle::CalculateQualityFunction",
                   dbg_verbosity);

    Index n_dual = IpData().curr_x()->Dim() + IpData().curr_s()->Dim();
    Index n_pri = IpData().curr_y_c()->Dim() + IpData().curr_y_d()->Dim();
    Index n_comp = IpData().curr_z_L()->Dim() + IpData().curr_z_U()->Dim() +
                   IpData().curr_v_L()->Dim() + IpData().curr_v_U()->Dim();

    // The scaling values have not yet been determined, compute them now
    if (dual_inf_scal_ < 0.) {
      // DELETEME
      assert(false && "Scaling in quality function not supported.");
      DBG_ASSERT(primal_inf_scal_ < 0.);
      DBG_ASSERT(compl_inf_scal_ < 0.);

      switch (quality_function_norm_) {
      case 1:
        dual_inf_scal_ = Max(1., IpCq().curr_grad_lag_x()->Asum() +
                             IpCq().curr_grad_lag_s()->Asum());

        primal_inf_scal_ = Max(1., IpCq().curr_c()->Asum() +
                               IpCq().curr_d_minus_s()->Asum());

        compl_inf_scal_ = Max(1., IpCq().curr_compl_x_L()->Asum() +
                              IpCq().curr_compl_x_U()->Asum() +
                              IpCq().curr_compl_s_L()->Asum() +
                              IpCq().curr_compl_s_U()->Asum());

        dual_inf_scal_ /= n_dual;
        if (n_pri>0) {
          primal_inf_scal_ /= n_pri;
        }
        DBG_ASSERT(n_comp>0);
        compl_inf_scal_ /= n_comp;
	break;
      case 2:
        dual_inf_scal_ = Max(1., pow(IpCq().curr_grad_lag_x()->Nrm2(), 2) +
                             pow(IpCq().curr_grad_lag_s()->Nrm2(), 2));

        primal_inf_scal_ = Max(1., pow(IpCq().curr_c()->Nrm2(), 2) +
                               pow(IpCq().curr_d_minus_s()->Nrm2(), 2));

        compl_inf_scal_ = Max(1., pow(IpCq().curr_compl_x_L()->Nrm2(), 2) +
                              pow(IpCq().curr_compl_x_U()->Nrm2(), 2) +
                              pow(IpCq().curr_compl_s_L()->Nrm2(), 2) +
                              pow(IpCq().curr_compl_s_U()->Nrm2(), 2));

	dual_inf_scal_ /= n_dual;
	if (n_pri>0) {
	  primal_inf_scal_ /= n_pri;
	}
	DBG_ASSERT(n_comp>0);
	compl_inf_scal_ /= n_comp;
	break;
      case 3:
        dual_inf_scal_ = Max(1., IpCq().curr_grad_lag_x()->Amax(),
                             IpCq().curr_grad_lag_s()->Amax());

        primal_inf_scal_ = Max(1., IpCq().curr_c()->Amax(),
                               IpCq().curr_d_minus_s()->Amax());

        compl_inf_scal_ = Max(1., Max(IpCq().curr_compl_x_L()->Amax(),
				      IpCq().curr_compl_x_U()->Amax(),
				      IpCq().curr_compl_s_L()->Amax(),
				      IpCq().curr_compl_s_U()->Amax()));
      }
    }

    // First compute the corresponding search direction
    SmartPtr<Vector> step_x_L = step_aff_x_L.MakeNew();
    SmartPtr<Vector> step_x_U = step_aff_x_U.MakeNew();
    SmartPtr<Vector> step_s_L = step_aff_s_L.MakeNew();
    SmartPtr<Vector> step_s_U = step_aff_s_U.MakeNew();
    SmartPtr<Vector> step_z_L = step_aff_z_L.MakeNew();
    SmartPtr<Vector> step_z_U = step_aff_z_U.MakeNew();
    SmartPtr<Vector> step_v_L = step_aff_v_L.MakeNew();
    SmartPtr<Vector> step_v_U = step_aff_v_U.MakeNew();

    step_x_L->Copy(step_aff_x_L);
    step_x_U->Copy(step_aff_x_U);
    step_s_L->Copy(step_aff_s_L);
    step_s_U->Copy(step_aff_s_U);
    step_z_L->Copy(step_aff_z_L);
    step_z_U->Copy(step_aff_z_U);
    step_v_L->Copy(step_aff_v_L);
    step_v_U->Copy(step_aff_v_U);

    step_x_L->Axpy(sigma, step_cen_x_L);
    step_x_U->Axpy(sigma, step_cen_x_U);
    step_s_L->Axpy(sigma, step_cen_s_L);
    step_s_U->Axpy(sigma, step_cen_s_U);
    step_z_L->Axpy(sigma, step_cen_z_L);
    step_z_U->Axpy(sigma, step_cen_z_U);
    step_v_L->Axpy(sigma, step_cen_v_L);
    step_v_U->Axpy(sigma, step_cen_v_U);

    // Compute the fraction-to-the-boundary step sizes
    // ToDo make sure we use the correct tau
    Number tau = 0.99;
    Number alpha_primal = IpCq().slack_frac_to_the_bound(tau,
                          *step_x_L,
                          *step_x_U,
                          *step_s_L,
                          *step_s_U);

    Number alpha_dual = IpCq().dual_frac_to_the_bound(tau,
                        *step_z_L,
                        *step_z_U,
                        *step_v_L,
                        *step_v_U);

    if (false) {
      if (alpha_dual < alpha_primal) {
        alpha_primal = alpha_dual;
      }
      else {
        alpha_dual = alpha_primal;
      }
    }

    // Additional reduction factor for the step size to ensure that
    // they are in the -infinity neighborhood
    // THIS DOESN'T WORK (YET?)
    Number beta = 1.;
    bool found_beta = false;
    Number xi; // centrality measure

    SmartPtr<Vector> slack_x_L = step_aff_z_L.MakeNew();
    SmartPtr<Vector> slack_x_U = step_aff_z_U.MakeNew();
    SmartPtr<Vector> slack_s_L = step_aff_v_L.MakeNew();
    SmartPtr<Vector> slack_s_U = step_aff_v_U.MakeNew();
    SmartPtr<Vector> z_L = step_aff_z_L.MakeNew();
    SmartPtr<Vector> z_U = step_aff_z_U.MakeNew();
    SmartPtr<Vector> v_L = step_aff_v_L.MakeNew();
    SmartPtr<Vector> v_U = step_aff_v_U.MakeNew();

    while (!found_beta) {

      slack_x_L->Copy(*IpCq().curr_slack_x_L());
      slack_x_U->Copy(*IpCq().curr_slack_x_U());
      slack_s_L->Copy(*IpCq().curr_slack_s_L());
      slack_s_U->Copy(*IpCq().curr_slack_s_U());
      slack_x_L->Axpy(alpha_primal, *step_x_L);
      slack_x_U->Axpy(alpha_primal, *step_x_U);
      slack_s_L->Axpy(alpha_primal, *step_s_L);
      slack_s_U->Axpy(alpha_primal, *step_s_U);

      z_L->Copy(*IpData().curr_z_L());
      z_U->Copy(*IpData().curr_z_U());
      v_L->Copy(*IpData().curr_v_L());
      v_U->Copy(*IpData().curr_v_U());
      z_L->Axpy(alpha_dual, *step_z_L);
      z_U->Axpy(alpha_dual, *step_z_U);
      v_L->Axpy(alpha_dual, *step_v_L);
      v_U->Axpy(alpha_dual, *step_v_U);

      slack_x_L->ElementWiseMultiply(*z_L);
      slack_x_U->ElementWiseMultiply(*z_U);
      slack_s_L->ElementWiseMultiply(*v_L);
      slack_s_U->ElementWiseMultiply(*v_U);

      DBG_PRINT_VECTOR(2, "compl_x_L", *slack_x_L);
      DBG_PRINT_VECTOR(2, "compl_x_U", *slack_x_U);
      DBG_PRINT_VECTOR(2, "compl_s_L", *slack_s_L);
      DBG_PRINT_VECTOR(2, "compl_s_U", *slack_s_U);

      xi = IpCq().CalcCentralityMeasure(*slack_x_L, *slack_x_U,
                                        *slack_s_L, *slack_s_U);

      // in order to make this work, we need to somehow tell the line
      // search what beta is
      Number xi_min_ = 0.;
      if (xi_min_>0.) {
        DBG_ASSERT(IpCq().curr_centrality_measure()>xi_min_);
        Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                       "beta = %8.2e alpha_p = %8.2e alpha_d = %8.2e xi = %8.2e\n", beta, alpha_primal, alpha_dual, xi);
        if (xi>=xi_min_) {
          found_beta = true;
        }
        else {
          Number redfact = 0.99;
          beta *= redfact;
          alpha_primal *= redfact;
          alpha_dual *= redfact;
        }
      }
      else {
        found_beta = true;
      }

    }

    if (beta<1.)
      IpData().Append_info_string("b");

    Number dual_inf;
    Number primal_inf;
    Number compl_inf;

    SmartPtr<Vector> dual_inf_x;
    SmartPtr<Vector> dual_inf_s;

    if (quality_function_dual_inf_==2) {
      dual_inf_x = IpCq().curr_grad_lag_x()->MakeNew();
      dual_inf_s = IpCq().curr_grad_lag_s()->MakeNew();
      dual_inf_x->Copy(*IpCq().curr_grad_lag_x());
      dual_inf_s->Copy(*IpCq().curr_grad_lag_s());
      dual_inf_x->Scal(1.-alpha_primal);
      dual_inf_s->Scal(1.-alpha_primal);
      IpNLP().Px_L()->MultVector(alpha_primal-alpha_dual,
				 *step_z_L, 1., *dual_inf_x);
      IpNLP().Px_U()->MultVector(-alpha_primal+alpha_dual,
				 *step_z_U, 1., *dual_inf_x);
      IpNLP().Pd_L()->MultVector(alpha_primal-alpha_dual,
				 *step_v_L, 1., *dual_inf_s);
      IpNLP().Pd_U()->MultVector(-alpha_primal+alpha_dual,
				 *step_v_U, 1., *dual_inf_s);
      if (dual_alpha_for_y_) {
	dual_inf_x->Axpy(alpha_dual-alpha_primal, *jac_cT_times_step_aff_y_c);  
	dual_inf_x->Axpy(sigma*(alpha_dual-alpha_primal), *jac_cT_times_step_cen_y_c);
	dual_inf_x->Axpy(alpha_dual-alpha_primal, *jac_dT_times_step_aff_y_d);  
	dual_inf_x->Axpy(sigma*(alpha_dual-alpha_primal), *jac_dT_times_step_cen_y_d);
	dual_inf_s->Axpy(alpha_primal-alpha_dual, step_aff_y_d);
	dual_inf_s->Axpy(sigma*(alpha_primal-alpha_dual), step_cen_y_d);
      }
    }

    switch (quality_function_norm_) {
    case 1:
      if (quality_function_dual_inf_==2) {
	dual_inf = dual_inf_x->Asum() + dual_inf_s->Asum();
      }
      else {
	dual_inf = (1.-alpha_dual)*(IpCq().curr_grad_lag_x()->Asum() +
				    IpCq().curr_grad_lag_s()->Asum());
      }

      primal_inf = (1.-alpha_primal)*(IpCq().curr_c()->Asum() +
                                      IpCq().curr_d_minus_s()->Asum());

      compl_inf = slack_x_L->Asum() + slack_x_U->Asum() +
                  slack_s_L->Asum() + slack_s_U->Asum();

      dual_inf /= n_dual;
      if (n_pri>0) {
        primal_inf /= n_pri;
      }
      DBG_ASSERT(n_comp>0);
      compl_inf /= n_comp;
      break;
    case 2:
      if (quality_function_dual_inf_==2) {
	dual_inf = pow(dual_inf_x->Nrm2(), 2) + pow(dual_inf_s->Nrm2(), 2);
      }
      else {
	dual_inf =
	  pow(1.-alpha_dual, 2)*(pow(IpCq().curr_grad_lag_x()->Nrm2(), 2) +
				 pow(IpCq().curr_grad_lag_s()->Nrm2(), 2));
      }
      primal_inf =
        pow(1.-alpha_primal, 2)*(pow(IpCq().curr_c()->Nrm2(), 2) +
                                 pow(IpCq().curr_d_minus_s()->Nrm2(), 2));
      compl_inf =
        pow(slack_x_L->Nrm2(), 2) + pow(slack_x_U->Nrm2(), 2) +
        pow(slack_s_L->Nrm2(), 2) + pow(slack_s_U->Nrm2(), 2);

      dual_inf /= n_dual;
      if (n_pri>0) {
        primal_inf /= n_pri;
      }
      DBG_ASSERT(n_comp>0);
      compl_inf /= n_comp;
      break;
    case 3:
      if (quality_function_dual_inf_==2) {
	dual_inf = Max(dual_inf_x->Amax(), dual_inf_s->Amax());
      }
      else {
	dual_inf =
	  (1.-alpha_dual)*Max(IpCq().curr_grad_lag_x()->Amax(),
			      IpCq().curr_grad_lag_s()->Amax());
      }
      primal_inf =
        (1.-alpha_primal, 2)*Max(IpCq().curr_c()->Amax(),
                                 IpCq().curr_d_minus_s()->Amax());
      compl_inf =
        Max(slack_x_L->Amax(), slack_x_U->Amax(),
	    slack_s_L->Amax(), slack_s_U->Amax());
      break;
    case 4:
      if (quality_function_dual_inf_==2) {
	dual_inf = sqrt(pow(dual_inf_x->Nrm2(), 2) + pow(dual_inf_s->Nrm2(), 2));
      }
      else {
	dual_inf =
	 (1.-alpha_dual)*sqrt(pow(IpCq().curr_grad_lag_x()->Nrm2(), 2) +
			      pow(IpCq().curr_grad_lag_s()->Nrm2(), 2));
      }
      primal_inf =
        (1.-alpha_primal, 2)*sqrt(pow(IpCq().curr_c()->Nrm2(), 2) +
				  pow(IpCq().curr_d_minus_s()->Nrm2(), 2));
      compl_inf =
        sqrt(pow(slack_x_L->Nrm2(), 2) + pow(slack_x_U->Nrm2(), 2) +
	     pow(slack_s_L->Nrm2(), 2) + pow(slack_s_U->Nrm2(), 2));

      dual_inf /= sqrt((Number)n_dual);
      if (n_pri>0) {
        primal_inf /= sqrt((Number)n_pri);
      }
      DBG_ASSERT(n_comp>0);
      compl_inf /= sqrt((Number)n_comp);
      break;
    default:
      DBG_ASSERT("Unknown value for quality_function_norm_");
    }

    // Scale the quantities
    dual_inf /= dual_inf_scal_;
    primal_inf /= primal_inf_scal_;
    compl_inf /= compl_inf_scal_;

    Number quality_function = dual_inf + primal_inf + compl_inf;

    switch (quality_function_centrality_) {
      case 0:
      //Nothing
      break;
      case 1:
      quality_function -= compl_inf*log(xi);
      break;
      case 2:
      quality_function += compl_inf/xi;
      case 3:
      quality_function += compl_inf/pow(xi,3);
      break;
      default:
      DBG_ASSERT("Unknown value for quality_function_centrality_");
    }

    Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                   "sigma = %8.2e d_inf = %18.12e p_inf = %18.12e cmpl = %18.12e q = %18.12e a_pri = %8.2e a_dual = %8.2e xi = %8.2e\n", sigma, dual_inf, primal_inf, compl_inf, quality_function, alpha_primal, alpha_dual, xi);


    return quality_function;
    //return compl_inf;
  }

  Number
  OptProbingMuOracle::PerformGoldenBisection
  (Number sigma_up,
   Number sigma_lo,
   Number tol,
   const Vector& step_aff_x_L,
   const Vector& step_aff_x_U,
   const Vector& step_aff_s_L,
   const Vector& step_aff_s_U,
   const Vector& step_aff_y_c,
   const Vector& step_aff_y_d,
   const Vector& step_aff_z_L,
   const Vector& step_aff_z_U,
   const Vector& step_aff_v_L,
   const Vector& step_aff_v_U,
   const Vector& step_cen_x_L,
   const Vector& step_cen_x_U,
   const Vector& step_cen_s_L,
   const Vector& step_cen_s_U,
   const Vector& step_cen_y_c,
   const Vector& step_cen_y_d,
   const Vector& step_cen_z_L,
   const Vector& step_cen_z_U,
   const Vector& step_cen_v_L,
   const Vector& step_cen_v_U,
   SmartPtr<const Vector> jac_cT_times_step_aff_y_c,
   SmartPtr<const Vector> jac_dT_times_step_aff_y_d,
   SmartPtr<const Vector> jac_cT_times_step_cen_y_c,
   SmartPtr<const Vector> jac_dT_times_step_cen_y_d
  )
  {
    Number sigma;
    Number sigma_up_in = sigma_up;
    Number sigma_lo_in = sigma_lo;
    Number gfac = (3.-sqrt(5.))/2.;
    Number sigma_mid1 = sigma_lo + gfac*(sigma_up-sigma_lo);
    Number sigma_mid2 = sigma_lo + (1.-gfac)*(sigma_up-sigma_lo);

    Number qmid1 = CalculateQualityFunction(sigma_mid1,
                                            step_aff_x_L,
                                            step_aff_x_U,
                                            step_aff_s_L,
                                            step_aff_s_U,
                                            step_aff_y_c,
                                            step_aff_y_d,
                                            step_aff_z_L,
                                            step_aff_z_U,
                                            step_aff_v_L,
                                            step_aff_v_U,
                                            step_cen_x_L,
                                            step_cen_x_U,
                                            step_cen_s_L,
                                            step_cen_s_U,
                                            step_cen_y_c,
                                            step_cen_y_d,
                                            step_cen_z_L,
                                            step_cen_z_U,
                                            step_cen_v_L,
                                            step_cen_v_U,
					    jac_cT_times_step_aff_y_c,
					    jac_dT_times_step_aff_y_d,
					    jac_cT_times_step_cen_y_c,
					    jac_dT_times_step_cen_y_d);
    Number qmid2 = CalculateQualityFunction(sigma_mid2,
                                            step_aff_x_L,
                                            step_aff_x_U,
                                            step_aff_s_L,
                                            step_aff_s_U,
                                            step_aff_y_c,
                                            step_aff_y_d,
                                            step_aff_z_L,
                                            step_aff_z_U,
                                            step_aff_v_L,
                                            step_aff_v_U,
                                            step_cen_x_L,
                                            step_cen_x_U,
                                            step_cen_s_L,
                                            step_cen_s_U,
                                            step_cen_y_c,
                                            step_cen_y_d,
                                            step_cen_z_L,
                                            step_cen_z_U,
                                            step_cen_v_L,
                                            step_cen_v_U,
					    jac_cT_times_step_aff_y_c,
					    jac_dT_times_step_aff_y_d,
					    jac_cT_times_step_cen_y_c,
					    jac_dT_times_step_cen_y_d);

    Index nbisections = 0;
    while ((sigma_up-sigma_lo)>=tol*sigma_up && nbisections<max_bisection_steps_) {
      //      printf("lo = %e mid1 = %e mid2 = %e up = %e\n",sigma_lo,sigma_mid1,sigma_mid2,sigma_up);
      nbisections++;
      if (qmid1 > qmid2) {
        sigma_lo = sigma_mid1;
        sigma_mid1 = sigma_mid2;
        qmid1 = qmid2;
        sigma_mid2 = sigma_lo + (1.-gfac)*(sigma_up-sigma_lo);
        qmid2 = CalculateQualityFunction(sigma_mid2,
                                         step_aff_x_L,
                                         step_aff_x_U,
                                         step_aff_s_L,
                                         step_aff_s_U,
					 step_aff_y_c,
					 step_aff_y_d,
                                         step_aff_z_L,
                                         step_aff_z_U,
                                         step_aff_v_L,
                                         step_aff_v_U,
                                         step_cen_x_L,
                                         step_cen_x_U,
                                         step_cen_s_L,
                                         step_cen_s_U,
					 step_cen_y_c,
					 step_cen_y_d,
                                         step_cen_z_L,
                                         step_cen_z_U,
                                         step_cen_v_L,
                                         step_cen_v_U,
					 jac_cT_times_step_aff_y_c,
					 jac_dT_times_step_aff_y_d,
					 jac_cT_times_step_cen_y_c,
					 jac_dT_times_step_cen_y_d);
      }
      else {
        sigma_up = sigma_mid2;
        sigma_mid2 = sigma_mid1;
        qmid2 = qmid1;
        sigma_mid1 = sigma_lo + gfac*(sigma_up-sigma_lo);
        qmid1 = CalculateQualityFunction(sigma_mid1,
                                         step_aff_x_L,
                                         step_aff_x_U,
                                         step_aff_s_L,
                                         step_aff_s_U,
					 step_aff_y_c,
					 step_aff_y_d,
                                         step_aff_z_L,
                                         step_aff_z_U,
                                         step_aff_v_L,
                                         step_aff_v_U,
                                         step_cen_x_L,
                                         step_cen_x_U,
                                         step_cen_s_L,
                                         step_cen_s_U,
					 step_cen_y_c,
					 step_cen_y_d,
                                         step_cen_z_L,
                                         step_cen_z_U,
                                         step_cen_v_L,
                                         step_cen_v_U,
					 jac_cT_times_step_aff_y_c,
					 jac_dT_times_step_aff_y_d,
					 jac_cT_times_step_cen_y_c,
					 jac_dT_times_step_cen_y_d);
      }
    }

    Number q;
    if (qmid1 < qmid2) {
      sigma = sigma_mid1;
      q = qmid1;
    }
    else {
      sigma = sigma_mid2;
      q = qmid2;
    }
    if (sigma_up == sigma_up_in) {
      Number qtmp = CalculateQualityFunction(sigma_up,
                                             step_aff_x_L,
                                             step_aff_x_U,
                                             step_aff_s_L,
                                             step_aff_s_U,
					     step_aff_y_c,
					     step_aff_y_d,
                                             step_aff_z_L,
                                             step_aff_z_U,
                                             step_aff_v_L,
                                             step_aff_v_U,
                                             step_cen_x_L,
                                             step_cen_x_U,
                                             step_cen_s_L,
                                             step_cen_s_U,
					     step_cen_y_c,
					     step_cen_y_d,
                                             step_cen_z_L,
                                             step_cen_z_U,
                                             step_cen_v_L,
                                             step_cen_v_U,
					     jac_cT_times_step_aff_y_c,
					     jac_dT_times_step_aff_y_d,
					     jac_cT_times_step_cen_y_c,
					     jac_dT_times_step_cen_y_d);
      if (qtmp < q) {
        sigma = sigma_up;
        q = qtmp;
      }
    }
    else if (sigma_lo == sigma_lo_in) {
      Number qtmp = CalculateQualityFunction(sigma_lo,
                                             step_aff_x_L,
                                             step_aff_x_U,
                                             step_aff_s_L,
                                             step_aff_s_U,
					     step_aff_y_c,
					     step_aff_y_d,
                                             step_aff_z_L,
                                             step_aff_z_U,
                                             step_aff_v_L,
                                             step_aff_v_U,
                                             step_cen_x_L,
                                             step_cen_x_U,
                                             step_cen_s_L,
                                             step_cen_s_U,
					     step_cen_y_c,
					     step_cen_y_d,
                                             step_cen_z_L,
                                             step_cen_z_U,
                                             step_cen_v_L,
                                             step_cen_v_U,
					     jac_cT_times_step_aff_y_c,
					     jac_dT_times_step_aff_y_d,
					     jac_cT_times_step_cen_y_c,
					     jac_dT_times_step_cen_y_d);
      if (qtmp < q) {
        sigma = sigma_lo;
        q = qtmp;
      }
    }

    return sigma;
  }


} // namespace Ipopt
