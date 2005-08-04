// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2004-11-12

#include "IpQualityFunctionMuOracle.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  DefineIpoptType(QualityFunctionMuOracle);

  QualityFunctionMuOracle::QualityFunctionMuOracle(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      MuOracle(),
      pd_solver_(pd_solver),

      tmp_step_x_L_(NULL),
      tmp_step_x_U_(NULL),
      tmp_step_s_L_(NULL),
      tmp_step_s_U_(NULL),
      tmp_step_z_L_(NULL),
      tmp_step_z_U_(NULL),
      tmp_step_v_L_(NULL),
      tmp_step_v_U_(NULL),

      tmp_slack_x_L_(NULL),
      tmp_slack_x_U_(NULL),
      tmp_slack_s_L_(NULL),
      tmp_slack_s_U_(NULL),
      tmp_z_L_(NULL),
      tmp_z_U_(NULL),
      tmp_v_L_(NULL),
      tmp_v_U_(NULL)
  {
    DBG_ASSERT(IsValid(pd_solver_));
  }

  QualityFunctionMuOracle::~QualityFunctionMuOracle()
  {}

  void QualityFunctionMuOracle::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "sigma_max",
      "Maximal value of centering parameter.",
      0.0, true, 1e2);
    roptions->AddStringOption4(
      "quality_function_norm_type",
      "Norm used for components of the quality function.",
      "2-norm-squared",
      "1-norm", "use the 1-norm (abs sum)",
      "2-norm-squared", "use the 2-norm squared (sum of squares)",
      "max-norm", "use the infinity norm (max)",
      "2-norm", "use 2-norm");
    roptions->AddStringOption4(
      "quality_function_centrality",
      "Determines whether a penalty term for centrality is included quality function.",
      "none",
      "none", "no penalty term is added",
      "log", "complementarity * the log of the centrality measure",
      "reciprocal", "complementarity * the reciprocal of the centrality measure",
      "cubed-reciprocal", "complementarity * the reciprocal of the centrality measure cubed",
      "This determines whether a term penalizing deviation from centrality "
      "with respect to complementarity is added the quality function.  The "
      "complementarity measure here is the xi in the Loqo update rule.");
    roptions->AddStringOption2(
      "quality_function_balancing_term",
      "Determines whether a balancing term for centrality is included in quality function.",
      "none",
      "none", "no balancing term is added",
      "cubic", "Max(0,Max(dual_ing,primal_inf)-compl)^3",
      "This determines whether a term penalizing stuations there the "
      "complementality is much smaller than dual and primal "
      "infeasibilities is added to the quality function.");
    roptions->AddLowerBoundedIntegerOption(
      "max_bisection_steps",
      "Maximal number of search steps during direct search procedure determining optimal centering parameter.",
      0, 4);
    roptions->AddBoundedNumberOption(
      "bisection_tol",
      "Tolerance for the bisection search procedure determining optimal centering parameter.",
      0.0, true, 1.0, true,
      1e-3);
  }


  bool QualityFunctionMuOracle::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    Index enum_int;

    options.GetNumericValue("sigma_max", sigma_max_, prefix);

    options.GetEnumValue("quality_function_norm_type", enum_int, prefix);
    quality_function_norm_ = NormEnum(enum_int);
    options.GetEnumValue("quality_function_centrality", enum_int, prefix);
    quality_function_centrality_ = CentralityEnum(enum_int);
    options.GetEnumValue("quality_function_balancing_term", enum_int, prefix);
    quality_function_balancing_term_ = BalancingTermEnum(enum_int);
    options.GetIntegerValue("max_bisection_steps", max_bisection_steps_, prefix);
    options.GetNumericValue("bisection_tol", bisection_tol_, prefix);

    return true;
  }

  Number QualityFunctionMuOracle::CalculateMu()
  {
    DBG_START_METH("QualityFunctionMuOracle::CalculateMu",
                   dbg_verbosity);

    ///////////////////////////////////////////////////////////////////////////
    // Reserve memory for temporary vectors used in CalculateQualityFunction //
    ///////////////////////////////////////////////////////////////////////////

    tmp_step_x_L_ = IpNLP().x_L()->MakeNew();
    tmp_step_x_U_ = IpNLP().x_U()->MakeNew();
    tmp_step_s_L_ = IpNLP().d_L()->MakeNew();
    tmp_step_s_U_ = IpNLP().d_U()->MakeNew();
    tmp_step_z_L_ = IpNLP().x_L()->MakeNew();
    tmp_step_z_U_ = IpNLP().x_U()->MakeNew();
    tmp_step_v_L_ = IpNLP().d_L()->MakeNew();
    tmp_step_v_U_ = IpNLP().d_U()->MakeNew();

    tmp_slack_x_L_ = IpNLP().x_L()->MakeNew();
    tmp_slack_x_U_ = IpNLP().x_U()->MakeNew();
    tmp_slack_s_L_ = IpNLP().d_L()->MakeNew();
    tmp_slack_s_U_ = IpNLP().d_U()->MakeNew();
    tmp_z_L_ = IpNLP().x_L()->MakeNew();
    tmp_z_U_ = IpNLP().x_U()->MakeNew();
    tmp_v_L_ = IpNLP().d_L()->MakeNew();
    tmp_v_U_ = IpNLP().d_U()->MakeNew();

    /////////////////////////////////////
    // Compute the affine scaling step //
    /////////////////////////////////////

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the affine step\n");
    // First get the right hand side
    SmartPtr<IteratesVector> rhs_aff = IpData().curr()->MakeNewIteratesVector(false);
    rhs_aff->Set_x(*IpCq().curr_grad_lag_x());
    rhs_aff->Set_s(*IpCq().curr_grad_lag_s());
    rhs_aff->Set_y_c(*IpCq().curr_c());
    rhs_aff->Set_y_d(*IpCq().curr_d_minus_s());
    rhs_aff->Set_z_L(*IpCq().curr_compl_x_L());
    rhs_aff->Set_z_U(*IpCq().curr_compl_x_U());
    rhs_aff->Set_v_L(*IpCq().curr_compl_s_L());
    rhs_aff->Set_v_U(*IpCq().curr_compl_s_U());

    // Get space for the affine scaling step
    SmartPtr<IteratesVector> step_aff = IpData().curr()->MakeNewIteratesVector(true);

    // Now solve the primal-dual system to get the step
    pd_solver_->Solve(-1.0, 0.0,
                      *rhs_aff,
                      *step_aff,
                      false           // want accurate solution here
                      // because we can use it to
                      // compute the overall search
                      // direction
                     );

    DBG_PRINT_VECTOR(2, "step_aff", *step_aff);

    /////////////////////////////////////
    // Compute the pure centering step //
    /////////////////////////////////////

    Number avrg_compl = IpCq().curr_avrg_compl();

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the centering step\n");
    // First get the right hand side
    SmartPtr<IteratesVector> rhs_cen = IpData().curr()->MakeNewIteratesVector(true);
    rhs_cen->x_NonConst()->AddOneVector(-avrg_compl,
                                        *IpCq().grad_kappa_times_damping_x(),
                                        0.);
    rhs_cen->s_NonConst()->AddOneVector(-avrg_compl,
                                        *IpCq().grad_kappa_times_damping_s(),
                                        0.);

    rhs_cen->y_c_NonConst()->Set(0.);
    rhs_cen->y_d_NonConst()->Set(0.);
    rhs_cen->z_L_NonConst()->Set(avrg_compl);
    rhs_cen->z_U_NonConst()->Set(avrg_compl);
    rhs_cen->v_L_NonConst()->Set(avrg_compl);
    rhs_cen->v_U_NonConst()->Set(avrg_compl);

    // Get space for the centering step
    SmartPtr<IteratesVector> step_cen = IpData().curr()->MakeNewIteratesVector(true);

    // Now solve the primal-dual system to get the step
    pd_solver_->Solve(1.0, 0.0,
                      *rhs_cen,
                      *step_cen,
                      false           // want accurate solution here
                      // because we can use it to
                      // compute the overall search
                      // direction
                     );

    DBG_PRINT_VECTOR(2, "step_cen", *step_cen);

    // We now compute the step for the slack variables.  This safes
    // time, because we then don't have to do this any more for each
    // evaluation of the quality function
    SmartPtr<Vector> step_aff_x_L = step_aff->z_L()->MakeNew();
    SmartPtr<Vector> step_aff_x_U = step_aff->z_U()->MakeNew();
    SmartPtr<Vector> step_aff_s_L = step_aff->v_L()->MakeNew();
    SmartPtr<Vector> step_aff_s_U = step_aff->v_U()->MakeNew();
    IpNLP().Px_L()->TransMultVector(1., *step_aff->x(), 0., *step_aff_x_L);
    IpNLP().Px_U()->TransMultVector(-1., *step_aff->x(), 0., *step_aff_x_U);
    IpNLP().Pd_L()->TransMultVector(1., *step_aff->s(), 0., *step_aff_s_L);
    IpNLP().Pd_U()->TransMultVector(-1., *step_aff->s(), 0., *step_aff_s_U);
    SmartPtr<Vector> step_cen_x_L = step_cen->z_L()->MakeNew();
    SmartPtr<Vector> step_cen_x_U = step_cen->z_U()->MakeNew();
    SmartPtr<Vector> step_cen_s_L = step_cen->v_L()->MakeNew();
    SmartPtr<Vector> step_cen_s_U = step_cen->v_U()->MakeNew();
    IpNLP().Px_L()->TransMultVector(1., *step_cen->x(), 0., *step_cen_x_L);
    IpNLP().Px_U()->TransMultVector(-1., *step_cen->x(), 0., *step_cen_x_U);
    IpNLP().Pd_L()->TransMultVector(1., *step_cen->s(), 0., *step_cen_s_L);
    IpNLP().Pd_U()->TransMultVector(-1., *step_cen->s(), 0., *step_cen_s_U);

    Number sigma;
    if (max_bisection_steps_>0) {
      // Now we do an search for the best centering parameter, that
      // gives us the lower value of a quality function, using golden
      // bisection

      Number sigma_up = Min(1., sigma_max_);
      Number sigma_lo = 1e-9/avrg_compl;
      sigma = PerformGoldenBisection(sigma_up, sigma_lo, bisection_tol_,
                                     //sigma = PerformGoldenBisectionLog(sigma_up, sigma_lo, bisection_tol_,
                                     *step_aff_x_L,
                                     *step_aff_x_U,
                                     *step_aff_s_L,
                                     *step_aff_s_U,
                                     *step_aff->y_c(),
                                     *step_aff->y_d(),
                                     *step_aff->z_L(),
                                     *step_aff->z_U(),
                                     *step_aff->v_L(),
                                     *step_aff->v_U(),
                                     *step_cen_x_L,
                                     *step_cen_x_U,
                                     *step_cen_s_L,
                                     *step_cen_s_U,
                                     *step_cen->y_c(),
                                     *step_cen->y_d(),
                                     *step_cen->z_L(),
                                     *step_cen->z_U(),
                                     *step_cen->v_L(),
                                     *step_cen->v_U());

      if (sigma_max_ > 1. && sigma >= 1.-2*bisection_tol_) {
        // It seems that the optimal value might be larger than one.
        sigma_up = sigma_max_;
        sigma_lo = sigma;
        sigma = PerformGoldenBisection(sigma_up, sigma_lo, bisection_tol_,
                                       //sigma = PerformGoldenBisectionLog(sigma_up, sigma_lo, bisection_tol_,
                                       *step_aff_x_L,
                                       *step_aff_x_U,
                                       *step_aff_s_L,
                                       *step_aff_s_U,
                                       *step_aff->y_c(),
                                       *step_aff->y_d(),
                                       *step_aff->z_L(),
                                       *step_aff->z_U(),
                                       *step_aff->v_L(),
                                       *step_aff->v_U(),
                                       *step_cen_x_L,
                                       *step_cen_x_U,
                                       *step_cen_s_L,
                                       *step_cen_s_U,
                                       *step_cen->y_c(),
                                       *step_cen->y_d(),
                                       *step_cen->z_L(),
                                       *step_cen->z_U(),
                                       *step_cen->v_L(),
                                       *step_cen->v_U());
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
                                 *step_aff->y_c(),
                                 *step_aff->y_d(),
                                 *step_aff->z_L(),
                                 *step_aff->z_U(),
                                 *step_aff->v_L(),
                                 *step_aff->v_U(),
                                 *step_cen_x_L,
                                 *step_cen_x_U,
                                 *step_cen_s_L,
                                 *step_cen_s_U,
                                 *step_cen->y_c(),
                                 *step_cen->y_d(),
                                 *step_cen->z_L(),
                                 *step_cen->z_U(),
                                 *step_cen->v_L(),
                                 *step_cen->v_U());
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
                                        *step_aff->y_c(),
                                        *step_aff->y_d(),
                                        *step_aff->z_L(),
                                        *step_aff->z_U(),
                                        *step_aff->v_L(),
                                        *step_aff->v_U(),
                                        *step_cen_x_L,
                                        *step_cen_x_U,
                                        *step_cen_s_L,
                                        *step_cen_s_U,
                                        *step_cen->y_c(),
                                        *step_cen->y_d(),
                                        *step_cen->z_L(),
                                        *step_cen->z_U(),
                                        *step_cen->v_L(),
                                        *step_cen->v_U());
      Index l_min = (Index)(-(log(avrg_compl)-log(1e-9))/log(base))-1;
      for (; l>=l_min; l--) {
        sigma = pow(base, l);
        Number q = CalculateQualityFunction(sigma,
                                            *step_aff_x_L,
                                            *step_aff_x_U,
                                            *step_aff_s_L,
                                            *step_aff_s_U,
                                            *step_aff->y_c(),
                                            *step_aff->y_d(),
                                            *step_aff->z_L(),
                                            *step_aff->z_U(),
                                            *step_aff->v_L(),
                                            *step_aff->v_U(),
                                            *step_cen_x_L,
                                            *step_cen_x_U,
                                            *step_cen_s_L,
                                            *step_cen_s_U,
                                            *step_cen->y_c(),
                                            *step_cen->y_d(),
                                            *step_cen->z_L(),
                                            *step_cen->z_U(),
                                            *step_cen->v_L(),
                                            *step_cen->v_U());
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

    // Store the affine search direction (in case it is needed in the
    // line search for a corrector step)
    IpData().set_delta_aff(step_aff);
    IpData().SetHaveAffineDeltas(true);

    // Now construct the overall search direction here
    SmartPtr<IteratesVector> step = IpData().curr()->MakeNewIteratesVector(true);
    step->AddTwoVectors(sigma, *step_cen, 1.0, *IpData().delta_aff(), 0.0);

    IpData().set_delta(step);
    IpData().SetHaveDeltas(true);

    ///////////////////////////////////////////////////////////////////////////
    // Release memory for temporary vectors used in CalculateQualityFunction //
    ///////////////////////////////////////////////////////////////////////////

    tmp_step_x_L_ = NULL;
    tmp_step_x_U_ = NULL;
    tmp_step_s_L_ = NULL;
    tmp_step_s_U_ = NULL;
    tmp_step_z_L_ = NULL;
    tmp_step_z_U_ = NULL;
    tmp_step_v_L_ = NULL;
    tmp_step_v_U_ = NULL;

    tmp_slack_x_L_ = NULL;
    tmp_slack_x_U_ = NULL;
    tmp_slack_s_L_ = NULL;
    tmp_slack_s_U_ = NULL;
    tmp_z_L_ = NULL;
    tmp_z_U_ = NULL;
    tmp_v_L_ = NULL;
    tmp_v_U_ = NULL;

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

  Number QualityFunctionMuOracle::CalculateQualityFunction
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
   const Vector& step_cen_v_U
  )
  {
    DBG_START_METH("QualityFunctionMuOracle::CalculateQualityFunction",
                   dbg_verbosity);

    Index n_dual = IpData().curr()->x()->Dim() + IpData().curr()->s()->Dim();
    Index n_pri = IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim();
    Index n_comp = IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim() +
                   IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim();

    tmp_step_x_L_->AddTwoVectors(1., step_aff_x_L, sigma, step_cen_x_L, 0.);
    tmp_step_x_U_->AddTwoVectors(1., step_aff_x_U, sigma, step_cen_x_U, 0.);
    tmp_step_s_L_->AddTwoVectors(1., step_aff_s_L, sigma, step_cen_s_L, 0.);
    tmp_step_s_U_->AddTwoVectors(1., step_aff_s_U, sigma, step_cen_s_U, 0.);
    tmp_step_z_L_->AddTwoVectors(1., step_aff_z_L, sigma, step_cen_z_L, 0.);
    tmp_step_z_U_->AddTwoVectors(1., step_aff_z_U, sigma, step_cen_z_U, 0.);
    tmp_step_v_L_->AddTwoVectors(1., step_aff_v_L, sigma, step_cen_v_L, 0.);
    tmp_step_v_U_->AddTwoVectors(1., step_aff_v_U, sigma, step_cen_v_U, 0.);

    // Compute the fraction-to-the-boundary step sizes
    Number tau = IpData().curr_tau();
    Number alpha_primal = IpCq().slack_frac_to_the_bound(tau,
                          *tmp_step_x_L_,
                          *tmp_step_x_U_,
                          *tmp_step_s_L_,
                          *tmp_step_s_U_);

    Number alpha_dual = IpCq().dual_frac_to_the_bound(tau,
                        *tmp_step_z_L_,
                        *tmp_step_z_U_,
                        *tmp_step_v_L_,
                        *tmp_step_v_U_);

    Number xi; // centrality measure

    tmp_slack_x_L_->AddTwoVectors(1., *IpCq().curr_slack_x_L(),
                                  alpha_primal, *tmp_step_x_L_, 0.);
    tmp_slack_x_U_->AddTwoVectors(1., *IpCq().curr_slack_x_U(),
                                  alpha_primal, *tmp_step_x_U_, 0.);
    tmp_slack_s_L_->AddTwoVectors(1., *IpCq().curr_slack_s_L(),
                                  alpha_primal, *tmp_step_s_L_, 0.);
    tmp_slack_s_U_->AddTwoVectors(1., *IpCq().curr_slack_s_U(),
                                  alpha_primal, *tmp_step_s_U_, 0.);

    tmp_z_L_->AddTwoVectors(1., *IpData().curr()->z_L(),
                            alpha_dual, *tmp_step_z_L_, 0.);
    tmp_z_U_->AddTwoVectors(1., *IpData().curr()->z_U(),
                            alpha_dual, *tmp_step_z_U_, 0.);
    tmp_v_L_->AddTwoVectors(1., *IpData().curr()->v_L(),
                            alpha_dual, *tmp_step_v_L_, 0.);
    tmp_v_U_->AddTwoVectors(1., *IpData().curr()->v_U(),
                            alpha_dual, *tmp_step_v_U_, 0.);

    tmp_slack_x_L_->ElementWiseMultiply(*tmp_z_L_);
    tmp_slack_x_U_->ElementWiseMultiply(*tmp_z_U_);
    tmp_slack_s_L_->ElementWiseMultiply(*tmp_v_L_);
    tmp_slack_s_U_->ElementWiseMultiply(*tmp_v_U_);

    DBG_PRINT_VECTOR(2, "compl_x_L", *tmp_slack_x_L_);
    DBG_PRINT_VECTOR(2, "compl_x_U", *tmp_slack_x_U_);
    DBG_PRINT_VECTOR(2, "compl_s_L", *tmp_slack_s_L_);
    DBG_PRINT_VECTOR(2, "compl_s_U", *tmp_slack_s_U_);

    xi = IpCq().CalcCentralityMeasure(*tmp_slack_x_L_, *tmp_slack_x_U_,
                                      *tmp_slack_s_L_, *tmp_slack_s_U_);

    Number dual_inf;
    Number primal_inf;
    Number compl_inf;

    switch (quality_function_norm_) {
      case NM_NORM_1:
      dual_inf = (1.-alpha_dual)*(IpCq().curr_grad_lag_x()->Asum() +
                                  IpCq().curr_grad_lag_s()->Asum());

      primal_inf = (1.-alpha_primal)*(IpCq().curr_c()->Asum() +
                                      IpCq().curr_d_minus_s()->Asum());

      compl_inf = tmp_slack_x_L_->Asum() + tmp_slack_x_U_->Asum() +
                  tmp_slack_s_L_->Asum() + tmp_slack_s_U_->Asum();

      dual_inf /= n_dual;
      if (n_pri>0) {
        primal_inf /= n_pri;
      }
      DBG_ASSERT(n_comp>0);
      compl_inf /= n_comp;
      break;
      case NM_NORM_2_SQUARED:
      dual_inf =
        pow(1.-alpha_dual, 2)*(pow(IpCq().curr_grad_lag_x()->Nrm2(), 2) +
                               pow(IpCq().curr_grad_lag_s()->Nrm2(), 2));
      primal_inf =
        pow(1.-alpha_primal, 2)*(pow(IpCq().curr_c()->Nrm2(), 2) +
                                 pow(IpCq().curr_d_minus_s()->Nrm2(), 2));
      compl_inf =
        pow(tmp_slack_x_L_->Nrm2(), 2) + pow(tmp_slack_x_U_->Nrm2(), 2) +
        pow(tmp_slack_s_L_->Nrm2(), 2) + pow(tmp_slack_s_U_->Nrm2(), 2);

      dual_inf /= n_dual;
      if (n_pri>0) {
        primal_inf /= n_pri;
      }
      DBG_ASSERT(n_comp>0);
      compl_inf /= n_comp;
      break;
      case NM_NORM_MAX:
      dual_inf =
        (1.-alpha_dual)*Max(IpCq().curr_grad_lag_x()->Amax(),
                            IpCq().curr_grad_lag_s()->Amax());
      primal_inf =
        (1.-alpha_primal)*Max(IpCq().curr_c()->Amax(),
                              IpCq().curr_d_minus_s()->Amax());
      compl_inf =
        Max(tmp_slack_x_L_->Amax(), tmp_slack_x_U_->Amax(),
            tmp_slack_s_L_->Amax(), tmp_slack_s_U_->Amax());
      break;
      case NM_NORM_2:
      dual_inf =
        (1.-alpha_dual)*sqrt(pow(IpCq().curr_grad_lag_x()->Nrm2(), 2) +
                             pow(IpCq().curr_grad_lag_s()->Nrm2(), 2));
      primal_inf =
        (1.-alpha_primal)*sqrt(pow(IpCq().curr_c()->Nrm2(), 2) +
                               pow(IpCq().curr_d_minus_s()->Nrm2(), 2));
      compl_inf =
        sqrt(pow(tmp_slack_x_L_->Nrm2(), 2) + pow(tmp_slack_x_U_->Nrm2(), 2) +
             pow(tmp_slack_s_L_->Nrm2(), 2) + pow(tmp_slack_s_U_->Nrm2(), 2));

      dual_inf /= sqrt((Number)n_dual);
      if (n_pri>0) {
        primal_inf /= sqrt((Number)n_pri);
      }
      DBG_ASSERT(n_comp>0);
      compl_inf /= sqrt((Number)n_comp);
      break;
      default:
      DBG_ASSERT(false && "Unknown value for quality_function_norm_");
    }

    Number quality_function = dual_inf + primal_inf + compl_inf;

    switch (quality_function_centrality_) {
      case CEN_NONE:
      //Nothing
      break;
      case CEN_LOG:
      quality_function -= compl_inf*log(xi);
      break;
      case CEN_RECIPROCAL:
      quality_function += compl_inf/xi;
      case CEN_CUBED_RECIPROCAL:
      quality_function += compl_inf/pow(xi,3);
      break;
      default:
      DBG_ASSERT("Unknown value for quality_function_centrality_");
    }

    switch (quality_function_balancing_term_) {
      case BT_NONE:
      //Nothing
      break;
      case BT_CUBIC:
      quality_function += pow(Max(0., Max(dual_inf,primal_inf)-compl_inf),3);
      break;
      default:
      DBG_ASSERT("Unknown value for quality_function_balancing term_");
    }

    Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                   "sigma = %8.2e d_inf = %18.12e p_inf = %18.12e cmpl = %18.12e q = %18.12e a_pri = %8.2e a_dual = %8.2e xi = %8.2e\n", sigma, dual_inf, primal_inf, compl_inf, quality_function, alpha_primal, alpha_dual, xi);


    return quality_function;
    //return compl_inf;
  }

  Number
  QualityFunctionMuOracle::PerformGoldenBisection
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
   const Vector& step_cen_v_U
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
                                            step_cen_v_U);
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
                                            step_cen_v_U);

    Index nbisections = 0;
    while ((sigma_up-sigma_lo)>=tol*sigma_up && nbisections<max_bisection_steps_) {
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
                                         step_cen_v_U);
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
                                         step_cen_v_U);
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
                                             step_cen_v_U);
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
                                             step_cen_v_U);
      if (qtmp < q) {
        sigma = sigma_lo;
        q = qtmp;
      }
    }

    return sigma;
  }

  /* AW: Tried search in the log space, but that was even worse than
     search in unscaled space */
  /*
  Number
  QualityFunctionMuOracle::PerformGoldenBisectionLog
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
   const Vector& step_cen_v_U
  )
  {
    Number log_sigma;
    Number log_sigma_up = log(sigma_up);
    Number log_sigma_lo = log(sigma_lo);

    Number log_sigma_up_in = log_sigma_up;
    Number log_sigma_lo_in = log_sigma_lo;
    Number gfac = (3.-sqrt(5.))/2.;
    Number log_sigma_mid1 = log_sigma_lo + gfac*(log_sigma_up-log_sigma_lo);
    Number log_sigma_mid2 = log_sigma_lo + (1.-gfac)*(log_sigma_up-log_sigma_lo);

    Number qmid1 = CalculateQualityFunction(exp(log_sigma_mid1),
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
                                            step_cen_v_U);
    Number qmid2 = CalculateQualityFunction(exp(log_sigma_mid2),
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
                                            step_cen_v_U);

    Index nbisections = 0;
    while ((log_sigma_up-log_sigma_lo)>=tol*log_sigma_up && nbisections<max_bisection_steps_) {
      nbisections++;
      if (qmid1 > qmid2) {
        log_sigma_lo = log_sigma_mid1;
        log_sigma_mid1 = log_sigma_mid2;
        qmid1 = qmid2;
        log_sigma_mid2 = log_sigma_lo + (1.-gfac)*(log_sigma_up-log_sigma_lo);
        qmid2 = CalculateQualityFunction(exp(log_sigma_mid2),
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
                                         step_cen_v_U);
      }
      else {
        log_sigma_up = log_sigma_mid2;
        log_sigma_mid2 = log_sigma_mid1;
        qmid2 = qmid1;
        log_sigma_mid1 = log_sigma_lo + gfac*(log_sigma_up-log_sigma_lo);
        qmid1 = CalculateQualityFunction(exp(log_sigma_mid1),
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
                                         step_cen_v_U);
      }
    }

    Number q;
    if (qmid1 < qmid2) {
      log_sigma = log_sigma_mid1;
      q = qmid1;
    }
    else {
      log_sigma = log_sigma_mid2;
      q = qmid2;
    }
    if (log_sigma_up == log_sigma_up_in) {
      Number qtmp = CalculateQualityFunction(exp(log_sigma_up),
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
                                             step_cen_v_U);
      if (qtmp < q) {
        log_sigma = log_sigma_up;
        q = qtmp;
      }
    }
    else if (log_sigma_lo == log_sigma_lo_in) {
      Number qtmp = CalculateQualityFunction(exp(log_sigma_lo),
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
                                             step_cen_v_U);
      if (qtmp < q) {
        log_sigma = log_sigma_lo;
        q = qtmp;
      }
    }

    return exp(log_sigma);
  }
  */


} // namespace Ipopt
