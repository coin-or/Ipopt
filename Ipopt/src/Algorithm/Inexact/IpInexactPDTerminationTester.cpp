// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-19

#include "IpInexactPDTerminationTester.hpp"
#include "IpBlas.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  InexactPDTerminationTester::InexactPDTerminationTester()
  {
    DBG_START_METH("InexactPDTerminationTester::InexactPDTerminationTester",dbg_verbosity);
  }

  InexactPDTerminationTester::~InexactPDTerminationTester()
  {
    DBG_START_METH("InexactPDTerminationTester::~InexactPDTerminationTester()",dbg_verbosity);
  }

  void InexactPDTerminationTester::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "tcc_psi",
      "Psi factor in Tangential Component Condition.",
      0.0, true,
      1e-1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "tcc_theta",
      "theta factor in Tangential Component Condition.",
      0.0, true,
      1e-12,
      "");
    roptions->AddLowerBoundedNumberOption(
      "tcc_theta_mu_exponent",
      "exponent for mu when multiplied with tcc_theta in Tangential Component Condition.",
      0., false,
      0.,
      "");
    roptions->AddLowerBoundedNumberOption(
      "tcc_zeta",
      "zeta factor in Tangential Component Condition.",
      0.0, true,
      1e-1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "tt_kappa1",
      "kappa1 factor in Termination Test 1 and 3.",
      0.0, true,
      1e-1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "tt_kappa2",
      "kappa2 factor in Termination Test 2.",
      0.0, true,
      1e-1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "tt_eps2",
      "eps2 factor in Termination Test 2.",
      0.0, true,
      1.,
      "");
    roptions->AddLowerBoundedNumberOption(
      "tt_eps3",
      "eps3 factor in Termination Test 3.",
      0.0, true,
      1.-1e-1,
      "");
    roptions->AddLowerBoundedNumberOption(
      "inexact_desired_pd_residual",
      "Desired relative residual tolerance for iterative solver during primal-dual step computation.",
      0.0, true, 1e-3,
      "");
    roptions->AddLowerBoundedIntegerOption(
      "inexact_desired_pd_residual_iter",
      "Number of iterations willing to be spent in obtaining desired primal-dual ration.",
      0, 100,
      "");
  }


  bool InexactPDTerminationTester::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("tcc_psi", tcc_psi_, prefix);
    options.GetNumericValue("tcc_theta", tcc_theta_, prefix);
    options.GetNumericValue("tcc_theta_mu_exponent", tcc_theta_mu_exponent_, prefix);
    options.GetNumericValue("tcc_zeta", tcc_zeta_, prefix);
    options.GetNumericValue("tt_kappa1", tt_kappa1_, prefix);
    options.GetNumericValue("tt_kappa2", tt_kappa2_, prefix);
    options.GetNumericValue("tt_eps2", tt_eps2_, prefix);
    options.GetNumericValue("tt_eps3", tt_eps3_, prefix);
    options.GetNumericValue("rho", rho_, prefix);
    options.GetNumericValue("inexact_desired_pd_residual",
                            inexact_desired_pd_residual_, prefix);
    options.GetIntegerValue("inexact_desired_pd_residual_iter",
                            inexact_desired_pd_residual_iter_, prefix);

    return true;
  }

  bool InexactPDTerminationTester::InitializeSolve()
  {
    DBG_START_METH("InexactPDTerminationTester::InitializeSolve",
                   dbg_verbosity);

    // calculate scaled Jacobian times normal step
    curr_Av_c_ = InexCq().curr_jac_times_normal_c();
    curr_Av_d_ = InexCq().curr_jac_times_normal_d();

    // compute the current infeasibility
    c_norm_ = IpCq().curr_primal_infeasibility(NORM_2);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: c_norm = %23.16e\n", c_norm_);

    // compute the linearized infeasibility at the normal step
    SmartPtr<const Vector> curr_c = IpCq().curr_c();
    SmartPtr<Vector> tmp1 = curr_c->MakeNew();
    tmp1->AddTwoVectors(1, *curr_c, 1., *curr_Av_c_, 0.);
    SmartPtr<const Vector> curr_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<Vector> tmp2 = curr_d_minus_s->MakeNew();
    tmp2->AddTwoVectors(1, *curr_d_minus_s, 1., *curr_Av_d_, 0.);
    c_plus_Av_norm_ = IpCq().CalcNormOfType(NORM_2, *tmp1, *tmp2);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: c_plus_Av_norm_ = %23.16e\n", c_plus_Av_norm_);
    curr_c_plus_Av_c_ = ConstPtr(tmp1);
    curr_c_plus_Av_d_ = ConstPtr(tmp2);

    // compute scaled norm of the normal step
    SmartPtr<const Vector> normal_x = InexData().normal_x();
    SmartPtr<const Vector> normal_s = InexData().normal_s();
    v_norm_scaled_ = InexCq().slack_scaled_norm(*normal_x, *normal_s);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: v_norm_scaled_ = %23.16e\n", v_norm_scaled_);

    // store the previous gradient and Jacobian information
    last_grad_barrier_obj_x_ = curr_grad_barrier_obj_x_;
    last_grad_barrier_obj_s_ = curr_grad_barrier_obj_s_;
    last_jac_c_ = curr_jac_c_;
    last_jac_d_ = curr_jac_d_;
    last_scaling_slacks_ = curr_scaling_slacks_;
    last_Av_norm_ = curr_Av_norm_;

    // get the current gradient and Jacobian information
    curr_grad_barrier_obj_x_ = IpCq().curr_grad_barrier_obj_x();
    curr_grad_barrier_obj_s_ = IpCq().curr_grad_barrier_obj_s();
    curr_jac_c_ = IpCq().curr_jac_c();
    curr_jac_d_ = IpCq().curr_jac_d();
    curr_scaling_slacks_ = InexCq().curr_scaling_slacks();

    // calculate \nabla phi(x_{k}) + A_{k}^Ty_k
    SmartPtr<const Vector> curr_jac_cT_times_curr_y_c =
      IpCq().curr_jac_cT_times_curr_y_c();
    SmartPtr<const Vector> curr_jac_cT_times_curr_y_d =
      IpCq().curr_jac_dT_times_curr_y_d();
    curr_nabla_phi_plus_ATy_x_ = curr_grad_barrier_obj_x_->MakeNewCopy();
    curr_nabla_phi_plus_ATy_x_->AddTwoVectors(1., *curr_jac_cT_times_curr_y_c,
        1., *curr_jac_cT_times_curr_y_d, 1.);
    curr_nabla_phi_plus_ATy_s_ = curr_grad_barrier_obj_s_->MakeNew();
    curr_nabla_phi_plus_ATy_s_->AddTwoVectors(1., *curr_grad_barrier_obj_s_,
        -1., *IpData().curr()->y_d(), 0.);
    // calculate norms appearing in termination tests
    curr_nabla_phi_plus_ATy_s_->ElementWiseMultiply(*curr_scaling_slacks_);
    curr_tt2_norm_ = IpCq().CalcNormOfType(NORM_2, *curr_nabla_phi_plus_ATy_x_,
                                           *curr_nabla_phi_plus_ATy_s_);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: curr_tt2_norm_ = %23.16e\n", curr_tt2_norm_);
    curr_Av_norm_ = IpCq().CalcNormOfType(NORM_2, *curr_Av_c_, *curr_Av_d_);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: curr_Av_norm_ = %23.16e\n", curr_Av_norm_);
    curr_tt1_norm_ = sqrt(pow(curr_tt2_norm_, 2) + pow(curr_Av_norm_, 2));
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: curr_tt1_norm_ = %23.16e\n", curr_tt1_norm_);

    if (IsValid(last_grad_barrier_obj_x_)) {
      // calculate \nabla phi(x_{k-1}) + A_{k-1}^Ty_k
      SmartPtr<Vector> last_nabla_phi_plus_ATy_x =
        last_grad_barrier_obj_x_->MakeNewCopy();
      last_jac_c_->TransMultVector(1., *IpData().curr()->y_c(),
                                   1., *last_nabla_phi_plus_ATy_x);
      last_jac_d_->TransMultVector(1., *IpData().curr()->y_d(),
                                   1., *last_nabla_phi_plus_ATy_x);
      SmartPtr<Vector> last_nabla_phi_plus_ATy_s =
        last_grad_barrier_obj_s_->MakeNew();
      last_nabla_phi_plus_ATy_s->AddTwoVectors(1., *last_grad_barrier_obj_s_,
          -1., *IpData().curr()->y_d(), 0.);
      last_nabla_phi_plus_ATy_s->ElementWiseMultiply(*last_scaling_slacks_);
      last_tt1_norm_ =
        IpCq().CalcNormOfType(NORM_2, *last_nabla_phi_plus_ATy_x,
                              *last_nabla_phi_plus_ATy_s);
      last_tt1_norm_ = sqrt(pow(last_tt1_norm_, 2) + pow(last_Av_norm_, 2));
    }
    else {
      last_tt1_norm_= 1e100;
    }
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: last_tt1_norm_ = %23.16e\n", last_tt1_norm_);

    // compute Wv (Hessian times normal step)
    curr_Wv_x_ = InexCq().curr_W_times_vec_x(*normal_x);
    curr_Wv_s_ = InexCq().curr_W_times_vec_s(*normal_s);

    // check if we need to test termination test 2
    Number ATc_norm = InexCq().curr_scaled_Ac_norm();
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: current ATc norm = %23.16e\n", ATc_norm);
    try_tt2_ = (ATc_norm <= tt_eps2_*curr_tt2_norm_);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: will %s try termination test 2.\n", try_tt2_ ? " not " : "");

    return true;
  }

  InexactPDTerminationTester::ETerminationTest
  InexactPDTerminationTester::
  TestTerminaion(Index ndim, const Number* sol, const Number* resid,
                 Index iter, Number norm2_rhs)
  {
    DBG_START_METH("InexactPDTerminationTester::TestTerminaion",
                   dbg_verbosity);

    last_iter_ = iter;

    ETerminationTest retval = CONTINUE;

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Starting PD Termination Tester for iteration %d.\n", iter);
    /*
    if (iter%5 != 4) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: immediately leaving tester for iteration %d.\n", iter);
      return retval;
    }
    */
    Number norm2_resid = IpBlasDnrm2(ndim, resid, 1);
    Number test_ratio = norm2_resid/norm2_rhs; // Min(norm2_resid/norm2_rhs, norm2_resid);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: test ratio %e (norm2_rhs = %e norm2_resid = %e).\n",
                   test_ratio, norm2_rhs, norm2_resid);
    if (iter < inexact_desired_pd_residual_iter_ &&
        test_ratio > inexact_desired_pd_residual_) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: immediately leaving tester with test ratio %e (norm2_rhs = %e norm2_resid = %e).\n",
                     test_ratio, norm2_rhs, norm2_resid);
      return retval;
    }

    SmartPtr<const Vector> sol_x;
    SmartPtr<const Vector> sol_s;
    SmartPtr<const Vector> sol_c;
    SmartPtr<const Vector> sol_d;
    GetVectors(ndim, sol, sol_x, sol_s, sol_c, sol_d);

    SmartPtr<const Vector> resid_x;
    SmartPtr<const Vector> resid_s;
    SmartPtr<const Vector> resid_c;
    SmartPtr<const Vector> resid_d;
    GetVectors(ndim, resid, resid_x, resid_s, resid_c, resid_d);

    // Get the tangential step and its scaled norm
    SmartPtr<const Vector> normal_x = InexData().normal_x();
    SmartPtr<Vector> tangential_x = sol_x->MakeNew();
    tangential_x->AddTwoVectors(1., *sol_x, -1, *normal_x, 0.);
    SmartPtr<const Vector> normal_s = InexData().normal_s();
    SmartPtr<Vector> tangential_s = sol_s->MakeNew();
    tangential_s->AddTwoVectors(1., *sol_s, -1, *normal_s, 0.);
    Number u_norm_scaled =
      InexCq().slack_scaled_norm(*tangential_x, *tangential_s);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: u_norm_scaled = %23.16e\n", u_norm_scaled);

    // Compute u^TWu
    SmartPtr<const Vector> Wu_x = InexCq().curr_W_times_vec_x(*tangential_x);
    SmartPtr<const Vector> Wu_s = InexCq().curr_W_times_vec_s(*tangential_s);
    Number uWu = Wu_x->Dot(*tangential_x) + Wu_s->Dot(*tangential_s);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: uWu = %23.16e\n", uWu);

    // Compute norm of c + Ad
    SmartPtr<Vector> c_plus_Ad_c = curr_c_plus_Av_c_->MakeNewCopy();
    curr_jac_c_->MultVector(1., *tangential_x, 1., *c_plus_Ad_c);
    SmartPtr<Vector> c_plus_Ad_d = curr_c_plus_Av_d_->MakeNew();
    c_plus_Ad_d->AddTwoVectors(1., *curr_c_plus_Av_d_, -1., *tangential_s, 0.);
    curr_jac_d_->MultVector(1., *tangential_x, 1., *c_plus_Ad_d);
    Number c_plus_Ad_norm =
      IpCq().CalcNormOfType(NORM_2, *c_plus_Ad_c, *c_plus_Ad_d);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: c_plus_Ad_norm  = %23.16e\n", c_plus_Ad_norm);

    // Compute norm of scaled residual rho
    SmartPtr<Vector> tmp = resid_s->MakeNewCopy();
    tmp->ElementWiseMultiply(*curr_scaling_slacks_);
    Number rho_norm = IpCq().CalcNormOfType(NORM_2, *resid_x, *tmp);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: rho_norm  = %23.16e\n", rho_norm);
    tmp = NULL;

    ////////// Check if the Tangential Component Condition is satisfied

    Number rhs = u_norm_scaled;
    Number lhs = tcc_psi_*v_norm_scaled_;
    // TODO: Find a good base value.  For now, we set this so that we
    // allow something on the order of the square root of the machine
    // precision...
    Number BasVal = Max(IpData().curr()->x()->Amax(), IpData().curr()->s()->Amax());
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TCC1 testing u_norm_scaled(=%23.16e) <= tcc_psi_*v_norm_scaled(=%23.16e) --> ", rhs, lhs);
    bool tcc1 = Compare_le(rhs, lhs, BasVal);

    if (tcc1) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
    }

    // check the first of the second pair of TCC (also needed in
    // Hessian update check)
    const Number mu = IpData().curr_mu();
    rhs = tcc_theta_*pow(mu,tcc_theta_mu_exponent_)*pow(u_norm_scaled, 2);
    lhs = 0.5*uWu;
    //const Number mach_eps_sqrt = pow(std::numeric_limits<Number>::epsilon(),0.5);
    //const Number mach_eps_sqrt = pow(std::numeric_limits<Number>::epsilon(),0.25);
    //BasVal = Max(IpData().curr()->x()->Amax(), IpData().curr()->s()->Amax())/mach_eps_sqrt;
    // check the second inequality of the tangential component condition
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TCC2a testing 0.5*uWu(=%23.16e) >= tcc_theta_*pow(mu,tcc_theta_mu_exponent_)*tangential_norm^2(=%23.16e) -->", lhs, rhs);
    bool tcc2a = Compare_le(rhs, lhs, BasVal);
    if (tcc2a) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
    }

    bool tcc = tcc1;
    if (!tcc && tcc2a) {
      // check the last condition
      rhs = 0.5*uWu + curr_grad_barrier_obj_x_->Dot(*tangential_x) +
            curr_grad_barrier_obj_s_->Dot(*tangential_s) +
            curr_Wv_x_->Dot(*tangential_x) + curr_Wv_s_->Dot(*tangential_s);
      lhs = tcc_zeta_*v_norm_scaled_;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TCC2b testing (grad_barr^Tu + v^TWu + 0.5*uWu)(=%23.16e) <= tcc_zeta_*v_norm(=%23.16e) -->", rhs, lhs);
      bool tcc2b = Compare_le(rhs, lhs, BasVal);
      if (tcc2b) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
      tcc = tcc2b;
    }
    if (tcc) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Tangential Component Condition satisfied\n");
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Tangential Component Condition violated\n");
    }

    ////////////  Termination Test 1
    bool tt1 = tcc;
    bool tt1_kappa1;
    if (tt1) {
      /////// Check residual condition for TT1
      rhs = rho_norm;
      lhs = tt_kappa1_*Min(curr_tt1_norm_, last_tt1_norm_);

      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT1 testing rho_norm(=%23.16e) <= kappa1*min(curr_tt1_norm_, last_tt1_norm_)(=%23.16e) -->", rhs, lhs);
      tt1_kappa1 = Compare_le(rhs, lhs, BasVal);
      tt1 = tt1_kappa1;
      if (tt1) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
    }

    if (tt1) {
      ///////  Check the model reduction condition
      Number curr_nu = InexData().curr_nu();
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: curr_nu = %23.16e\n", curr_nu);
      Number delta_m = -(curr_grad_barrier_obj_x_->Dot(*sol_x) +
                         curr_grad_barrier_obj_s_->Dot(*sol_s));
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: -grad_barr^Td = %23.16e\n", delta_m);
      delta_m += curr_nu*(c_norm_ - c_plus_Ad_norm);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: delta_m = %23.16e\n", delta_m);
      lhs = delta_m;

      Number sigma = rho_*tt_eps3_; // this is sigma in Model Reduction Cond.
      rhs = Max(0.5*uWu, tcc_theta_*pow(u_norm_scaled, 2)) +
            curr_nu*sigma*(c_norm_ - c_plus_Av_norm_);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "MRC testing delta_m(=%23.16e) >= max(0.5*uWu,tcc_theta_*u_norm^2) + curr_nu*sigma*(c_norm_ - c_plus_Av_norm_)(=%23.16e) -->", lhs, rhs);
      tt1 = Compare_le(rhs, lhs, BasVal);
      if (tt1) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
    }

    if (tt1) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 1 satisfied.\n");
      return TEST_1_SATISFIED;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 1 not satisfied.\n");
    }

    //////////////// Termination Test 3
    bool tt3 = tcc;
    if (tt3) {
      //////// Check residual condition for TT3
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT3 with residual condition from TT1 -->");
      tt3 = tt1_kappa1;
      if (tt3) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
    }
    if (tt3) {
      //////// Check linearized feasibility condition for TT3
      rhs = c_norm_-c_plus_Ad_norm;
      lhs = tt_eps3_*(c_norm_-c_plus_Av_norm_);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT3 testing (c_norm_-c_plus_Ad_norm)(=%23.16e) >= eps3(c_norm_-c_plus_Av_norm_)(=%23.16e) -->", rhs, lhs);
      tt3 = Compare_le(lhs, rhs, BasVal);
      if (tt3) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
    }

    if (tt3) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 3 satisfied.\n");
      return TEST_3_SATISFIED;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 3 not satisfied.\n");
    }

    //////////////// Termination Test 2
    if (try_tt2_) {
      DBG_PRINT_VECTOR(2, "curr_nabla_phi_plus_ATy_x_", *curr_nabla_phi_plus_ATy_x_);
      DBG_PRINT_VECTOR(2, "curr_nabla_phi_plus_ATy_s_", *curr_nabla_phi_plus_ATy_s_);
      DBG_PRINT_VECTOR(2, "sol_c", *sol_c);
      DBG_PRINT_VECTOR(2, "sol_d", *sol_d);
      SmartPtr<Vector> sol_d_scaled = sol_d->MakeNewCopy();
      sol_d_scaled->ElementWiseMultiply(*curr_scaling_slacks_);
      DBG_PRINT_VECTOR(2, "sol_d_scaled", *sol_d_scaled);
      SmartPtr<Vector> nabla_phi_plus_ATydelta_x =
        curr_nabla_phi_plus_ATy_x_->MakeNewCopy();
      curr_jac_c_->TransMultVector(1., *sol_c,
                                   1., *curr_nabla_phi_plus_ATy_x_);
      curr_jac_d_->TransMultVector(1., *sol_d,
                                   1., *curr_nabla_phi_plus_ATy_x_);
      SmartPtr<Vector> nabla_phi_plus_ATydelta_s =
        curr_nabla_phi_plus_ATy_s_->MakeNew();
      nabla_phi_plus_ATydelta_s->AddTwoVectors(1., *curr_nabla_phi_plus_ATy_s_,
          -1., *sol_d_scaled, 0.);
      Number nabla_phi_plus_ATydelta_norm =
        IpCq().CalcNormOfType(NORM_2, *nabla_phi_plus_ATydelta_x,
                              *nabla_phi_plus_ATydelta_s);
      /////// Check residual condition for TT1
      lhs = nabla_phi_plus_ATydelta_norm;
      rhs = tt_kappa2_*Min(curr_tt2_norm_, last_tt1_norm_);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT2 testing ||gamma+A^T(y+delta)||(=%23.16e) <= kappa2*min(curr_tt2_norm_, last_tt1_norm_)(=%23.16e) -->", lhs, rhs);
      bool tt2 = Compare_le(lhs, rhs, BasVal);
      if (tt2) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
      if (tt2) {
        return TEST_2_SATISFIED;
      }
    }

    // Check if the Hessian should be modified
    if (tcc1 || tcc2a) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Hessian Modification not requested.\n");
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Hessian Modification requested.\n");
      return MODIFY_HESSIAN;
    }

    return retval;
  }

  void
  InexactPDTerminationTester::Clear()
  {
    DBG_START_METH("InexactPDTerminationTester::Clear",
                   dbg_verbosity);

    last_grad_barrier_obj_x_ = NULL;
    last_grad_barrier_obj_s_ = NULL;
    last_jac_c_ = NULL;
    last_jac_d_ = NULL;
    last_scaling_slacks_ = NULL;

    curr_Av_c_ = NULL;
    curr_Av_d_ = NULL;
    curr_c_plus_Av_c_ = NULL;
    curr_c_plus_Av_d_ = NULL;
    curr_Wv_x_ = NULL;
    curr_Wv_s_ = NULL;
    curr_nabla_phi_plus_ATy_x_ = NULL;
    curr_nabla_phi_plus_ATy_s_ = NULL;
  }

} // namespace Ipopt
