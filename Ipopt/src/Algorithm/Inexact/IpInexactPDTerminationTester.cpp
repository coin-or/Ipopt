// Copyright (C) 2008, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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
      1e-3,
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
      0, 1,
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

    std::string inexact_linear_system_scaling;
    options.GetStringValue("inexact_linear_system_scaling",
                           inexact_linear_system_scaling, prefix);
    if (inexact_linear_system_scaling=="slack-based") {
      requires_scaling_ = true;
    }
    else {
      requires_scaling_ = false;
    }

    return true;
  }

  bool InexactPDTerminationTester::InitializeSolve()
  {
    DBG_START_METH("InexactPDTerminationTester::InitializeSolve",
                   dbg_verbosity);

    bool compute_normal = InexData().compute_normal();

    // compute the current infeasibility
    c_norm_ = IpCq().curr_primal_infeasibility(NORM_2);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: c_norm = %23.16e\n", c_norm_);

    SmartPtr<const Vector> normal_x = InexData().normal_x();
    SmartPtr<const Vector> normal_s = InexData().normal_s();
    if (compute_normal) {
      // calculate Jacobian times normal step (no scaling relevant in this space)
      curr_Av_c_ = InexCq().curr_jac_times_normal_c();
      curr_Av_d_ = InexCq().curr_jac_times_normal_d();

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

      // compute norm of the normal step in the scaled space
      v_norm_scaled_ = InexCq().slack_scaled_norm(*normal_x, *normal_s);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: v_norm_scaled_ = %23.16e\n", v_norm_scaled_);

      // compute Wv (Hessian times normal step) in the unscaled space
      curr_Wv_x_ = InexCq().curr_W_times_vec_x(*normal_x);
      curr_Wv_s_ = InexCq().curr_W_times_vec_s(*normal_s);
    }
    else {
      curr_Av_c_ = NULL;
      curr_Av_d_ = NULL;
      c_plus_Av_norm_ = -1.;
      v_norm_scaled_ = -1.;
      curr_Wv_x_ = NULL;
      curr_Wv_s_ = NULL;
    }

    // store the previous gradient and Jacobian information
    SmartPtr<const Vector> last_grad_barrier_obj_x = curr_grad_barrier_obj_x_;
    SmartPtr<const Vector> last_grad_barrier_obj_s = curr_grad_barrier_obj_s_;
    SmartPtr<const Matrix> last_jac_c = curr_jac_c_;
    SmartPtr<const Matrix> last_jac_d = curr_jac_d_;
    SmartPtr<const Vector> last_scaling_slacks = curr_scaling_slacks_;
    last_Av_norm_ = curr_Av_norm_;

    // get the current gradient and Jacobian information
    curr_grad_barrier_obj_x_ = IpCq().curr_grad_barrier_obj_x();
    curr_grad_barrier_obj_s_ = IpCq().curr_grad_barrier_obj_s(); // (unscaled)
    curr_jac_c_ = IpCq().curr_jac_c();
    curr_jac_d_ = IpCq().curr_jac_d();
    curr_scaling_slacks_ = InexCq().curr_scaling_slacks();

    // calculate \nabla phi(x_{k}) + A_{k}^Ty_k (in scaled space)
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
    curr_nabla_phi_plus_ATy_s_->ElementWiseMultiply(*curr_scaling_slacks_);

    // calculate norms appearing in termination tests
    curr_tt2_norm_ = IpCq().CalcNormOfType(NORM_2, *curr_nabla_phi_plus_ATy_x_,
                                           *curr_nabla_phi_plus_ATy_s_);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: curr_tt2_norm_ = %23.16e\n", curr_tt2_norm_);
    if (compute_normal) {
      curr_Av_norm_ = IpCq().CalcNormOfType(NORM_2, *curr_Av_c_, *curr_Av_d_);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: curr_Av_norm_ = %23.16e\n", curr_Av_norm_);
      curr_tt1_norm_ = sqrt(pow(curr_tt2_norm_, 2) + pow(curr_Av_norm_, 2));
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT: curr_tt1_norm_ = %23.16e\n", curr_tt1_norm_);
    }
    else {
      curr_Av_norm_ = -1.;
      curr_tt1_norm_ = sqrt(pow(curr_tt2_norm_, 2) + pow(c_norm_, 2));
    }

    if (compute_normal && IsValid(last_grad_barrier_obj_x)) {
      // calculate \nabla phi(x_{k-1}) + A_{k-1}^Ty_k (in scaled space)
      SmartPtr<Vector> last_nabla_phi_plus_ATy_x =
        last_grad_barrier_obj_x->MakeNewCopy();
      last_jac_c->TransMultVector(1., *IpData().curr()->y_c(),
                                  1., *last_nabla_phi_plus_ATy_x);
      last_jac_d->TransMultVector(1., *IpData().curr()->y_d(),
                                  1., *last_nabla_phi_plus_ATy_x);
      SmartPtr<Vector> last_nabla_phi_plus_ATy_s =
        last_grad_barrier_obj_s->MakeNew();
      last_nabla_phi_plus_ATy_s->AddTwoVectors(1., *last_grad_barrier_obj_s,
          -1., *IpData().curr()->y_d(), 0.);
      last_nabla_phi_plus_ATy_s->ElementWiseMultiply(*last_scaling_slacks);
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

    // check if we need to test termination test 2
    Number ATc_norm = InexCq().curr_scaled_Ac_norm();
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: current ATc norm = %23.16e\n", ATc_norm);
    try_tt2_ = (ATc_norm <= tt_eps2_*curr_tt2_norm_);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: will%s try termination test 2.\n", try_tt2_ ? "" : " not");

    return true;
  }

  InexactPDTerminationTester::ETerminationTest
  InexactPDTerminationTester::
  TestTermination(Index ndim, const Number* sol, const Number* resid,
                  Index iter, Number norm2_rhs)
  {
    DBG_START_METH("InexactPDTerminationTester::TestTermination",
                   dbg_verbosity);

    bool compute_normal = InexData().compute_normal();

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

    DBG_PRINT_VECTOR(2, "sol_x", *sol_x);
    DBG_PRINT_VECTOR(2, "sol_s", *sol_s);
    DBG_PRINT_VECTOR(2, "sol_c", *sol_c);
    DBG_PRINT_VECTOR(2, "sol_d", *sol_d);

    SmartPtr<const Vector> resid_x;
    SmartPtr<const Vector> resid_s;
    SmartPtr<const Vector> resid_c;
    SmartPtr<const Vector> resid_d;
    GetVectors(ndim, resid, resid_x, resid_s, resid_c, resid_d);

    if (requires_scaling_) {
      SmartPtr<const Vector> scaling_vec = curr_scaling_slacks_;
      SmartPtr<Vector> tmp = sol_s->MakeNewCopy();
      tmp->ElementWiseMultiply(*scaling_vec);
      sol_s = ConstPtr(tmp);
      tmp = resid_s->MakeNewCopy();
      tmp->ElementWiseDivide(*scaling_vec);
      resid_s = ConstPtr(tmp);
    }

    //// Set algorithm
    if (!compute_normal) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "RUNNING TERMINATION TESTS FOR INEXACT NEWTON (which means u=d)\n");
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "RUNNING TERMINATION TESTS FOR INEXACT NEWTON - TRUST REGION\n");
    }

    // Get the tangential step and its scaled norm
    SmartPtr<const Vector> tangential_x;
    SmartPtr<const Vector> tangential_s;
    if (compute_normal) {
      SmartPtr<const Vector> normal_x = InexData().normal_x();
      SmartPtr<Vector> tmp = sol_x->MakeNew();
      tmp->AddTwoVectors(1., *sol_x, -1, *normal_x, 0.);
      tangential_x = ConstPtr(tmp);
      SmartPtr<const Vector> normal_s = InexData().normal_s();
      tmp = sol_s->MakeNew();
      tmp->AddTwoVectors(1., *sol_s, -1, *normal_s, 0.);
      tangential_s = ConstPtr(tmp);
    }
    else {
      tangential_x = sol_x;
      tangential_s = sol_s;
    }
    Number u_norm_scaled =
      InexCq().slack_scaled_norm(*tangential_x, *tangential_s);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: u_norm_scaled = %23.16e\n", u_norm_scaled);

    // Compute u^TWu
    DBG_PRINT_VECTOR(2, "tangential_x", *tangential_x);
    DBG_PRINT_VECTOR(2, "tangential_s", *tangential_s);
    SmartPtr<const Vector> Wu_x = InexCq().curr_W_times_vec_x(*tangential_x);
    SmartPtr<const Vector> Wu_s = InexCq().curr_W_times_vec_s(*tangential_s);
    DBG_PRINT_VECTOR(2, "Wu_x", *Wu_x);
    DBG_PRINT_VECTOR(2, "Wu_s", *Wu_s);
    Number uWu = Wu_x->Dot(*tangential_x) + Wu_s->Dot(*tangential_s);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: uWu = %23.16e\n", uWu);

    // Compute norm of c + Ad
    SmartPtr<Vector> c_plus_Ad_c = IpCq().curr_c()->MakeNewCopy();
    curr_jac_c_->MultVector(1., *sol_x, 1., *c_plus_Ad_c);
    SmartPtr<const Vector> curr_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<Vector> c_plus_Ad_d = curr_d_minus_s->MakeNew();
    c_plus_Ad_d->AddTwoVectors(1., *curr_d_minus_s, -1., *sol_s, 0.);
    curr_jac_d_->MultVector(1., *sol_x, 1., *c_plus_Ad_d);
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

    // TODO: AW wants to discuss with Frank
    Number Upsilon = -1.;
    Number Nu = -1.;
    if (!compute_normal) {
#if 0
      // Compute Nu = ||A*u||^2/||A||^2
      SmartPtr<const Vector> curr_Au_c = IpCq().curr_jac_c_times_vec(*tangential_x);
      SmartPtr<Vector> curr_Au_d = sol_s->MakeNew();
      curr_Au_d->AddTwoVectors(1., *IpCq().curr_jac_d_times_vec(*tangential_x), -1., *tangential_s, 0.);
      Number Nu = IpCq().CalcNormOfType(NORM_2, *curr_Au_c, *curr_Au_d);
      Nu = pow(Nu,2);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: ||A*u||^2 = %23.16e\n", Nu);
      Number A_norm2 = InexCq().curr_scaled_A_norm2();
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: ||A||^2 = %23.16e\n", A_norm2);
#endif
      Nu = 0;//Nu/A_norm2;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: Nu = ||A*u||^2/||A||^2 = %23.16e\n", Nu);

      // Compute Upsilon = ||u||^2 - Nu
      Upsilon = u_norm_scaled*u_norm_scaled - Nu;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: Upsilon = ||u||^2 - ||A*u||^2/||A||^2 = %23.16e\n", Upsilon);
    }

    // Base value, something on the order of square root of machine epsilon; TODO: find a better base value
    Number BasVal = Max(IpData().curr()->x()->Amax(), IpData().curr()->s()->Amax());

    // Check tangential component condition, part 1
    Number lhs;
    Number rhs;
    if (!compute_normal) {
      lhs = Upsilon;
      rhs = pow(tcc_psi_,2)*Nu;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TCC1 testing Upsilon(=%23.16e) <= (tcc_psi_^2)*Nu(=%23.16e) --> ", lhs, rhs);
    }
    else {
      lhs = u_norm_scaled;
      rhs = tcc_psi_*v_norm_scaled_;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TCC1 testing u_norm_scaled(=%23.16e) <= tcc_psi_*v_norm_scaled(=%23.16e) --> ", lhs, rhs);
    }
    bool tcc1 = Compare_le(lhs, rhs, BasVal);
    if (tcc1) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
    }

    // Check tangential component condition, part 2a
    const Number mu = IpData().curr_mu();
    rhs = 0.5*uWu;
    if (!compute_normal) {
      lhs = tcc_theta_*pow(mu,tcc_theta_mu_exponent_)*Upsilon;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TCC2a testing 0.5*uWu(=%23.16e) >= tcc_theta_*pow(mu,tcc_theta_mu_exponent_)*Upsilon(=%23.16e) -->", rhs, lhs);
    }
    else {
      lhs = tcc_theta_*pow(mu,tcc_theta_mu_exponent_)*pow(u_norm_scaled, 2);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TCC2a testing 0.5*uWu(=%23.16e) >= tcc_theta_*pow(mu,tcc_theta_mu_exponent_)*u_norm^2(=%23.16e) -->", rhs, lhs);
    }
    bool tcc2a = Compare_le(lhs, rhs, BasVal);
    if (tcc2a) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
    }

    // Check tangential component condition, part 2b (only in MIPS)
    bool tcc = tcc1;
    if (!tcc && tcc2a) {
      if (!compute_normal) {
        tcc = tcc2a;
      }
      else {
        lhs = 0.5*uWu + curr_grad_barrier_obj_x_->Dot(*tangential_x) + curr_grad_barrier_obj_s_->Dot(*tangential_s) + curr_Wv_x_->Dot(*tangential_x) + curr_Wv_s_->Dot(*tangential_s);
        rhs = tcc_zeta_*v_norm_scaled_;
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "TCC2b testing (grad_barr^Tu + v^TWu + 0.5*uWu)(=%23.16e) <= tcc_zeta_*v_norm(=%23.16e) -->", lhs, rhs);
        tcc = Compare_le(lhs, rhs, BasVal);
        if (tcc) {
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
        }
        else {
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
        }
      }
    }

    // Check tangential component condition
    if (tcc) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Tangential Component Condition satisfied\n");
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Tangential Component Condition violated\n");
    }

    // Check termination test 1, residual condition
    bool tt1 = tcc;
    bool tt1_kappa1 = tcc;
    if (!compute_normal) {
      // Compute scaled norm of entire residual in case there is no step
      // decomposition.  In that case, c_plus_Ad_norm should indeed be
      // the same as what resid_c and resid_d woulod give (TODO:
      // check?!?)
      Number resid_norm = sqrt(pow(rho_norm, 2) + pow(c_plus_Ad_norm, 2));
      lhs = resid_norm;
      rhs = tt_kappa1_*curr_tt1_norm_;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT1 testing resid_norm(=%23.16e) <= tt_kappa1_*curr_tt1_norm_(=%23.16e) --> ", lhs, rhs);
      tt1 = Compare_le(lhs, rhs, BasVal);
    }
    else if (tt1) {
      lhs = rho_norm;
      rhs = tt_kappa1_*Min(curr_tt1_norm_, last_tt1_norm_);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TT1 testing rho_norm(=%23.16e) <= kappa1*min(curr_tt1_norm_, last_tt1_norm_)(=%23.16e) -->", lhs, rhs);
      tt1_kappa1 = Compare_le(lhs, rhs, BasVal);
      tt1 = tt1_kappa1;
    }
    if (tt1) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
    }

    // Check termination test 1, model reduction condition
    bool model_reduction = false;
    if (!compute_normal || tt1) {
      Number curr_nu = InexData().curr_nu();
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: curr_nu = %23.16e\n", curr_nu);
      Number delta_m = -(curr_grad_barrier_obj_x_->Dot(*sol_x) + curr_grad_barrier_obj_s_->Dot(*sol_s));
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: -grad_barr^Td = %23.16e\n", delta_m);
      delta_m += curr_nu*(c_norm_ - c_plus_Ad_norm);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: delta_m = %23.16e\n", delta_m);
      rhs = delta_m;
      Number sigma = rho_*tt_eps3_;
      if (!compute_normal) {
        lhs = Max(0.5*uWu, tcc_theta_*Upsilon) + sigma*curr_nu*Max(c_norm_, c_plus_Ad_norm - c_norm_);
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "MRC testing delta_m(=%23.16e) >= max(0.5*uWu,tcc_theta_*Upsilon) + sigma*nu*max(c_norm_, c_plus_Ad_norm - c_norm_)(=%23.16e) -->", rhs, lhs);
        model_reduction = Compare_le(lhs, rhs, BasVal);
        if (model_reduction) {
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
        }
        else {
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
        }
        if (tt1) {
          tt1 = model_reduction;
        }
      }
      else {
        lhs = Max(0.5*uWu, tcc_theta_*pow(u_norm_scaled, 2)) + sigma*curr_nu*(c_norm_ - c_plus_Av_norm_);
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "MRC testing delta_m(=%23.16e) >= max(0.5*uWu,tcc_theta_*u_norm^2) + sigma*nu*(c_norm_ - c_plus_Av_norm_)(=%23.16e) -->", rhs, lhs);
        tt1 = Compare_le(lhs, rhs, BasVal);
        if (tt1) {
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
        }
        else {
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
        }
      }
    }

    // Check termination test 1
    if (tt1) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 1 satisfied.\n");
      return TEST_1_SATISFIED;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 1 not satisfied.\n");
    }

    // Check termination test 3, residual condition
    bool tt3 = tcc;
    if (tt3) {
      if (!compute_normal) {
        lhs = rho_norm;
        rhs = tt_kappa1_*c_norm_;
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "TT3 testing rho_norm(=%23.16e) <= tt_kappa1_*c_norm_(=%23.16e) -->", lhs, rhs);
        tt3 = Compare_le(lhs, rhs, BasVal);
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "TT3 with residual condition from TT1 -->");
        tt3 = tt1_kappa1;
      }
      if (tt3) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
    }

    // Check termination test 3, linearized feasibility condition
    if (tt3) {
      if (!compute_normal) {
        lhs = c_plus_Ad_norm;
        rhs = tt_kappa1_*c_norm_;
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "TT3 testing c_plus_Ad_norm(=%23.16e) <= tt_kappa1_*c_norm(=%23.16e) -->", lhs, rhs);
        tt3 = Compare_le(lhs, rhs, BasVal);
      }
      else {
        lhs = tt_eps3_*(c_norm_ - c_plus_Av_norm_);
        rhs = c_norm_ - c_plus_Ad_norm;
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "TT3 testing (c_norm_-c_plus_Ad_norm)(=%23.16e) >= eps3*(c_norm_-c_plus_Av_norm_)(=%23.16e) -->", rhs, lhs);
        tt3 = Compare_le(lhs, rhs, BasVal);
      }
      if (tt3) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
    }

    // Check termination test 3
    if (tt3) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 3 satisfied.\n");
      return TEST_3_SATISFIED;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 3 not satisfied.\n");
    }

    // Check termination test 2
    bool tt2 = try_tt2_;
    if (tt2) {
      DBG_PRINT_VECTOR(2, "curr_nabla_phi_plus_ATy_x_", *curr_nabla_phi_plus_ATy_x_);
      DBG_PRINT_VECTOR(2, "curr_nabla_phi_plus_ATy_s_", *curr_nabla_phi_plus_ATy_s_);
      DBG_PRINT_VECTOR(2, "sol_c", *sol_c);
      DBG_PRINT_VECTOR(2, "sol_d", *sol_d);
      SmartPtr<Vector> sol_d_scaled = sol_d->MakeNewCopy();
      sol_d_scaled->ElementWiseMultiply(*curr_scaling_slacks_);
      DBG_PRINT_VECTOR(2, "sol_d_scaled", *sol_d_scaled);
      SmartPtr<Vector> nabla_phi_plus_ATydelta_x = curr_nabla_phi_plus_ATy_x_->MakeNewCopy();
      curr_jac_c_->TransMultVector(1., *sol_c, 1., *curr_nabla_phi_plus_ATy_x_);
      curr_jac_d_->TransMultVector(1., *sol_d, 1., *curr_nabla_phi_plus_ATy_x_);
      SmartPtr<Vector> nabla_phi_plus_ATydelta_s = curr_nabla_phi_plus_ATy_s_->MakeNew();
      nabla_phi_plus_ATydelta_s->AddTwoVectors(1., *curr_nabla_phi_plus_ATy_s_, -1., *sol_d_scaled, 0.);
      Number nabla_phi_plus_ATydelta_norm = IpCq().CalcNormOfType(NORM_2, *nabla_phi_plus_ATydelta_x, *nabla_phi_plus_ATydelta_s);
      lhs = nabla_phi_plus_ATydelta_norm;
      if (!compute_normal) {
        rhs = tt_kappa2_*curr_tt2_norm_;
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "TT2 testing ||gamma+A^T(y+delta)||(=%23.16e) <= kappa2*curr_tt2_norm_(=%23.16e) -->", lhs, rhs);
      }
      else {
        rhs = tt_kappa2_*Min(curr_tt2_norm_, last_tt1_norm_);
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "TT2 testing ||gamma+A^T(y+delta)||(=%23.16e) <= kappa2*min(curr_tt2_norm_, last_tt1_norm_)(=%23.16e) -->", lhs, rhs);
      }
      tt2 = Compare_le(lhs, rhs, BasVal);
      if (tt2) {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
      }
    }

    // Check termination test 2
    if (tt2) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 2 satisfied.\n");
      return TEST_2_SATISFIED;
    }
    else if (try_tt2_) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Termination Test 2 not satisfied.\n");
    }

    // Check if the Hessian should be modified
    if (tcc1 || tcc2a) {// || (!compute_normal && model_reduction)) {
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

    curr_Av_c_ = NULL;
    curr_Av_d_ = NULL;
    curr_Wv_x_ = NULL;
    curr_Wv_s_ = NULL;
    curr_nabla_phi_plus_ATy_x_ = NULL;
    curr_nabla_phi_plus_ATy_s_ = NULL;
  }

} // namespace Ipopt
