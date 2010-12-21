// Copyright (C) 2008, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-09

#include "IpInexactPDSolver.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include "IpIterativeSolverTerminationTester.hpp"

extern Ipopt::IterativeSolverTerminationTester::ETerminationTest test_result_;

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  InexactPDSolver::InexactPDSolver(AugSystemSolver& augSysSolver,
                                   PDPerturbationHandler& perturbHandler)
      :
      augSysSolver_(&augSysSolver),
      perturbHandler_(&perturbHandler)
  {
    DBG_START_METH("InexactPDSolver::InexactPDSolver",dbg_verbosity);
  }

  InexactPDSolver::~InexactPDSolver()
  {
    DBG_START_METH("InexactPDSolver::~InexactPDSolver()",dbg_verbosity);
  }

  void InexactPDSolver::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddStringOption2(
      "modify_hessian_with_slacks",
      "Hessian modification strategy for slack part",
      "no",
      "no", "add multiple of identity",
      "yes", "add multiple of slacks squared inverse",
      "");
    roptions->AddLowerBoundedIntegerOption(
      "inexact_regularization_ls_count_trigger",
      "Threshold on line search count in previous iteration to trigger "
      "Hessian regularization.",
      1, 1,
      "If the ls count in the previous iteration is larger than this value, "
      "the Hessian will be regularized.");
  }


  bool InexactPDSolver::InitializeImpl(const OptionsList& options,
                                       const std::string& prefix)
  {
    options.GetNumericValue("tcc_psi", tcc_psi_, prefix);
    options.GetNumericValue("tcc_theta", tcc_theta_, prefix);
    options.GetNumericValue("tcc_theta_mu_exponent", tcc_theta_mu_exponent_, prefix);
    options.GetBoolValue("modify_hessian_with_slacks",
                         modify_hessian_with_slacks_, prefix);
    options.GetIntegerValue("inexact_regularization_ls_count_trigger",
                            inexact_regularization_ls_count_trigger_, prefix);

    std::string linear_solver;
    options.GetStringValue("linear_solver", linear_solver, prefix);
    is_pardiso_ = (linear_solver=="pardiso");

    if (!augSysSolver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                   options, prefix)) {
      return false;
    }

    last_info_ls_count_ = 0;

    return perturbHandler_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                       options, prefix);
  }

  bool InexactPDSolver::Solve(const IteratesVector& rhs,
                              IteratesVector& sol)
  {
    DBG_START_METH("InexactPDSolver::Solve",dbg_verbosity);

    // Timing of PDSystem solver starts here
    IpData().TimingStats().PDSystemSolverTotal().Start();

    DBG_PRINT_VECTOR(2, "rhs_x", *rhs.x());
    DBG_PRINT_VECTOR(2, "rhs_s", *rhs.s());
    DBG_PRINT_VECTOR(2, "rhs_c", *rhs.y_c());
    DBG_PRINT_VECTOR(2, "rhs_d", *rhs.y_d());
    DBG_PRINT_VECTOR(2, "rhs_vL", *rhs.v_L());
    DBG_PRINT_VECTOR(2, "rhs_vU", *rhs.v_U());
    DBG_PRINT_VECTOR(2, "sol_x in", *sol.x());
    DBG_PRINT_VECTOR(2, "sol_s in", *sol.s());
    DBG_PRINT_VECTOR(2, "sol_c in", *sol.y_c());
    DBG_PRINT_VECTOR(2, "sol_d in", *sol.y_d());
    DBG_PRINT_VECTOR(2, "sol_vL in", *sol.v_L());
    DBG_PRINT_VECTOR(2, "sol_vU in", *sol.v_U());

    // Receive data about matrix
    SmartPtr<const Vector> x = IpData().curr()->x();
    SmartPtr<const Vector> s = IpData().curr()->s();
    SmartPtr<const SymMatrix> W = IpData().W();
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<const Vector> v_L = IpData().curr()->v_L();
    SmartPtr<const Vector> v_U = IpData().curr()->v_U();
    SmartPtr<const Vector> slack_s_L = IpCq().curr_slack_s_L();
    SmartPtr<const Vector> slack_s_U = IpCq().curr_slack_s_U();
    SmartPtr<const Vector> sigma_s = IpCq().curr_sigma_s();
    DBG_PRINT_VECTOR(2, "Sigma_s", *sigma_s);

    // Compute the right hand side for the augmented system formulation
    SmartPtr<Vector> augRhs_s = rhs.s()->MakeNewCopy();
    Pd_L->AddMSinvZ(1.0, *slack_s_L, *rhs.v_L(), *augRhs_s);
    Pd_U->AddMSinvZ(-1.0, *slack_s_U, *rhs.v_U(), *augRhs_s);

    bool notDone = true;

    // Get the very first perturbation values from the perturbation
    // Handler
    Number delta_x;
    Number delta_s;
    Number delta_c;
    Number delta_d;
    perturbHandler_->ConsiderNewSystem(delta_x, delta_s, delta_c, delta_d);

    if (IpData().info_ls_count() > inexact_regularization_ls_count_trigger_) {
      perturbHandler_->PerturbForWrongInertia(delta_x, delta_s,
                                              delta_c, delta_d);
      if (last_info_ls_count_ > inexact_regularization_ls_count_trigger_) {
        perturbHandler_->PerturbForWrongInertia(delta_x, delta_s,
                                                delta_c, delta_d);
      }
    }
    last_info_ls_count_ = IpData().info_ls_count();

    ESymSolverStatus retval;
    Index count = 0;

    SmartPtr<const Vector> normal_x = InexData().normal_x();
    SmartPtr<const Vector> normal_s = InexData().normal_s();

    while (notDone) {
      count++;

      // Set the perturbation values in the Data object
      IpData().setPDPert(delta_x, delta_s, delta_c, delta_d);

      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Doing solve with perturbation parameters: delta_x=%e delta_s=%e\n                         delta_c=%e delta_d=%e\n",
                     delta_x, delta_s, delta_c, delta_d);
      if (delta_s>0. && modify_hessian_with_slacks_) {
        SmartPtr<const Vector> curr_scaling_slacks =
          InexCq().curr_scaling_slacks();
        SmartPtr<Vector> shifted_slacks = curr_scaling_slacks->MakeNewCopy();
        shifted_slacks->ElementWiseMultiply(*curr_scaling_slacks);
        shifted_slacks->ElementWiseReciprocal();
        const Number curr_mu = IpData().curr_mu();
        shifted_slacks->AddOneVector(1., *sigma_s, curr_mu*delta_s);
        retval = augSysSolver_->Solve(GetRawPtr(W), 1.0, NULL, delta_x,
                                      GetRawPtr(shifted_slacks), 0., GetRawPtr(J_c), NULL,
                                      delta_c, GetRawPtr(J_d), NULL, delta_d,
                                      *rhs.x(), *augRhs_s, *rhs.y_c(), *rhs.y_d(),
                                      *sol.x_NonConst(), *sol.s_NonConst(),
                                      *sol.y_c_NonConst(), *sol.y_d_NonConst(),
                                      false, 0);
      }
      else {
        retval = augSysSolver_->Solve(GetRawPtr(W), 1.0, NULL, delta_x,
                                      GetRawPtr(sigma_s), delta_s, GetRawPtr(J_c), NULL,
                                      delta_c, GetRawPtr(J_d), NULL, delta_d,
                                      *rhs.x(), *augRhs_s, *rhs.y_c(), *rhs.y_d(),
                                      *sol.x_NonConst(), *sol.s_NonConst(),
                                      *sol.y_c_NonConst(), *sol.y_d_NonConst(),
                                      false, 0);
      }
      if (retval==SYMSOLVER_SINGULAR) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "System seems singular.\n");
        if (InexData().compute_normal()) {
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "  We are already using the decomposition, now perturb the system.\n");
          bool pert_return = perturbHandler_->PerturbForSingularity(delta_x, delta_s,
                             delta_c, delta_d);
          if (!pert_return) {
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                           "PerturbForWrongInertia can't be done for singular.\n");
            IpData().TimingStats().PDSystemSolverTotal().End();
            return false;
          }
        }
        else {
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "  Switch to using the decomposition.\n");
          InexData().set_next_compute_normal(true);
          IpData().Append_info_string("@");
        }
      }
      else if (retval==SYMSOLVER_WRONG_INERTIA) {
        bool pert_return = perturbHandler_->PerturbForWrongInertia(delta_x, delta_s,
                           delta_c, delta_d);
        if (!pert_return) {
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "PerturbForWrongInertia can't be done for Hessian modification.\n");
          IpData().TimingStats().PDSystemSolverTotal().End();
          return false;
        }
      }
      else if (retval==SYMSOLVER_SUCCESS) {

        SmartPtr<const Vector> tangential_x;
        SmartPtr<const Vector> tangential_s;

        if (InexData().compute_normal()) {
          // Compute the tangetial part of the step from the overall step
          SmartPtr<Vector> tmp = normal_x->MakeNew();
          tmp->AddTwoVectors(1., *sol.x(), -1., *normal_x, 0.);
          tangential_x = ConstPtr(tmp);
          tmp = normal_s->MakeNew();
          tmp->AddTwoVectors(1., *sol.s(), -1., *normal_s, 0.);
          tangential_s = ConstPtr(tmp);
          // output
          if (Jnlst().ProduceOutput(J_MOREVECTOR, J_SOLVE_PD_SYSTEM)) {
            Jnlst().Printf(J_MOREVECTOR, J_SOLVE_PD_SYSTEM,
                           "Trial tangential step (without slack scaling):\n");
            tangential_x->Print(Jnlst(), J_MOREVECTOR, J_SOLVE_PD_SYSTEM, "tangential_x");
            tangential_s->Print(Jnlst(), J_MOREVECTOR, J_SOLVE_PD_SYSTEM, "tangential_s");
          }
        }
        else {
          tangential_x = sol.x();
          tangential_s = sol.s();
        }
        InexData().set_tangential_x(tangential_x);
        InexData().set_tangential_s(tangential_s);

        if (!is_pardiso_) {
          // check if we need to modify the system
          bool modify_hessian = HessianRequiresChange();
          if (modify_hessian) {
            bool pert_return = perturbHandler_->PerturbForWrongInertia(delta_x, delta_s,
                               delta_c, delta_d);
            if (!pert_return) {
              Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                             "PerturbForWrongInertia can't be done for Hessian modification.\n");
              IpData().TimingStats().PDSystemSolverTotal().End();
              return false;
            }
            retval = SYMSOLVER_WRONG_INERTIA;
          }
        }
        else {
          char buf[32];
          Snprintf(buf, 31, " TT=%d ", test_result_);
          IpData().Append_info_string(buf);
          if (test_result_ == IterativeSolverTerminationTester::CONTINUE) {
            if (InexData().compute_normal()) {
              Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                             "Termination tester not satisfied!!! Pretend singular\n");
              bool pert_return = perturbHandler_->PerturbForSingularity(delta_x, delta_s,
                                 delta_c, delta_d);
              if (!pert_return) {
                Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                               "PerturbForWrongInertia can't be done for singular.\n");
                IpData().TimingStats().PDSystemSolverTotal().End();
                return false;
              }
            }
            else {
              InexData().set_next_compute_normal(true);
              IpData().Append_info_string("@");
            }
          }
        }
        if (retval==SYMSOLVER_SUCCESS) {
          notDone = false;
        }
      }
      else {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Bad return code from augmented system solver = %d.\n",
                       retval);
        IpData().TimingStats().PDSystemSolverTotal().End();
        return false;
      }

    } // while (notDone)

    // Some output
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of trial factorizations performed: %d\n",
                   count);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Final perturbation parameters: delta_x=%e delta_s=%e\n                         delta_c=%e delta_d=%e\n",
                   delta_x, delta_s, delta_c, delta_d);

    // TODO: FRANK, how should we handle the multiplier updates for the slack duals??
#if 1
    // Compute the remaining sol Vectors
    Pd_L->SinvBlrmZMTdBr(-1., *slack_s_L, *rhs.v_L(), *v_L, *sol.s(), *sol.v_L_NonConst());
    Pd_U->SinvBlrmZMTdBr(1., *slack_s_U, *rhs.v_U(), *v_U, *sol.s(), *sol.v_U_NonConst());
#else
    // Compute the step for v_l and v_u simply from the y_d step.
    // This means that we are not really solving the linear system,
    // but those steps are consistent even if delta_s is positive
    Pd_L->TransMultVector(-1., *sol.y_d(), 0., *sol.v_L_NonConst());
    Pd_U->TransMultVector(1., *sol.y_d(), 0., *sol.v_U_NonConst());
#endif

    // Get space for the residual
    SmartPtr<IteratesVector> resid = sol.MakeNewIteratesVector(true);

    ComputeResiduals(*W, *J_c, *J_d, *Pd_L, *Pd_U, *v_L, *v_U,
                     *slack_s_L, *slack_s_U, *sigma_s, rhs, sol, *resid);

    DBG_PRINT_VECTOR(2, "sol", sol);
    IpData().TimingStats().PDSystemSolverTotal().End();

    return true;
  }

  void InexactPDSolver::ComputeResiduals(
    const SymMatrix& W,
    const Matrix& J_c,
    const Matrix& J_d,
    const Matrix& Pd_L,
    const Matrix& Pd_U,
    const Vector& v_L,
    const Vector& v_U,
    const Vector& slack_s_L,
    const Vector& slack_s_U,
    const Vector& sigma_s,
    const IteratesVector& rhs,
    const IteratesVector& res,
    IteratesVector& resid)
  {
    DBG_START_METH("InexactPDSolver::ComputeResiduals", dbg_verbosity);

    DBG_PRINT_VECTOR(2, "res", res);
    IpData().TimingStats().ComputeResiduals().Start();

    // Get the current sizes of the perturbation factors
    Number delta_x;
    Number delta_s;
    Number delta_c;
    Number delta_d;
    perturbHandler_->CurrentPerturbation(delta_x, delta_s, delta_c, delta_d);

    SmartPtr<Vector> tmp;

    // x
    W.MultVector(1., *res.x(), 0., *resid.x_NonConst());
    J_c.TransMultVector(1., *res.y_c(), 1., *resid.x_NonConst());
    J_d.TransMultVector(1., *res.y_d(), 1., *resid.x_NonConst());
    resid.x_NonConst()->AddTwoVectors(delta_x, *res.x(), -1., *rhs.x(), 1.);

    // s
    Pd_U.MultVector(1., *res.v_U(), 0., *resid.s_NonConst());
    Pd_L.MultVector(-1., *res.v_L(), 1., *resid.s_NonConst());
    resid.s_NonConst()->AddTwoVectors(-1., *res.y_d(), -1., *rhs.s(), 1.);
    if (delta_s!=0.) {
      resid.s_NonConst()->Axpy(delta_s, *res.s());
    }

    // c
    J_c.MultVector(1., *res.x(), 0., *resid.y_c_NonConst());
    resid.y_c_NonConst()->AddTwoVectors(-delta_c, *res.y_c(), -1., *rhs.y_c(), 1.);

    // d
    J_d.MultVector(1., *res.x(), 0., *resid.y_d_NonConst());
    resid.y_d_NonConst()->AddTwoVectors(-1., *res.s(), -1., *rhs.y_d(), 1.);
    if (delta_d!=0.) {
      resid.y_d_NonConst()->Axpy(-delta_d, *res.y_d());
    }

    // vL
    resid.v_L_NonConst()->Copy(*res.v_L());
    resid.v_L_NonConst()->ElementWiseMultiply(slack_s_L);
    tmp = v_L.MakeNew();
    Pd_L.TransMultVector(1., *res.s(), 0., *tmp);
    tmp->ElementWiseMultiply(v_L);
    resid.v_L_NonConst()->AddTwoVectors(1., *tmp, -1., *rhs.v_L(), 1.);

    // vU
    resid.v_U_NonConst()->Copy(*res.v_U());
    resid.v_U_NonConst()->ElementWiseMultiply(slack_s_U);
    tmp = v_U.MakeNew();
    Pd_U.TransMultVector(1., *res.s(), 0., *tmp);
    tmp->ElementWiseMultiply(v_U);
    resid.v_U_NonConst()->AddTwoVectors(-1., *tmp, -1., *rhs.v_U(), 1.);

    DBG_PRINT_VECTOR(2, "resid", resid);

    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_LINEAR_ALGEBRA)) {
      resid.Print(Jnlst(), J_MOREVECTOR, J_LINEAR_ALGEBRA, "resid");
    }

    if (Jnlst().ProduceOutput(J_MOREDETAILED, J_LINEAR_ALGEBRA)) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_x  %e\n", resid.x()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_s  %e\n", resid.s()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_c  %e\n", resid.y_c()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_d  %e\n", resid.y_d()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_vL %e\n", resid.v_L()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_vU %e\n", resid.v_U()->Amax());
    }
    IpData().TimingStats().ComputeResiduals().End();
  }

  bool
  InexactPDSolver::HessianRequiresChange()
  {
    // This code should be in sync with InexactPDTerminationTester
    bool compute_normal = InexData().compute_normal();

    SmartPtr<const Vector> normal_x = InexData().normal_x();
    SmartPtr<const Vector> normal_s = InexData().normal_s();
    SmartPtr<const Vector> tangential_x = InexData().tangential_x();
    SmartPtr<const Vector> tangential_s = InexData().tangential_s();

    Number u_norm_scaled =
      InexCq().slack_scaled_norm(*tangential_x, *tangential_s);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: u_norm_scaled = %23.16e\n", u_norm_scaled);

    Number Upsilon;
    Number Nu;
    Number v_norm_scaled = -1.;
    if (compute_normal) {
      v_norm_scaled = InexCq().slack_scaled_norm(*normal_x, *normal_s);
    }
    else {
      Nu = 0;//Nu/A_norm2;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: Nu = ||A*u||^2/||A||^2 = %23.16e\n", Nu);

      // Compute Upsilon = ||u||^2 - Nu
      Upsilon = u_norm_scaled - Nu;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "TT: Upsilon = ||u||^2 - ||A*u||^2/||A||^2 = %23.16e\n", Upsilon);
    }

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
      rhs = tcc_psi_*v_norm_scaled;
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "TCC1 testing u_norm_scaled(=%23.16e) <= tcc_psi_*v_norm_scaled(=%23.16e) --> ", lhs, rhs);
    }
    bool tcc1 = Compare_le(lhs, rhs, BasVal);
    if (tcc1) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "satisfied\n");
      return false;
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
    }

    // Compute u^TWu
    SmartPtr<const Vector> Wu_x = InexCq().curr_W_times_vec_x(*tangential_x);
    SmartPtr<const Vector> Wu_s = InexCq().curr_W_times_vec_s(*tangential_s);
    Number uWu = Wu_x->Dot(*tangential_x) + Wu_s->Dot(*tangential_s);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "TT: uWu = %23.16e\n", uWu);
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
      return false;
    }
    else {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA, "violated\n");
    }

    return true;
  }

} // namespace Ipopt
