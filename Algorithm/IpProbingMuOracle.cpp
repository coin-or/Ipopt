// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpProbingMuOracle.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  static const Index dbg_verbosity = 0;

  ProbingMuOracle::ProbingMuOracle(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      MuOracle(),
      pd_solver_(pd_solver)
  {
    DBG_ASSERT(IsValid(pd_solver_));
  }

  ProbingMuOracle::~ProbingMuOracle()
  {}

  bool ProbingMuOracle::InitializeImpl(const OptionsList& options,
                                       const std::string& prefix)
  {
    // Check for the algorithm options
    Number value;
    if (options.GetNumericValue("sigma_max", value, prefix)) {
      ASSERT_EXCEPTION(value > 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"sigma_max\": This value must be positive.");
      sigma_max_ = value;
    }
    else {
      sigma_max_ = 1.;
    }

    // The following line is only here so that
    // IpoptCalculatedQuantities::CalculateSafeSlack and the first
    // output line have something to work with
    IpData().Set_mu(1.);

    return true;
  }

  Number ProbingMuOracle::CalculateMu()
  {
    DBG_START_METH("ProbingMuOracle::CalculateMu",
                   dbg_verbosity);

    /////////////////////////////////////
    // Compute the affine scaling step //
    /////////////////////////////////////

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the affine step\n");
    // First get the right hand side
    SmartPtr<const Vector> rhs_grad_lag_x  = IpCq().curr_grad_lag_x();
    SmartPtr<const Vector> rhs_grad_lag_s  = IpCq().curr_grad_lag_s();
    SmartPtr<const Vector> rhs_c = IpCq().curr_c();
    SmartPtr<const Vector> rhs_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<const Vector> rhs_compl_x_L = IpCq().curr_compl_x_L();
    SmartPtr<const Vector> rhs_compl_x_U = IpCq().curr_compl_x_U();
    SmartPtr<const Vector> rhs_compl_s_L = IpCq().curr_compl_s_L();
    SmartPtr<const Vector> rhs_compl_s_U = IpCq().curr_compl_s_U();

    // Get space for the affine scaling step
    SmartPtr<Vector> step_x = rhs_grad_lag_x->MakeNew();
    SmartPtr<Vector> step_s = rhs_grad_lag_s->MakeNew();
    SmartPtr<Vector> step_y_c = rhs_c->MakeNew();
    SmartPtr<Vector> step_y_d = rhs_d_minus_s->MakeNew();
    SmartPtr<Vector> step_z_L = rhs_compl_x_L->MakeNew();
    SmartPtr<Vector> step_z_U = rhs_compl_x_U->MakeNew();
    SmartPtr<Vector> step_v_L = rhs_compl_s_L->MakeNew();
    SmartPtr<Vector> step_v_U = rhs_compl_s_U->MakeNew();

    // Now solve the primal-dual system to get the step
    pd_solver_->Solve(-1.0, 0.0,
                      *rhs_grad_lag_x,
                      *rhs_grad_lag_s,
                      *rhs_c,
                      *rhs_d_minus_s,
                      *rhs_compl_x_L,
                      *rhs_compl_x_U,
                      *rhs_compl_s_L,
                      *rhs_compl_s_U,
                      *step_x,
                      *step_s,
                      *step_y_c,
                      *step_y_d,
                      *step_z_L,
                      *step_z_U,
                      *step_v_L,
                      *step_v_U,
                      true           // don't need high accuracy
                     );

    DBG_PRINT_VECTOR(2, "step_x", *step_x);
    DBG_PRINT_VECTOR(2, "step_s", *step_s);
    DBG_PRINT_VECTOR(2, "step_y_c", *step_y_c);
    DBG_PRINT_VECTOR(2, "step_y_d", *step_y_d);
    DBG_PRINT_VECTOR(2, "step_z_L", *step_z_L);
    DBG_PRINT_VECTOR(2, "step_z_U", *step_z_U);
    DBG_PRINT_VECTOR(2, "step_v_L", *step_v_L);
    DBG_PRINT_VECTOR(2, "step_v_U", *step_v_U);

    /////////////////////////////////////////////////////////////
    // Use Mehrotra's formula to compute the barrier parameter //
    /////////////////////////////////////////////////////////////

    // First compute the fraction-to-the-boundary step sizes
    Number alpha_primal_aff = IpCq().primal_frac_to_the_bound(1.0,
                              *step_x,
                              *step_s);

    Number alpha_dual_aff = IpCq().dual_frac_to_the_bound(1.0,
                            *step_z_L,
                            *step_z_U,
                            *step_v_L,
                            *step_v_U);

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  The affine maximal step sizes are\n"
                   "   alpha_primal_aff = %23.16e\n"
                   "   alpha_primal_aff = %23.16e\n",
                   alpha_primal_aff,
                   alpha_dual_aff);

    // now compute the average complementarity at the affine step
    Number mu_aff = CalculateAffineMu(alpha_primal_aff, alpha_dual_aff,
                                      *step_x, *step_s,
                                      *step_z_L, *step_z_U,
                                      *step_v_L, *step_v_U);
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  The average complementariy at the affine step is %23.16e\n",
                   mu_aff);

    // get the current average complementarity
    Number mu_curr = IpCq().curr_avrg_compl();
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  The average complementariy at the current point is %23.16e\n",
                   mu_curr);
    DBG_ASSERT(mu_curr>0.);

    // Apply Mehrotra's rule
    Number sigma = pow((mu_aff/mu_curr),3);
    // Make sure, sigma is not too large
    sigma = Min(sigma, sigma_max_);

    Number mu = sigma*mu_curr;

    // DELETEME
    char ssigma[40];
    sprintf(ssigma, "sigma=%e", sigma);
    IpData().Append_info_string(ssigma);
    sprintf(ssigma, "xi=%e", IpCq().curr_centrality_measure());
    IpData().Append_info_string(ssigma);

    return mu;
  }

  Number ProbingMuOracle::CalculateAffineMu
  (
    Number alpha_primal,
    Number alpha_dual,
    const Vector& step_x,
    const Vector& step_s,
    const Vector& step_z_L,
    const Vector& step_z_U,
    const Vector& step_v_L,
    const Vector& step_v_U)
  {
    // Get the current values of the slack variables and bound multipliers
    SmartPtr<const Vector> slack_x_L = IpCq().curr_slack_x_L();
    SmartPtr<const Vector> slack_x_U = IpCq().curr_slack_x_U();
    SmartPtr<const Vector> slack_s_L = IpCq().curr_slack_s_L();
    SmartPtr<const Vector> slack_s_U = IpCq().curr_slack_s_U();

    SmartPtr<const Vector> z_L = IpData().curr_z_L();
    SmartPtr<const Vector> z_U = IpData().curr_z_U();
    SmartPtr<const Vector> v_L = IpData().curr_v_L();
    SmartPtr<const Vector> v_U = IpData().curr_v_U();

    SmartPtr<Vector> tmp_slack;
    SmartPtr<Vector> tmp_mult;
    SmartPtr<const Matrix> P;
    Index ncomp = 0;
    Number sum =0.;

    // For each combination of slack and multiplier, compute the new
    // values and their dot products.

    // slack_x_L
    if (slack_x_L->Dim()>0) {
      ncomp += slack_x_L->Dim();

      P = IpNLP().Px_L();
      tmp_slack = slack_x_L->MakeNew();
      tmp_slack->Copy(*slack_x_L);
      P->TransMultVector(alpha_primal, step_x, 1.0, *tmp_slack);

      tmp_mult = z_L->MakeNew();
      tmp_mult->Copy(*z_L);
      tmp_mult->Axpy(alpha_dual, step_z_L);

      sum += tmp_slack->Dot(*tmp_mult);
    }

    // slack_x_U
    if (slack_x_U->Dim()>0) {
      ncomp += slack_x_U->Dim();

      P = IpNLP().Px_U();
      tmp_slack = slack_x_U->MakeNew();
      tmp_slack->Copy(*slack_x_U);
      P->TransMultVector(-alpha_primal, step_x, 1.0, *tmp_slack);

      tmp_mult = z_U->MakeNew();
      tmp_mult->Copy(*z_U);
      tmp_mult->Axpy(alpha_dual, step_z_U);

      sum += tmp_slack->Dot(*tmp_mult);
    }

    // slack_s_L
    if (slack_s_L->Dim()>0) {
      ncomp += slack_s_L->Dim();

      P = IpNLP().Pd_L();
      tmp_slack = slack_s_L->MakeNew();
      tmp_slack->Copy(*slack_s_L);
      P->TransMultVector(alpha_primal, step_s, 1.0, *tmp_slack);

      tmp_mult = v_L->MakeNew();
      tmp_mult->Copy(*v_L);
      tmp_mult->Axpy(alpha_dual, step_v_L);

      sum += tmp_slack->Dot(*tmp_mult);
    }

    // slack_s_U
    if (slack_s_U->Dim()>0) {
      ncomp += slack_s_U->Dim();

      P = IpNLP().Pd_U();
      tmp_slack = slack_s_U->MakeNew();
      tmp_slack->Copy(*slack_s_U);
      P->TransMultVector(-alpha_primal, step_s, 1.0, *tmp_slack);

      tmp_mult = v_U->MakeNew();
      tmp_mult->Copy(*v_U);
      tmp_mult->Axpy(alpha_dual, step_v_U);

      sum += tmp_slack->Dot(*tmp_mult);
    }

    DBG_ASSERT(ncomp>0);

    return sum/((Number)ncomp);
  }

} // namespace Ipopt
