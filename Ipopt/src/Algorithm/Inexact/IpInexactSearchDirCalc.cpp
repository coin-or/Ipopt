// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-08-31

#include "IpInexactSearchDirCalc.hpp"

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  InexactSearchDirCalculator::
  InexactSearchDirCalculator(SmartPtr<InexactNormalStepCalculator> normal_step_calculator,
                             SmartPtr<InexactPDSolver> inexact_pd_solver)
      :
      normal_step_calculator_(normal_step_calculator),
      inexact_pd_solver_(inexact_pd_solver)
  {}

  InexactSearchDirCalculator::~InexactSearchDirCalculator()
  {}

  void InexactSearchDirCalculator::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "local_inf_Ac_tol",
      "Termination tolerance for local infeasibility (scaled ||Ac||).",
      0.0, true,
      1e-8,
      "");
  }

  bool InexactSearchDirCalculator::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("local_inf_Ac_tol", local_inf_Ac_tol_, prefix);

    bool retval = inexact_pd_solver_->Initialize(Jnlst(), IpNLP(), IpData(),
                  IpCq(), options, prefix);
    if (!retval) return false;
    return normal_step_calculator_->Initialize(Jnlst(), IpNLP(), IpData(),
           IpCq(), options, prefix);
  }

  bool
  InexactSearchDirCalculator::ComputeSearchDirection()
  {
    DBG_START_METH("InexactSearchDirCalculator::ComputeSearchDirection",
                   dbg_verbosity);

    // First check if the iterates have converged to a locally
    // infeasible point
    Number curr_scaled_Ac_norm = InexCq().curr_scaled_Ac_norm();
    Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                   "curr_scaled_Ac_norm = %e\n", curr_scaled_Ac_norm);
    Number curr_inf = IpCq().curr_primal_infeasibility(NORM_2);
    // ToDo work on termination criteria
    if (curr_scaled_Ac_norm <= local_inf_Ac_tol_ && curr_inf > 1e-4) {
      THROW_EXCEPTION(LOCALLY_INFEASIBLE,
                      "The scaled norm of Ac is satisfying tolerance");
    }

    SmartPtr<Vector> normal_x;
    SmartPtr<Vector> normal_s;

    bool retval =
      normal_step_calculator_->ComputeNormalStep(normal_x, normal_s);
    if (!retval) return false;

    // output
    if (Jnlst().ProduceOutput(J_VECTOR, J_SOLVE_PD_SYSTEM)) {
      Jnlst().Printf(J_VECTOR, J_SOLVE_PD_SYSTEM,
                     "Normal step (without slack scaling):\n");
      normal_x->Print(Jnlst(), J_VECTOR, J_SOLVE_PD_SYSTEM, "normal_x");
      normal_s->Print(Jnlst(), J_VECTOR, J_SOLVE_PD_SYSTEM, "normal_s");
    }

    // Now we set up the primal-dual system for computing the
    // tangential step and the search direction for the multipliers.
    // This is taken from IpPDSearchDirCal.cpp (rev 549).
    // We do not need entries for the variable bound multipliers

    SmartPtr<const IteratesVector> curr = IpData().curr();
    SmartPtr<IteratesVector> rhs = curr->MakeNewContainer();
    SmartPtr<Vector> tmp = curr->x()->MakeNew();
    tmp->AddOneVector(-1, *IpCq().curr_grad_lag_with_damping_x(), 0.);
    rhs->Set_x(*tmp);
    tmp = curr->s()->MakeNew();
    tmp->AddOneVector(-1., *IpCq().curr_grad_lag_with_damping_s(), 0.);
    rhs->Set_s(*tmp);
    tmp = curr->v_L()->MakeNew();
    tmp->AddOneVector(-1., *IpCq().curr_relaxed_compl_s_L(), 0.);
    rhs->Set_v_L(*tmp);
    tmp = curr->v_U()->MakeNew();
    tmp->AddOneVector(-1., *IpCq().curr_relaxed_compl_s_U(), 0.);
    rhs->Set_v_U(*tmp);

    // the right hand side for the constraints comes from the normal step
    rhs->Set_y_c(*IpCq().curr_jac_c_times_vec(*normal_x));
    tmp = normal_s->MakeNew();
    tmp->AddTwoVectors(1., *IpCq().curr_jac_d_times_vec(*normal_x),
                       -1., *normal_s, 0.);
    rhs->Set_y_d(*tmp);

    InexData().set_normal_x(normal_x);
    InexData().set_normal_s(normal_s);

    SmartPtr<IteratesVector> delta = rhs->MakeNewIteratesVector();
    retval = inexact_pd_solver_->Solve(*rhs, *delta);

    if (retval) {
      // Store the search directions in the IpData object
      IpData().set_delta(delta);
    }

    return retval;
  }

} // namespace Ipopt
