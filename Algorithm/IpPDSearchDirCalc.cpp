// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter                 IBM    2005-10-13
//               derived from IpIpoptAlg.cpp

#include "IpPDSearchDirCalc.hpp"

namespace Ipopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  PDSearchDirCalculator::PDSearchDirCalculator(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      pd_solver_(pd_solver)
  {
    DBG_START_FUN("PDSearchDirCalculator::PDSearchDirCalculator",
                  dbg_verbosity);
    DBG_ASSERT(IsValid(pd_solver_));
  }

  PDSearchDirCalculator::~PDSearchDirCalculator()
  {
    DBG_START_FUN("PDSearchDirCalculator::~PDSearchDirCalculator()",
                  dbg_verbosity);
  }

  bool PDSearchDirCalculator::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return pd_solver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                  options, prefix);
  }

  void PDSearchDirCalculator::ComputeSearchDirection()
  {
    DBG_START_METH("PDSearchDirCalculator::ComputeSearchDirection",
                   dbg_verbosity);

    bool improve_solution = false;
    if (IpData().HaveDeltas()) {
      improve_solution = true;
    }

    SmartPtr<IteratesVector> rhs = IpData().curr()->MakeNewContainer();
    rhs->Set_x(*IpCq().curr_grad_lag_with_damping_x());
    rhs->Set_s(*IpCq().curr_grad_lag_with_damping_s());
    rhs->Set_y_c(*IpCq().curr_c());
    rhs->Set_y_d(*IpCq().curr_d_minus_s());
    rhs->Set_z_L(*IpCq().curr_relaxed_compl_x_L());
    rhs->Set_z_U(*IpCq().curr_relaxed_compl_x_U());
    rhs->Set_v_L(*IpCq().curr_relaxed_compl_s_L());
    rhs->Set_v_U(*IpCq().curr_relaxed_compl_s_U());

    DBG_PRINT_VECTOR(2, "rhs", *rhs);

    // Get space for the search direction
    SmartPtr<IteratesVector> delta =
      IpData().curr()->MakeNewIteratesVector(true);

    if (improve_solution) {
      // We can probably avoid copying and scaling...
      delta->AddOneVector(-1., *IpData().delta(), 0.);
    }

    bool allow_inexact = false;
    pd_solver_->Solve(-1.0, 0.0, *rhs, *delta, allow_inexact,
                      improve_solution);

    // Store the search directions in the IpData object
    IpData().set_delta(delta);
  }

} // namespace Ipopt
