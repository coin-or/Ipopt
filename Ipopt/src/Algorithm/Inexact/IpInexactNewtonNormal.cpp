// Copyright (C) 2008, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-08-31

#include "IpInexactNewtonNormal.hpp"
#include "IpSymLinearSolver.hpp"

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  InexactNewtonNormalStep::InexactNewtonNormalStep(SmartPtr<AugSystemSolver> aug_solver)
      :
      aug_solver_(aug_solver)
  {}

  InexactNewtonNormalStep::~InexactNewtonNormalStep()
  {}

  void InexactNewtonNormalStep::RegisterOptions(SmartPtr<RegisteredOptions> reg_options)
  {}

  bool InexactNewtonNormalStep::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return aug_solver_->Initialize(Jnlst(), IpNLP(), IpData(),
                                   IpCq(), options, prefix);
  }

  bool
  InexactNewtonNormalStep::ComputeNewtonNormalStep(Vector& newton_x,
      Vector& newton_s)

  {
    DBG_START_METH("InexactNewtonNormalStep::ComputeNormalNewtonStep",
                   dbg_verbosity);

    // Get the entires for the augmented system matrix

    // TODO: Make it possible to provide no Hessian!!!
    SmartPtr<const SymMatrix> zeroW = IpNLP().uninitialized_h();
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();
    // D_s is S^2 where S is the scaling factors from the slacks
    SmartPtr<const Vector> curr_scaling_slacks = InexCq().curr_scaling_slacks();
    SmartPtr<Vector> D_s = curr_scaling_slacks->MakeNewCopy();
    D_s->ElementWiseMultiply(*curr_scaling_slacks);
    D_s->ElementWiseReciprocal();

    // Get the entires for the right hand side
    SmartPtr<const Vector> curr_c = IpCq().curr_c();
    SmartPtr<const Vector> curr_d_minus_s = IpCq().curr_d_minus_s();
    SmartPtr<Vector> rhs_x = IpData().curr()->x()->MakeNew();
    rhs_x->Set(0.);
    SmartPtr<Vector> rhs_s = IpData().curr()->s()->MakeNew();
    rhs_s->Set(0.);

    // Get the space for the solution
    SmartPtr<Vector> sol_c = curr_c->MakeNew();
    SmartPtr<Vector> sol_d = curr_d_minus_s->MakeNew();

    ESymSolverStatus retval =
      aug_solver_->Solve(GetRawPtr(zeroW), 0., NULL, 1., GetRawPtr(D_s), 0.,
                         GetRawPtr(J_c), NULL, 0., GetRawPtr(J_d), NULL, 0.,
                         *rhs_x, *rhs_s, *curr_c, *curr_d_minus_s,
                         newton_x, newton_s, *sol_c, *sol_d, false, 0);

    if (retval==SYMSOLVER_SINGULAR) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                     "Resolving Newton step system with c-d perturbation.\n");
      retval = aug_solver_->Solve(GetRawPtr(zeroW), 0., NULL, 1., GetRawPtr(D_s), 0.,
                                  GetRawPtr(J_c), NULL, 1e-8, GetRawPtr(J_d), NULL, 1e-8,
                                  *rhs_x, *rhs_s, *curr_c, *curr_d_minus_s,
                                  newton_x, newton_s, *sol_c, *sol_d, false, 0);
    }

    if (retval!=SYMSOLVER_SUCCESS) return false;

    // return the step in the slack-scaled space
    newton_s.ElementWiseDivide(*curr_scaling_slacks);

    newton_x.Scal(-1.);
    newton_s.Scal(-1.);

    return true;
  }

} // namespace Ipopt
