// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date    : 2009-04-06
// 
// Purpose : This is the same as IpSensitivityCalculator.hpp


#include "AsSimpleBacksolver.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  SimpleBacksolver::SimpleBacksolver(SmartPtr<PDSystemSolver> pd_solver)
    :
    pd_solver_(pd_solver)
  {
    DBG_START_METH("SimpleBacksolver::SimpleBacksolver", dbg_verbosity);
  }


  bool SimpleBacksolver::InitializeImpl(const OptionsList& options,
				   const std::string& prefix)
  {
    DBG_START_METH("SimpleBackSolver::InitializeImpl",dbg_verbosity);
    return true;
  }

  bool SimpleBacksolver::Solve(SmartPtr<IteratesVector> delta_lhs, SmartPtr<const IteratesVector> delta_rhs)
  {
    DBG_START_METH("SimpleBacksolver::Solve(IteratesVector,IteratesVector)", dbg_verbosity);
    bool retval;

    bool allow_inexact = false;
    bool improve_solution = false;

    retval = pd_solver_->Solve(1.0, 0.0, *delta_rhs, *delta_lhs, allow_inexact,
			       improve_solution);

    return retval;
  }
} // end namespace
