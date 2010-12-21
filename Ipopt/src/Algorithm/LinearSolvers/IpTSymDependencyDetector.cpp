// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2007-04-17

#include "IpoptConfig.h"
#include "IpTSymDependencyDetector.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  TSymDependencyDetector::
  TSymDependencyDetector(TSymLinearSolver& tsym_linear_solver)
      :
      tsym_linear_solver_(&tsym_linear_solver)
  {}

  void TSymDependencyDetector::
  RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    // Nothing so far
  }

  bool TSymDependencyDetector::InitializeImpl(
    const OptionsList& options,
    const std::string& prefix)
  {
    ASSERT_EXCEPTION(tsym_linear_solver_->ProvidesDegeneracyDetection(),
                     OPTION_INVALID,
                     "Selected linear solver does not support dependency detection");
    return tsym_linear_solver_->ReducedInitialize(Jnlst(), options, prefix);
  }

  bool TSymDependencyDetector::DetermineDependentRows(
    Index n_rows, Index n_cols, Index n_jac_nz, Number* jac_c_vals,
    Index* jac_c_iRow, Index* jac_c_jCol, std::list<Index>& c_deps)
  {
    DBG_START_METH("TSymDependencyDetector::DetermineDependentRows",
                   dbg_verbosity);

    ESymSolverStatus retval =
      tsym_linear_solver_->DetermineDependentRows(n_rows, n_cols, n_jac_nz,
          jac_c_vals, jac_c_iRow,
          jac_c_jCol, c_deps);

    return (retval == SYMSOLVER_SUCCESS);
  }

} // namespace Ipopt
