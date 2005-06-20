// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpOptErrorConvCheck.hpp"

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  DefineIpoptType(OptimalityErrorConvergenceCheck);

  OptimalityErrorConvergenceCheck::OptimalityErrorConvergenceCheck()
  {}

  OptimalityErrorConvergenceCheck::~OptimalityErrorConvergenceCheck()
  {}

  void OptimalityErrorConvergenceCheck::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedIntegerOption("maxiter", "Maximum number of iterations to allow",
                                           0, 1000);
  }

  bool
  OptimalityErrorConvergenceCheck::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetIntegerValue("maxiter", max_iterations_, prefix);

    return true;
  }

  ConvergenceCheck::ConvergenceStatus OptimalityErrorConvergenceCheck::CheckConvergence()
  {
    DBG_START_METH("OptimalityErrorConvergenceCheck::CheckConvergence", dbg_verbosity);
    // maybe we should throw exceptions here instead?

    if (IpData().iter_count() >= max_iterations_) {
      return ConvergenceCheck::MAXITER_EXCEEDED;
    }

    Number overall_error = IpCq().curr_nlp_error();
    Number dual_inf = IpCq().curr_dual_infeasibility(NORM_MAX);
    Number primal_inf = IpCq().curr_primal_infeasibility(NORM_MAX);
    Number compl_inf = IpCq().curr_complementarity(0., NORM_MAX);
    DBG_PRINT((1,"overall_error = %8.2e\n",overall_error));
    if (overall_error <= IpData().tol() &&
        dual_inf <= IpData().dual_inf_tol() &&
        primal_inf <= IpData().primal_inf_tol() &&
        compl_inf <= IpData().compl_inf_tol()) {
      return ConvergenceCheck::CONVERGED;
    }

    return ConvergenceCheck::CONTINUE;
  }
} // namespace Ipopt
