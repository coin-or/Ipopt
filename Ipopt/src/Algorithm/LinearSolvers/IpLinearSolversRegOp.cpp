// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpoptConfig.h"
#include "IpLinearSolversRegOp.hpp"
#include "IpRegOptions.hpp"
#include "IpTSymLinearSolver.hpp"

#include "IpMa27TSolverInterface.hpp"
#include "IpMa57TSolverInterface.hpp"
#include "IpMa28TDependencyDetector.hpp"
#include "IpPardisoSolverInterface.hpp"
#ifdef COIN_HAS_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif
#ifdef HAVE_WSMP
# include "IpWsmpSolverInterface.hpp"
#endif

namespace Ipopt
{

  void RegisterOptions_LinearSolvers(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Linear Solver");
    TSymLinearSolver::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("MA27 Linear Solver");
    Ma27TSolverInterface::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("MA57 Linear Solver");
    Ma57TSolverInterface::RegisterOptions(roptions);

#ifdef COIN_HAS_MUMPS

    roptions->SetRegisteringCategory("Mumps Linear Solver");
    MumpsSolverInterface::RegisterOptions(roptions);
#endif

    roptions->SetRegisteringCategory("Pardiso Linear Solver");
    PardisoSolverInterface::RegisterOptions(roptions);

#ifdef HAVE_WSMP

    roptions->SetRegisteringCategory("WSMP Linear Solver");
    WsmpSolverInterface::RegisterOptions(roptions);
#endif

    roptions->SetRegisteringCategory("MA28 Linear Solver");
    Ma28TDependencyDetector::RegisterOptions(roptions);

    roptions->SetRegisteringCategory("Uncategorized");
  }

} // namespace Ipopt
