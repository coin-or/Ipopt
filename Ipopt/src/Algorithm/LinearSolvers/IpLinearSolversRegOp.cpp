// Copyright (C) 2005, 2006 International Business Machines and others.
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

#ifdef HAVE_MA27
# include "IpMa27TSolverInterface.hpp"
#endif
#ifdef HAVE_MA57
# include "IpMa57TSolverInterface.hpp"
#endif
#ifdef COIN_HAS_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif
#ifdef HAVE_PARDISO
# include "IpPardisoSolverInterface.hpp"
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
#ifdef HAVE_MA27

    roptions->SetRegisteringCategory("MA27 Linear Solver");
    Ma27TSolverInterface::RegisterOptions(roptions);
#endif

#ifdef HAVE_MA57

    roptions->SetRegisteringCategory("MA57 Linear Solver");
    Ma57TSolverInterface::RegisterOptions(roptions);
#endif

#ifdef COIN_HAS_MUMPS

    roptions->SetRegisteringCategory("Mumps Linear Solver");
    MumpsSolverInterface::RegisterOptions(roptions);
#endif

#ifdef HAVE_PARDISO

    roptions->SetRegisteringCategory("Pardiso Linear Solver");
    PardisoSolverInterface::RegisterOptions(roptions);
#endif

#ifdef HAVE_WSMP

    roptions->SetRegisteringCategory("WSMP Linear Solver");
    WsmpSolverInterface::RegisterOptions(roptions);
#endif

    roptions->SetRegisteringCategory("Uncategorized");
  }

} // namespace Ipopt
