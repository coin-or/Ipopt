// Copyright (C) 2005, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpoptConfig.h"
#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif
#include "IpLinearSolversRegOp.hpp"
#include "IpRegOptions.hpp"
#include "IpTSymLinearSolver.hpp"

#include "IpMa27TSolverInterface.hpp"
#include "IpMa57TSolverInterface.hpp"
#include "IpMa77SolverInterface.hpp"
#include "IpMa86SolverInterface.hpp"
#include "IpMa97SolverInterface.hpp"
#include "IpMa28TDependencyDetector.hpp"
#include "IpPardisoSolverInterface.hpp"
#ifdef IPOPT_HAS_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif
#ifdef IPOPT_HAS_WSMP
# include "IpWsmpSolverInterface.hpp"
# include "IpIterativeWsmpSolverInterface.hpp"
#endif

namespace Ipopt
{

void RegisterOptions_LinearSolvers(
   const SmartPtr<RegisteredOptions>& roptions
)
{
   roptions->SetRegisteringCategory("Linear Solver");
   TSymLinearSolver::RegisterOptions(roptions);
#if defined(COINHSL_HAS_MA27) || defined(IPOPT_HAS_LINEARSOLVERLOADER)
   roptions->SetRegisteringCategory("MA27 Linear Solver");
   Ma27TSolverInterface::RegisterOptions(roptions);
#endif
#if defined(COINHSL_HAS_MA57) || defined(IPOPT_HAS_LINEARSOLVERLOADER)
   roptions->SetRegisteringCategory("MA57 Linear Solver");
   Ma57TSolverInterface::RegisterOptions(roptions);
#endif

#ifndef IPOPT_SINGLE // don't build these if using single precision
#if defined(COINHSL_HAS_MA77) || defined(IPOPT_HAS_LINEARSOLVERLOADER)
   roptions->SetRegisteringCategory("MA77 Linear Solver");
   Ma77SolverInterface::RegisterOptions(roptions);
#endif
#if defined(COINHSL_HAS_MA86) || defined(IPOPT_HAS_LINEARSOLVERLOADER)
   roptions->SetRegisteringCategory("MA86 Linear Solver");
   Ma86SolverInterface::RegisterOptions(roptions);
#endif
#if defined(COINHSL_HAS_MA97) || defined(IPOPT_HAS_LINEARSOLVERLOADER)
   roptions->SetRegisteringCategory("MA97 Linear Solver");
   Ma97SolverInterface::RegisterOptions(roptions);
#endif
#endif // IPOPT_SINGLE

#ifdef IPOPT_HAS_MUMPS
   roptions->SetRegisteringCategory("Mumps Linear Solver");
   MumpsSolverInterface::RegisterOptions(roptions);
#endif

#if defined(IPOPT_HAS_PARDISO) || defined(IPOPT_HAS_LINEARSOLVERLOADER)
   roptions->SetRegisteringCategory("Pardiso Linear Solver");
   PardisoSolverInterface::RegisterOptions(roptions);
#endif

#ifdef IPOPT_HAS_WSMP
   roptions->SetRegisteringCategory("WSMP Linear Solver");
   WsmpSolverInterface::RegisterOptions(roptions);
   IterativeWsmpSolverInterface::RegisterOptions(roptions);
#endif

#if defined(COINHSL_HAS_MA28) && defined(F77_FUNC)
   roptions->SetRegisteringCategory("MA28 Linear Solver");
   Ma28TDependencyDetector::RegisterOptions(roptions);
#endif

   roptions->SetRegisteringCategory("Uncategorized");
}

} // namespace Ipopt
