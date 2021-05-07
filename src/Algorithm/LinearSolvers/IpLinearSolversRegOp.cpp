// Copyright (C) 2005, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpoptConfig.h"
#include "IpLinearSolversRegOp.hpp"
#include "IpLinearSolvers.h"
#include "IpRegOptions.hpp"
#include "IpTSymLinearSolver.hpp"

#include "IpMa27TSolverInterface.hpp"
#include "IpMa57TSolverInterface.hpp"
#include "IpMa77SolverInterface.hpp"
#include "IpMa86SolverInterface.hpp"
#include "IpMa97SolverInterface.hpp"
#include "IpMa28TDependencyDetector.hpp"
#include "IpPardisoSolverInterface.hpp"
#ifdef IPOPT_HAS_PARDISO_MKL
#include "IpPardisoMKLSolverInterface.hpp"
#endif
#ifdef IPOPT_HAS_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif
#ifdef IPOPT_HAS_SPRAL
# include "IpSpralSolverInterface.hpp"
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

   IpoptLinearSolver availablesolvers = IpoptGetAvailableLinearSolvers(false);

#ifndef IPOPT_INT64
   if( availablesolvers & IPOPTLINEARSOLVER_MA27 )
   {
      roptions->SetRegisteringCategory("MA27 Linear Solver");
      Ma27TSolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA57 )
   {
      roptions->SetRegisteringCategory("MA57 Linear Solver");
      Ma57TSolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA77 )
   {
      roptions->SetRegisteringCategory("MA77 Linear Solver");
      Ma77SolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA86 )
   {
      roptions->SetRegisteringCategory("MA86 Linear Solver");
      Ma86SolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA97 )
   {
      roptions->SetRegisteringCategory("MA97 Linear Solver");
      Ma97SolverInterface::RegisterOptions(roptions);
   }
#endif

#ifdef IPOPT_HAS_MUMPS
   if( availablesolvers & IPOPTLINEARSOLVER_MUMPS )
   {
      roptions->SetRegisteringCategory("Mumps Linear Solver");
      MumpsSolverInterface::RegisterOptions(roptions);
   }
#endif

#ifndef IPOPT_INT64
   if( availablesolvers & IPOPTLINEARSOLVER_PARDISO )
   {
      roptions->SetRegisteringCategory("Pardiso (pardiso-project.org) Linear Solver");
      PardisoSolverInterface::RegisterOptions(roptions);
   }
#endif

#ifdef IPOPT_HAS_PARDISO_MKL
   if( availablesolvers & IPOPTLINEARSOLVER_PARDISOMKL )
   {
      roptions->SetRegisteringCategory("Pardiso (MKL) Linear Solver");
      PardisoMKLSolverInterface::RegisterOptions(roptions);
   }
#endif

#ifdef IPOPT_HAS_SPRAL
   if( availablesolvers & IPOPTLINEARSOLVER_SPRAL )
   {
      roptions->SetRegisteringCategory("SPRAL Linear Solver");
      SpralSolverInterface::RegisterOptions(roptions);
   }
#endif

#ifdef IPOPT_HAS_WSMP
   if( availablesolvers & IPOPTLINEARSOLVER_WSMP )
   {
      roptions->SetRegisteringCategory("WSMP Linear Solver");
      WsmpSolverInterface::RegisterOptions(roptions);
      IterativeWsmpSolverInterface::RegisterOptions(roptions);
   }
#endif

#if ((defined(COINHSL_HAS_MA28) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA28S) && defined(IPOPT_SINGLE))) && defined(F77_FUNC)
   roptions->SetRegisteringCategory("MA28 Linear Solver");
   Ma28TDependencyDetector::RegisterOptions(roptions);
#endif
}

} // namespace Ipopt
