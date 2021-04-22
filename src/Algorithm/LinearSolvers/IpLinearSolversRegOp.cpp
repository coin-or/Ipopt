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
   roptions->SetRegisteringCategory("Linear Solver", 360000);
   TSymLinearSolver::RegisterOptions(roptions);

   IpoptLinearSolver availablesolvers = IpoptGetAvailableLinearSolvers(false);

   if( availablesolvers & IPOPTLINEARSOLVER_MA27 )
   {
      roptions->SetRegisteringCategory("MA27 Linear Solver", 199000);
      Ma27TSolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA57 )
   {
      roptions->SetRegisteringCategory("MA57 Linear Solver", 198000);
      Ma57TSolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA77 )
   {
      roptions->SetRegisteringCategory("MA77 Linear Solver", 197000);
      Ma77SolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA86 )
   {
      roptions->SetRegisteringCategory("MA86 Linear Solver", 196000);
      Ma86SolverInterface::RegisterOptions(roptions);
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA97 )
   {
      roptions->SetRegisteringCategory("MA97 Linear Solver", 195000);
      Ma97SolverInterface::RegisterOptions(roptions);
   }

#ifdef IPOPT_HAS_MUMPS
   if( availablesolvers & IPOPTLINEARSOLVER_MUMPS )
   {
      roptions->SetRegisteringCategory("Mumps Linear Solver", 160000);
      MumpsSolverInterface::RegisterOptions(roptions);
   }
#endif

   if( availablesolvers & IPOPTLINEARSOLVER_PARDISO )
   {
      roptions->SetRegisteringCategory("Pardiso (pardiso-project.org) Linear Solver", 190000);
      PardisoSolverInterface::RegisterOptions(roptions);
   }

#ifdef IPOPT_HAS_PARDISO_MKL
   if( availablesolvers & IPOPTLINEARSOLVER_PARDISOMKL )
   {
      roptions->SetRegisteringCategory("Pardiso (MKL) Linear Solver", 189000);
      PardisoMKLSolverInterface::RegisterOptions(roptions);
   }
#endif

#ifdef IPOPT_HAS_SPRAL
   if( availablesolvers & IPOPTLINEARSOLVER_SPRAL )
   {
      roptions->SetRegisteringCategory("SPRAL Linear Solver", 180000);
      SpralSolverInterface::RegisterOptions(roptions);
   }
#endif

#ifdef IPOPT_HAS_WSMP
   if( availablesolvers & IPOPTLINEARSOLVER_WSMP )
   {
      roptions->SetRegisteringCategory("WSMP Linear Solver", 170000);
      WsmpSolverInterface::RegisterOptions(roptions);
      IterativeWsmpSolverInterface::RegisterOptions(roptions);
   }
#endif

#if ((defined(COINHSL_HAS_MA28) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA28S) && defined(IPOPT_SINGLE))) && defined(F77_FUNC)
   roptions->SetRegisteringCategory("MA28 Linear Solver", 150000);
   Ma28TDependencyDetector::RegisterOptions(roptions);
#endif
}

} // namespace Ipopt
