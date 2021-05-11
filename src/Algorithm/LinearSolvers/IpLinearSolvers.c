// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include "IpLinearSolvers.h"
#include "IpoptConfig.h"
#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif

IpoptLinearSolver IpoptGetAvailableLinearSolvers(
   int buildinonly
)
{
   IpoptLinearSolver solvers = 0u;

#ifndef IPOPT_INT64
#if (defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA27S)) || (!defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA27))
   solvers |= IPOPTLINEARSOLVER_MA27;
#endif

#if (defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA57S)) || (!defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA57))
   solvers |= IPOPTLINEARSOLVER_MA57;
#endif

#if (defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA77S)) || (!defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA77))
   solvers |= IPOPTLINEARSOLVER_MA77;
#endif

#if (defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA86S)) || (!defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA86))
   solvers |= IPOPTLINEARSOLVER_MA86;
#endif

#if (defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA97S)) || (!defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MA97))
   solvers |= IPOPTLINEARSOLVER_MA97;
#endif

#if (defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MC19S)) || (!defined(IPOPT_SINGLE) && defined(COINHSL_HAS_MC19))
   solvers |= IPOPTLINEARSOLVER_MC19;
#endif

#ifdef PARDISO_LIB
   solvers |= IPOPTLINEARSOLVER_PARDISO;
#endif
#endif

#ifdef IPOPT_HAS_PARDISO_MKL
   solvers |= IPOPTLINEARSOLVER_PARDISOMKL;
#endif

#if !defined(IPOPT_SINGLE) && defined(IPOPT_HAS_SPRAL)
   solvers |= IPOPTLINEARSOLVER_SPRAL;
#endif

#if !defined(IPOPT_SINGLE) && defined(IPOPT_HAS_WSMP)
   solvers |= IPOPTLINEARSOLVER_WSMP;
#endif

#ifdef IPOPT_HAS_MUMPS
   solvers |= IPOPTLINEARSOLVER_MUMPS;
#endif

#if defined(IPOPT_HAS_LINEARSOLVERLOADER)
   if( !buildinonly )
   {
#ifndef IPOPT_INT64
      solvers |= IPOPTLINEARSOLVER_MA27;
      solvers |= IPOPTLINEARSOLVER_MA57;
      solvers |= IPOPTLINEARSOLVER_MA77;
      solvers |= IPOPTLINEARSOLVER_MA86;
      solvers |= IPOPTLINEARSOLVER_MA97;
      solvers |= IPOPTLINEARSOLVER_MC19;
      solvers |= IPOPTLINEARSOLVER_PARDISO;
#endif
   }
#else
   (void) buildinonly;
#endif

   return solvers;
}
