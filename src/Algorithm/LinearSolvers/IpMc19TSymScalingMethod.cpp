// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpoptConfig.h"
#include "IpMc19TSymScalingMethod.hpp"
#include "IpTypes.h"

#include <cmath>

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif

#if (defined(COINHSL_HAS_MC19) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MC19S) && defined(IPOPT_SINGLE))
#ifdef IPOPT_SINGLE
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name,NAME)
#else
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name ## d,NAME ## D)
#endif

/** Prototypes for MC19's Fortran subroutines */
extern "C"
{
   IPOPT_DECL_MC19A(IPOPT_HSL_FUNCP(mc19a, MC19A));
}
#else
#ifdef IPOPT_SINGLE
#define HSLFUNCNAMESUFFIX ""
#else
#define HSLFUNCNAMESUFFIX "d"
#endif
#endif

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

static IPOPT_DECL_MC19A(*user_mc19a) = NULL;

/// set MC19 function to use for every instantiation of this class
void Mc19TSymScalingMethod::SetFunctions(
   IPOPT_DECL_MC19A(*mc19a)
)
{
   DBG_ASSERT(mc19a != NULL);

   user_mc19a = mc19a;
}

IPOPT_DECL_MC19A(*Mc19TSymScalingMethod::GetMC19A())
{
   return user_mc19a;
}

bool Mc19TSymScalingMethod::InitializeImpl(
   const OptionsList& /*options*/,
   const std::string& /*prefix*/
)
{
   if( user_mc19a != NULL )
   {
      // someone set MC19 functions via setFunctions - prefer these
      mc19a = user_mc19a;
   }
   else
   {
#if (defined(COINHSL_HAS_MC19) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MC19S) && defined(IPOPT_SINGLE))
      // use HSL function that should be available in linked HSL library
      mc19a = &::IPOPT_HSL_FUNCP(mc19a, MC19A);
#else
      // try to load HSL function from a shared library at runtime
      DBG_ASSERT(IsValid(hslloader));

      mc19a = (IPOPT_DECL_MC19A(*))hslloader->loadSymbol("mc19a" HSLFUNCNAMESUFFIX);
#endif
   }

   DBG_ASSERT(mc19a != NULL);

   return true;
}

bool Mc19TSymScalingMethod::ComputeSymTScalingFactors(
   Index         n,
   Index         nnz,
   const Index*  airn,
   const Index*  ajcn,
   const Number* a,
   Number*       scaling_factors
)
{
   DBG_START_METH("Mc19TSymScalingMethod::ComputeSymTScalingFactors",
                  dbg_verbosity);

   if( DBG_VERBOSITY() >= 2 )
   {
      for( Index i = 0; i < nnz; i++ )
      {
         DBG_PRINT((2, "%5d A[%5d,%5d] = %23.15e\n", i, airn[i], ajcn[i], a[i]));
      }
   }
   // First copy the symmetric matrix into an unsymmetric (MA28)
   // format matrix
   Index* AIRN2 = new Index[2 * nnz];
   Index* AJCN2 = new Index[2 * nnz];
   Number* A2 = new Number[2 * nnz];
   Index nnz2 = 0;
   for( Index i = 0; i < nnz; i++ )
   {
      // ToDo decide if small values in A2 should be set to 0,
      // tolerance should probably be based on maximal element in A
      if( airn[i] == ajcn[i] )
      {
         AIRN2[nnz2] = airn[i];
         AJCN2[nnz2] = ajcn[i];
         A2[nnz2] = a[i];
         nnz2++;
      }
      else
      {
         AIRN2[nnz2] = airn[i];
         AJCN2[nnz2] = ajcn[i];
         A2[nnz2] = a[i];
         nnz2++;

         AIRN2[nnz2] = ajcn[i];
         AJCN2[nnz2] = airn[i];
         A2[nnz2] = a[i];
         nnz2++;
      }
   }

   if( DBG_VERBOSITY() >= 3 )
   {
      for( Index i = 0; i < nnz2; i++ )
      {
         DBG_PRINT((3, "A2[%5d] = %23.15e\n", i, A2[i]));
      }
   }

   // Call MC19 to get the scaling factors (for the matrix in the
   // general format)
   float* R = new float[n];
   float* C = new float[n];
   float* W = new float[5 * n];
   mc19a(&n, &nnz2, A2, AIRN2, AJCN2, R, C, W);
   delete[] W;

   if( DBG_VERBOSITY() >= 3 )
   {
      for( Index i = 0; i < n; i++ )
      {
         DBG_PRINT((3, "R[%5d] = %23.15e  C[%5d] = %23.15e\n",
                    i, R[i], i, C[i]));
      }
   }

   // Get the symmetric scaling factors as mean of the general ones
   // If some of the entries in A2 are too large, the scaling factors
   // are NaN.  Here, we check for this and return no scaling factors
   // if that is the case
   Number sum = 0.;
   Number smax = 0.;
   for( Index i = 0; i < n; i++ )
   {
      scaling_factors[i] = std::exp((Number) ((R[i] + C[i]) / 2.));
      sum += scaling_factors[i];
      smax = Max(smax, scaling_factors[i]);
   }
   if( !IsFiniteNumber(sum) || smax > 1e40 )
   {
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "Scaling factors are invalid - setting them all to 1.\n");
      for( Index i = 0; i < n; i++ )
      {
         scaling_factors[i] = 1.;
      }
   }

   if( DBG_VERBOSITY() >= 2 )
   {
      for( Index i = 0; i < n; i++ )
      {
         DBG_PRINT((2, "scaling_factors[%5d] = %23.15e\n",
                    i, scaling_factors[i]));
      }
   }

   // Clean up
   delete[] C;
   delete[] R;
   delete[] A2;
   delete[] AIRN2;
   delete[] AJCN2;

   return true;
}

} // namespace Ipopt
