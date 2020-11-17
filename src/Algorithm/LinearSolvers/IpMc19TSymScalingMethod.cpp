// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpoptConfig.h"

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#else
/* if we build for the Linear Solver loader, then use normal C-naming style */
#define IPOPT_HSL_FUNC(name,NAME) name
#endif

// if we do not have MC19 in HSL or the linear solver loader, then we want to build the MC19 interface
#if defined(COINHSL_HAS_MC19) || defined(IPOPT_HAS_LINEARSOLVERLOADER)

#include "IpMc19TSymScalingMethod.hpp"

#include <cmath>

/** Prototypes for MC19's Fortran subroutines */
extern "C"
{
#ifdef IPOPT_SINGLE
   void IPOPT_HSL_FUNC(mc19a, MC19A)(
      ipfint* N,
      ipfint* NZ,
      float* A,
      ipfint* IRN,
      ipfint* ICN,
      float*  R,
      float*  C,
      float*  W
   );
#else
// here we assume that float corresponds to Fortran's single precision
   void IPOPT_HSL_FUNC(mc19ad, MC19AD)(
      ipfint* N,
      ipfint* NZ,
      double* A,
      ipfint* IRN,
      ipfint* ICN,
      float*  R,
      float*  C,
      float*  W
   );
#endif
}

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

bool Mc19TSymScalingMethod::InitializeImpl(
   const OptionsList& /*options*/,
   const std::string& /*prefix*/
)
{
   return true;
}

bool Mc19TSymScalingMethod::ComputeSymTScalingFactors(
   Index         n,
   Index         nnz,
   const ipfint* airn,
   const ipfint* ajcn,
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
   ipfint* AIRN2 = new ipfint[2 * nnz];
   ipfint* AJCN2 = new ipfint[2 * nnz];
   Number* A2 = new Number[2 * nnz];
   ipfint nnz2 = 0;
   for( Index i = 0; i < nnz; i++ )
   {
      if( airn[i] == ajcn[i] )
      {
         AIRN2[nnz2] = airn[i];
         AJCN2[nnz2] = ajcn[i];
         /*
          // ToDo decide if there should be a cut-off for small values in A
          // probably based on maximal element in A
          // DELETEME
          if (fabs(a[i])<1e-10) {
          A2[nnz2] = 0.;
          }
          else {
          A2[nnz2] = a[i];
          }
          */
         A2[nnz2] = a[i];
         nnz2++;
      }
      else
      {
         AIRN2[nnz2] = airn[i];
         AJCN2[nnz2] = ajcn[i];
         /*
          // DELETEME
          if (fabs(a[i])<1e-10) {
          A2[nnz2] = 0.;
          }
          else {
          A2[nnz2] = a[i];
          }
          */
         A2[nnz2] = a[i];
         nnz2++;
         AIRN2[nnz2] = ajcn[i];
         AJCN2[nnz2] = airn[i];
         /*
          // DELETEME
          if (fabs(a[i])<1e-10) {
          A2[nnz2] = 0.;
          }
          else {
          A2[nnz2] = a[i];
          }
          */
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
#ifdef IPOPT_SINGLE
   IPOPT_HSL_FUNC(mc19a, MC19A)(&n, &nnz2, A2, AIRN2, AJCN2, R, C, W);
#else
   IPOPT_HSL_FUNC(mc19ad, MC19AD)(&n, &nnz2, A2, AIRN2, AJCN2, R, C, W);
#endif
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
      scaling_factors[i] = exp((Number) ((R[i] + C[i]) / 2.));
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

#endif /* COINHSL_HAS_MC19 or IPOPT_HAS_LINEARSOLVERLOADER */
