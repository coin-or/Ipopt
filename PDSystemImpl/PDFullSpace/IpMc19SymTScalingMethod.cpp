// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpMc19SymTScalingMethod.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

/** Prototypes for MA27's Fortran subroutines */
extern "C"
{
  // here we assume that float corresponds to Fortran's single
  // precision
  void F77_FUNC(mc19ad,MM19AD)(ipfint *N, ipfint *NZ, double* A, ipfint *IRN,
                               ipfint* ICN, float* R, float* C, float* W);
}

namespace Ipopt
{

  static const Index dbg_verbosity = 0;

  bool Mc19SymTScalingMethod::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return true;
  }

  bool Mc19SymTScalingMethod::ComputeSymTScalingFactors(Index n,
      Index nnz,
      const ipfint* airn,
      const ipfint* ajcn,
      const double* a,
      double* scaling_factors)
  {
    DBG_START_METH("Mc19SymTScalingMethod::ComputeSymTScalingFactors",
                   dbg_verbosity);

    // First copy the symmetric matrix into an unsymmetric (MA28)
    // format matrix
    ipfint* AIRN2 = new ipfint[2*nnz];
    ipfint* AJCN2 = new ipfint[2*nnz];
    double* A2 = new double[2*nnz];
    ipfint nnz2=0;
    for (Index i=0; i<nnz; i++) {
      if (airn[i]==ajcn[i]) {
        AIRN2[nnz2] = airn[i];
        AJCN2[nnz2] = ajcn[i];
        A2[nnz2] = a[i];
        nnz2++;
      }
      else {
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

    if (DBG_VERBOSITY()>=3) {
      for (Index i=0; i<nnz2; i++) {
        DBG_PRINT((3, "A2[%5d] = %23.15e\n", i, A2[i]));
      }
    }

    // Call MC19 to get the scaling factors (for the matrix in the
    // general format)
    float* R = new float[n];
    float* C = new float[n];
    float* W = new float[5*n];
    F77_FUNC(mc19ad,MM19AD)(&n, &nnz2, A2, AIRN2, AJCN2, R, C, W);
    delete[] W;

    if (DBG_VERBOSITY()>=3) {
      for (Index i=0; i<n; i++) {
        DBG_PRINT((3, "R[%5d] = %23.15e  C[%5d] = %23.15e\n",
                   i, R[i], C[i]));
      }
    }

    // Get the symmetric scaling factors as mean of the general ones
    // If some of the entries in A2 are too large, the scaling factors
    // are NaN.  Here, we check for this and return no scaling factors
    // if that is the case
    Number sum=0.;
    for (Index i=0; i<n; i++) {
      scaling_factors[i] = exp((R[i]+C[i])/2.);
      sum += scaling_factors[i];
    }
    if (!FiniteNumber(sum)) {
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "Scaling factors are invalid - setting them all to 1.\n");
      for (Index i=0; i<n; i++) {
        scaling_factors[i] = 1.;
      }
    }

    if (DBG_VERBOSITY()>=2) {
      for (Index i=0; i<n; i++) {
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
